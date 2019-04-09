#[macro_use]
extern crate serde_json;
extern crate config;
extern crate libflate;
extern crate mediawiki;
extern crate reqwest;

use bio::io::{gaf, gff};
use config::{Config, File};
use libflate::gzip::Decoder;
use std::collections::HashMap;
use std::env;
use std::{error::Error, fmt};
use wikibase::entity_container::*;
use wikibase::entity_diff::*;
use wikibase::*;

const SPECIES_CONFIG_FILE: &str = "https://www.genedb.org/data/datasets.json";
const TMHMM_Q: &str = "Q61895944";

#[derive(Debug, Clone)]
pub struct GeneDBotError {}

// TODO use this!
impl Error for GeneDBotError {}

impl fmt::Display for GeneDBotError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "Oh no")
    }
}

#[derive(Debug, Clone)]
pub struct GeneDBotConfig {
    pub abbreviation: String,
    pub common_name: String,
    pub genus: String,
    pub species: String,
    pub strain: String,
    pub taxon_id: String,
    pub version: String,
    pub wikidata_id: String,
}

impl GeneDBotConfig {
    pub fn new_from_json(j: &serde_json::Value) -> Self {
        Self {
            abbreviation: j["abbreviation"].as_str().unwrap_or("").to_string(),
            common_name: j["common_name"].as_str().unwrap_or("").to_string(),
            genus: j["genus"].as_str().unwrap_or("").to_string(),
            species: j["species"].as_str().unwrap_or("").to_string(),
            strain: j["strain"].as_str().unwrap_or("").to_string(),
            taxon_id: j["taxon_id"].as_str().unwrap_or("").to_string(),
            version: j["version"].as_str().unwrap_or("").to_string(),
            wikidata_id: j["wikidata_id"].as_str().unwrap_or("").to_string(),
        }
    }
}

#[derive(Debug, Clone)]
pub struct GeneDBot {
    gff: HashMap<String, bio::io::gff::Record>,
    gaf: Vec<bio::io::gaf::Record>,
    config: GeneDBotConfig,
    species_key: String,
    genomic_assembly_q: String,
    ec: EntityContainer,
    api: mediawiki::api::Api,
    chr2q: HashMap<String, String>,
    genedb2q: HashMap<String, String>,
    alternate_gene_subclasses: HashMap<String, String>,
}

impl GeneDBot {
    pub fn new() -> Self {
        Self {
            gff: HashMap::new(),
            gaf: vec![],
            config: GeneDBotConfig::new_from_json(&json!({})),
            species_key: "".to_string(),
            genomic_assembly_q: "".to_string(),
            ec: EntityContainer::new(),
            api: mediawiki::api::Api::new("https://www.wikidata.org/w/api.php").unwrap(),
            chr2q: HashMap::new(),
            genedb2q: HashMap::new(),
            alternate_gene_subclasses: vec![
                ("tRNA", "Q201448"),
                ("rRNA", "Q215980"),
                ("pseudogene", "Q277338"),
                ("snoRNA", "Q284416"),
                ("ncRNA", "Q427087"),
                ("snRNA", "Q284578"),
            ]
            .iter()
            .map(|x| (x.0.to_string(), x.1.to_string()))
            .collect(),
        }
    }

    pub fn species_q(&self) -> String {
        self.config.wikidata_id.to_owned()
    }

    pub fn load_config_file(&mut self, species_key: &str) -> Result<(), reqwest::Error> {
        self.set_species(species_key);
        let config: serde_json::Value = reqwest::get(SPECIES_CONFIG_FILE)?.json()?;
        config
            .as_object()
            .unwrap()
            .iter()
            .for_each(|(_group, array)| {
                array
                    .as_array()
                    .unwrap()
                    .iter()
                    .filter(|x| x["abbreviation"].as_str().unwrap() == species_key)
                    .for_each(|x| self.config = GeneDBotConfig::new_from_json(x))
            });
        Ok(())
    }

    fn set_species(&mut self, species_key: &str) {
        self.species_key = species_key.into();
    }

    pub fn get_text_from_url(&self, url: &str) -> Result<String, reqwest::Error> {
        Ok(reqwest::get(url)?.text()?)
    }

    pub fn get_json_from_url(&self, url: &str) -> Result<serde_json::Value, reqwest::Error> {
        Ok(reqwest::get(url)?.json()?)
    }

    pub fn gff_url(&self) -> String {
        format!(
            "ftp://ftp.sanger.ac.uk/pub/genedb/releases/latest/{}/{}.gff.gz",
            self.species_key, self.species_key
        )
    }

    pub fn gaf_url(&self) -> String {
        format!(
            "ftp://ftp.sanger.ac.uk/pub/genedb/releases/latest/{}/{}.gaf",
            self.species_key, self.species_key
        )
    }

    pub fn load_gff_file_from_url(&mut self, url: &str) -> Result<(), Box<Error>> {
        let mut res = reqwest::get(url)?;
        let decoder = Decoder::new(&mut res)?;
        let mut reader = gff::Reader::new(decoder, gff::GffType::GFF3);
        for element in reader.records() {
            match element {
                Ok(e) => {
                    if e.attributes().contains_key("ID") {
                        let id = e.attributes()["ID"].clone();
                        self.gff.insert(id, e);
                    }
                }
                _ => continue,
            }
        }
        //println!("{:?}", &self.gff);
        Ok(())
    }

    pub fn load_gaf_file_from_url(&mut self, url: &str) -> Result<(), Box<Error>> {
        let res = reqwest::get(url)?;
        let mut reader = gaf::Reader::new(res, gaf::GafType::GAF2);
        for element in reader.records() {
            match element {
                Ok(e) => {
                    self.gaf.push(e);
                }
                _ => continue,
            }
        }
        //println!("{:?}", &self.gaf);
        Ok(())
    }

    pub fn load_gff_file(&mut self) -> Result<(), Box<Error>> {
        self.load_gff_file_from_url(self.gff_url().as_str())
    }

    pub fn load_gaf_file(&mut self) -> Result<(), Box<Error>> {
        self.load_gaf_file_from_url(self.gaf_url().as_str())
    }

    fn references(&self) -> Vec<Reference> {
        // TODO
        vec![]
    }

    fn create_genomic_assembly(&mut self) -> Result<String, Box<Error>> {
        let species_q = self.species_q();
        let species_i = self.ec.load_entity(&self.api, species_q.clone())?;
        let err1 = Err(From::from(format!(
            "'{}' has no P225 string value",
            &species_q
        )));
        let taxon_name = match species_i.values_for_property("P225").first() {
            Some(v) => match v {
                Value::StringValue(s) => s.clone(),
                _ => return err1,
            },
            None => return err1,
        };

        let references = self.references();
        let qualifiers = vec![];
        let mut new_item = Entity::new_empty_item();
        new_item.set_label(LocaleString::new("en", &(taxon_name + " reference genome")));
        new_item.add_claim(Statement::new_normal(
            Snak::new_item("P279", "Q7307127"),
            qualifiers.clone(),
            references.clone(),
        ));
        new_item.add_claim(Statement::new_normal(
            Snak::new_item("P703", species_q.as_str()),
            qualifiers.clone(),
            references.clone(),
        ));

        let params = EntityDiffParams::all();
        let diff = EntityDiff::new(&Entity::new_empty_item(), &new_item, &params);
        let res =
            EntityDiff::apply_diff(&mut self.api, &diff, EditTarget::New("item".to_string()))?;
        Ok(EntityDiff::get_entity_id(&res).unwrap())
    }

    fn find_genomic_assembly(&mut self) -> Result<(), Box<Error>> {
        let sparql = format!(
            "SELECT ?q {{ ?q wdt:P279 wd:Q7307127 ; wdt:P703 wd:{} }}",
            self.species_q()
        );
        let res = self.api.sparql_query(&sparql)?;
        let candidates = self.api.entities_from_sparql_result(&res, "q");
        self.genomic_assembly_q = match candidates.len() {
            0 => self.create_genomic_assembly()?,
            1 => candidates[0].to_owned(),
            _ => {
                return Err(From::from(format!(
                    "More than one genetic assembly for {}",
                    &sparql
                )));
            }
        };
        Ok(())
    }

    fn load_basic_items_chr(&mut self) -> Result<(), Box<Error>> {
        let sparql = format!(
            "SELECT ?q {{ ?q wdt:P31 wd:Q37748 ; wdt:P703 wd:{} }}",
            &self.species_q()
        );
        let res = self.api.sparql_query(&sparql)?;
        let mut items = self.api.entities_from_sparql_result(&res, "q");
        items.push(self.species_q());
        items.push(TMHMM_Q.to_string());
        self.ec.load_entities(&self.api, &items)?;

        items
            .iter()
            .for_each(|q| match self.ec.get_entity(q.to_string()) {
                Some(i) => {
                    if i.has_target_entity("P31", "Q37748") {
                        match i.label_in_locale("en") {
                            Some(label) => {
                                self.chr2q.insert(label.to_string(), q.to_string());
                            }
                            None => {}
                        }
                    }
                }
                None => {}
            });
        Ok(())
    }

    fn parent_taxon_q(&self) -> Option<String> {
        let species_q = self.species_q();
        let species_i = self.ec.get_entity(species_q.clone())?;
        match species_i.values_for_property("P171").first() {
            Some(v) => match v {
                Value::Entity(entity) => Some(entity.id().to_string()),
                _ => None,
            },
            None => None,
        }
    }

    pub fn load_basic_items_genes(&mut self) -> Result<(), Box<Error>> {
        let mut species_list = vec![self.species_q()];
        match self.parent_taxon_q() {
            Some(parent_q) => species_list.push(parent_q),
            None => {}
        }
        let species_list = format!(" VALUES ?species {{ wd:{} }} ", species_list.join(" wd:"));

        let alternate_gene_subclasses_values: Vec<String> = self
            .alternate_gene_subclasses
            .iter()
            .map(|(_, v)| v.to_owned())
            .collect();
        let gene_p31 = format!(
            " VALUES ?gene_types {{ wd:Q7187 wd:{} }} ",
            alternate_gene_subclasses_values.join(" wd:")
        );

        let sparql = format!("SELECT DISTINCT ?q ?genedb {{ {} . {} . ?q wdt:P31 ?gene_types  ; wdt:P703 ?species ; wdt:P3382 ?genedb }}",&species_list,&gene_p31) ;
        let res = self.api.sparql_query(&sparql)?;
        self.genedb2q = res["results"]["bindings"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|b| b["q"]["value"].as_str().is_some())
            .filter(|b| b["genedb"]["value"].as_str().is_some())
            .map(|b| {
                let entity_url = b["q"]["value"].as_str().unwrap();
                let q = self
                    .api
                    .extract_entity_from_uri(entity_url)
                    .unwrap()
                    .to_string();
                let genedb_id = b["genedb"]["value"].as_str().unwrap().to_string();
                (genedb_id, q)
            })
            .collect();
        Ok(())
    }

    pub fn load_basic_items(&mut self) -> Result<(), Box<Error>> {
        self.load_basic_items_chr()?;
        self.load_basic_items_genes()?;
        Ok(())
    }

    pub fn init(&mut self) -> Result<(), Box<Error>> {
        //self.load_gff_file().expect("Can't load GFF file");
        //self.load_gaf_file().expect("Can't load GAF file");
        let species_q = self.species_q();
        let _species_i = self.ec.load_entity(&self.api, species_q)?;
        self.find_genomic_assembly()?;
        println!("Genomic assembly: {}", self.genomic_assembly_q);
        self.load_basic_items()?;
        Ok(())
    }
}

fn main() {
    let args: Vec<_> = env::args().collect();
    if args.len() < 2 {
        panic!("Argument (species key) required\n");
    }
    let species_key = &args[1];

    let mut settings = Config::default();
    settings.merge(File::with_name("test.ini")).unwrap();
    let lgname = settings.get_str("user.user").unwrap();
    let lgpass = settings.get_str("user.pass").unwrap();

    let mut bot = GeneDBot::new();
    bot.api.login(lgname, lgpass).unwrap();
    bot.load_config_file(species_key)
        .expect("Can't load config file");
    //println!("{:?}", &bot.config);
    bot.init().unwrap();
}
