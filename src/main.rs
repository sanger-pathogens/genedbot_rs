#[macro_use]
extern crate serde_json;
extern crate chrono;
extern crate config;
extern crate libflate;
extern crate mediawiki;
extern crate reqwest;

use bio::io::{gaf, gff};
use chrono::Local;
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
    protein_genedb2q: HashMap<String, String>,
    evidence_codes_labels: HashMap<String, String>,
    evidence_codes: HashMap<String, String>,
    alternate_gene_subclasses: HashMap<String, String>,
    pub specific_genes_only: Option<Vec<String>>,
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
            protein_genedb2q: HashMap::new(),
            evidence_codes_labels: HashMap::new(),
            evidence_codes: HashMap::new(),
            specific_genes_only: None,
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

    fn sparql_result_to_pairs(
        &self,
        j: &serde_json::Value,
        k1: &str,
        k2: &str,
    ) -> Vec<(String, String)> {
        j.as_array()
            .unwrap()
            .iter()
            .filter(|b| b[k1]["value"].as_str().is_some())
            .filter(|b| b[k2]["value"].as_str().is_some())
            .map(|b| {
                let v1 = b[k1]["value"].as_str().unwrap();
                let v1 = self
                    .api
                    .extract_entity_from_uri(v1)
                    .unwrap_or(v1.to_string())
                    .to_string();
                let v2 = b[k2]["value"].as_str().unwrap();
                let v2 = self
                    .api
                    .extract_entity_from_uri(v2)
                    .unwrap_or(v2.to_string())
                    .to_string();
                (v1, v2)
            })
            .collect()
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

        // Genes
        let sparql = format!("SELECT DISTINCT ?q ?genedb {{ {} . {} . ?q wdt:P31 ?gene_types ; wdt:P703 ?species ; wdt:P3382 ?genedb }}",&species_list,&gene_p31) ;
        let res = self.api.sparql_query(&sparql)?;
        self.genedb2q = self
            .sparql_result_to_pairs(&res["results"]["bindings"], "genedb", "q")
            .into_iter()
            .collect();

        // Proteins
        let sparql = format!("SELECT DISTINCT ?q ?genedb {{ {} . ?q wdt:P31 wd:Q8054 ; wdt:P703 ?species ; wdt:P3382 ?genedb }}",&species_list) ;
        let res = self.api.sparql_query(&sparql)?;
        self.protein_genedb2q = self
            .sparql_result_to_pairs(&res["results"]["bindings"], "genedb", "q")
            .into_iter()
            .collect();
        Ok(())
    }

    pub fn load_evidence_codes(&mut self) -> Result<(), Box<Error>> {
        let sparql = "SELECT DISTINCT ?q ?qLabel ?qAltLabel { ?q wdt:P31 wd:Q23173209 SERVICE wikibase:label { bd:serviceParam wikibase:language \"en\" } }" ;
        let res = self.api.sparql_query(&sparql)?;
        res["results"]["bindings"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|b| b["q"]["value"].as_str().is_some())
            .filter(|b| b["qLabel"]["value"].as_str().is_some())
            .for_each(|b| {
                let q = b["q"]["value"].as_str().unwrap();
                let q = self.api.extract_entity_from_uri(q).unwrap().to_string();
                let label = b["qLabel"]["value"].as_str().unwrap_or("").to_string();
                let alt_label = b["qAltLabel"]["value"].as_str().unwrap_or("").to_string();
                self.evidence_codes_labels
                    .insert(self.normalize_evidence_label(&alt_label), q.clone());
                self.evidence_codes.insert(label, q.clone());
            });
        Ok(())
    }

    fn normalize_evidence_label(&self, s: &String) -> String {
        s.trim().to_lowercase()
    }

    fn get_gene_entities_to_process(&self) -> Vec<String> {
        match &self.specific_genes_only {
            Some(_genes) => vec![],
            None => self
                .genedb2q
                .iter()
                .map(|(_, v)| v.to_owned())
                .chain(self.protein_genedb2q.iter().map(|(_, v)| v.to_owned()))
                .collect(),
        }
    }

    fn get_gene_ids_to_process(&self) -> Vec<String> {
        match &self.specific_genes_only {
            Some(genes) => genes.clone(),
            None => self.genedb2q.iter().map(|(q, _)| q.to_owned()).collect(),
        }
    }

    pub fn load_basic_items_entities(&mut self) -> Result<(), Box<Error>> {
        let mut items_to_load: Vec<String> = self.get_gene_entities_to_process().to_vec();
        self.evidence_codes
            .iter()
            .map(|(_, v)| v.to_owned())
            .for_each(|s| items_to_load.push(s));
        //println!("Loading {} items", items_to_load.len());
        self.ec.load_entities(&self.api, &items_to_load)?;
        Ok(())
    }

    pub fn load_basic_items(&mut self) -> Result<(), Box<Error>> {
        self.load_basic_items_chr()?;
        self.load_basic_items_genes()?;
        self.load_evidence_codes()?;
        self.load_basic_items_entities()?;
        Ok(())
    }

    fn get_entity_id_for_genedb_id(&self, id: &String) -> Option<String> {
        match self.genedb2q.get(id) {
            Some(q) => Some(q.to_string()),
            None => match self.protein_genedb2q.get(id) {
                Some(q) => Some(q.to_string()),
                None => None,
            },
        }
    }

    fn get_entity_for_genedb_id(&mut self, id: &String) -> Option<&wikibase::Entity> {
        match self.get_entity_id_for_genedb_id(&id.to_string()) {
            Some(q) => self.ec.load_entity(&self.api, q).ok(),
            None => None,
        }
    }

    fn get_or_create_chromosome_entity(&mut self, id: &str) -> Option<&String> {
        match self.chr2q.get(id) {
            Some(q) => Some(q),
            None => {
                println!("TODO create new item for chromosome '{}'", id);
                None
            }
        }
    }

    fn process_gene(&mut self, genedb_id: String) {
        let _gene_q = match self.get_entity_id_for_genedb_id(&genedb_id) {
            Some(q) => q,
            None => return,
        };
        let gff = match self.gff.get(&genedb_id) {
            Some(gff) => gff.clone(),
            None => return,
        };
        let gene_type = match gff.feature_type() {
            "gene" => ("gene", "Q7187"),
            "pseudogene" => ("pseudogene", "Q277338"),
            other => {
                println!("Gene {} has unknown type {}", &genedb_id, other);
                return;
            }
        };
        dbg!(&gff);
        let mut item = Entity::new_empty_item();
        let _item_to_diff = match self.get_entity_for_genedb_id(&genedb_id) {
            Some(i) => i.clone(),
            None => Entity::new_empty_item(),
        };
        let chr_q = self
            .get_or_create_chromosome_entity(gff.seqname())
            .unwrap()
            .to_owned();

        // Labels and aliases
        match gff.attributes().get("Name") {
            Some(name) => {
                item.set_label(LocaleString::new("en", name));
                item.add_alias(LocaleString::new("en", &self.fix_alias_name(&genedb_id)));
            }
            None => item.set_label(LocaleString::new("en", &genedb_id)),
        };
        vec!["previous_systematic_id", "synonym", "alias"]
            .iter()
            .for_each(|key| {
                match gff.attributes().get(&key.to_string()) {
                    Some(ids) => ids.split(',').for_each(|id| {
                        item.add_alias(LocaleString::new("en", &self.fix_alias_name(id)))
                    }),
                    None => {}
                };
            });

        // Statements
        let today = Local::now();
        let reference = Reference::new(vec![
            Snak::new_item("P248", "Q5531047"),
            Snak::new_time(
                "P813",
                &format!("{}", today.format("+%Y-%m-%dT00:00:00Z")),
                11,
            ),
        ]);
        let ga_quals = vec![
            Snak::new_item("P659", &self.genomic_assembly_q),
            Snak::new_item("P1057", &chr_q),
        ];
        /*
        let strand_q = match gff.strand() {
            "+" => "Q22809680",
            _ => "Q22809711",
        };
        */

        // Instance of
        item.add_claim(Statement::new_normal(
            Snak::new_item("P31", gene_type.1),
            vec![],
            vec![reference.clone()],
        ));

        // Found in:Species
        item.add_claim(Statement::new_normal(
            Snak::new_item("P703", &self.species_q()),
            vec![],
            vec![reference.clone()],
        ));

        // Chromosome
        item.add_claim(Statement::new_normal(
            Snak::new_item("P1057", &chr_q),
            vec![],
            vec![reference.clone()],
        ));

        /*
        $gene_i->addClaim ( $gene_i->newClaim('P2548',$gene_i->newItem($strand_q) , [$refs] , $ga_quals ) ) ; # Strand
        */

        // GeneDB ID
        item.add_claim(Statement::new_normal(
            Snak::new_string("P3382", &genedb_id),
            vec![],
            vec![reference.clone()],
        ));

        // Genomic start
        item.add_claim(Statement::new_normal(
            Snak::new_string("P644", &gff.start().to_string()),
            ga_quals.clone(),
            vec![reference.clone()],
        ));

        // Genomic end
        item.add_claim(Statement::new_normal(
            Snak::new_string("P645", &gff.end().to_string()),
            ga_quals.clone(),
            vec![reference.clone()],
        ));

        println!("{}", serde_json::to_string_pretty(&item).unwrap());
    }

    fn fix_alias_name(&self, name: &str) -> String {
        name.split(";").next().unwrap().to_string()
    }

    pub fn run(&mut self) -> Result<(), Box<Error>> {
        self.get_gene_ids_to_process()
            .iter()
            .for_each(|genedb_id| self.process_gene(genedb_id.to_string()));
        Ok(())
    }

    pub fn init(&mut self) -> Result<(), Box<Error>> {
        self.load_gff_file().expect("Can't load GFF file");
        self.load_gaf_file().expect("Can't load GAF file");
        let species_q = self.species_q();
        let _species_i = self.ec.load_entity(&self.api, species_q)?;
        self.find_genomic_assembly()?;
        self.load_basic_items()?;
        Ok(())
    }
}

fn main() {
    let args: Vec<_> = env::args().collect();
    if args.len() < 2 {
        println!("Argument (species key) required\n");
        return;
    }
    let species_key = &args[1];

    let mut settings = Config::default();
    settings.merge(File::with_name("test.ini")).unwrap();
    let lgname = settings.get_str("user.user").unwrap();
    let lgpass = settings.get_str("user.pass").unwrap();

    let mut bot = GeneDBot::new();
    if args.len() == 3 {
        bot.specific_genes_only = Some(vec![args[2].to_string()]);
    }
    bot.api.login(lgname, lgpass).unwrap();
    //println!("is bot: {}", bot.api.user().is_bot());
    bot.load_config_file(species_key)
        .expect("Can't load config file");
    //println!("{:?}", &bot.config);
    bot.init().unwrap();
    bot.run().unwrap();
}
