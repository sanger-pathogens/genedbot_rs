#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate serde_json;
extern crate chrono;
extern crate config;
extern crate libflate;
extern crate mediawiki;
extern crate regex;
extern crate reqwest;

use bio::io::{gaf, gff};
use chrono::Local;
use config::{Config, File};
use libflate::gzip::Decoder;
use percent_encoding::percent_decode;
use regex::Regex;
use std::collections::{HashMap, HashSet};
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

//#[derive(Debug, Clone)]
//pub struct Literature {}
type Literature = String;

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
    gaf: HashMap<String, Vec<bio::io::gaf::Record>>,
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
    other_types: HashMap<String, HashMap<String, Option<bio::io::gff::Record>>>,
    orth_genedb2q: HashMap<String, String>,
    orth_genedb2q_taxon: HashMap<String, String>,
    parent2child: HashMap<String, Vec<(String, String)>>,
    paper2q: HashMap<(String, String), String>,
    aspects: HashMap<String, String>,
    go_term2q: HashMap<String, String>,
    pub specific_genes_only: Option<Vec<String>>,
}

impl GeneDBot {
    pub fn new() -> Self {
        Self {
            gff: HashMap::new(),
            gaf: HashMap::new(),
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
            other_types: HashMap::new(),
            orth_genedb2q: HashMap::new(),
            orth_genedb2q_taxon: HashMap::new(),
            parent2child: HashMap::new(),
            paper2q: HashMap::new(),
            go_term2q: HashMap::new(),
            specific_genes_only: None,
            aspects: vec![("P", "P682"), ("F", "P680"), ("C", "P681")]
                .iter()
                .map(|x| (x.0.to_string(), x.1.to_string()))
                .collect(),
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
            "ftp://ftp.sanger.ac.uk/pub/genedb/releases/latest/{}/{}.gaf.gz",
            self.species_key, self.species_key
        )
    }

    pub fn load_gff_file_from_url(&mut self, url: &str) -> Result<(), Box<Error>> {
        let mut orth_ids: HashSet<String> = HashSet::new();
        let mut res = reqwest::get(url)?;
        let decoder = Decoder::new(&mut res)?;
        let mut reader = gff::Reader::new(decoder, gff::GffType::GFF3);
        for element in reader.records() {
            match element {
                Ok(e) => {
                    self.process_gff_element(&e, &mut orth_ids);
                }
                _ => continue,
            }
        }

        if self.gff.is_empty() {
            panic!("Can't get GFF data from {}", url);
        }

        self.load_orthologs(orth_ids).unwrap(); //TODO
        Ok(())
    }

    fn load_orthologs(
        &mut self,
        mut orth_ids: HashSet<String>,
    ) -> Result<(), Box<::std::error::Error>> {
        if orth_ids.is_empty() {
            return Ok(());
        }

        if orth_ids.len() < 1000 {
            // Usually for testing
            let orth_ids: Vec<String> = orth_ids.drain().collect();
            for chunk in orth_ids.chunks(100) {
                let sparql = format!("SELECT ?q ?genedb ?taxon {{ VALUES ?genedb {{'{}'}} . ?q wdt:P3382 ?genedb ; wdt:P703 ?taxon }}",chunk.join("' '"));
                let sparql_result = self.api.sparql_query(&sparql)?;
                for b in sparql_result["results"]["bindings"].as_array().unwrap() {
                    let q = match b["q"]["value"].as_str() {
                        Some(s) => self.api.extract_entity_from_uri(s).unwrap(),
                        None => continue,
                    };
                    let taxon_q = match b["taxon"]["value"].as_str() {
                        Some(s) => self.api.extract_entity_from_uri(s).unwrap(),
                        None => continue,
                    };
                    let genedb = match b["genedb"]["value"].as_str() {
                        Some(s) => s.to_string(),
                        None => continue,
                    };
                    self.orth_genedb2q.insert(genedb.to_string(), q);
                    self.orth_genedb2q_taxon.insert(genedb.to_string(), taxon_q);
                }
            }
        } else {
            // Retrieven 'em all and let HashSet sort 'em out...
            let sparql = "SELECT ?q ?genedb ?taxon { ?q wdt:P3382 ?genedb ; wdt:P703 ?taxon }";
            let sparql_result = self.api.sparql_query(&sparql)?;
            for b in sparql_result["results"]["bindings"].as_array().unwrap() {
                let genedb = match b["genedb"]["value"].as_str() {
                    Some(s) => s.to_string(),
                    None => continue,
                };
                if !orth_ids.contains(&genedb) {
                    continue;
                }
                let q = match b["q"]["value"].as_str() {
                    Some(s) => self.api.extract_entity_from_uri(s).unwrap(),
                    None => continue,
                };
                let taxon_q = match b["taxon"]["value"].as_str() {
                    Some(s) => self.api.extract_entity_from_uri(s).unwrap(),
                    None => continue,
                };
                self.orth_genedb2q.insert(genedb.to_string(), q);
                self.orth_genedb2q_taxon.insert(genedb.to_string(), taxon_q);
            }
        }
        Ok(())
    }

    fn fix_id(&self, id: &str) -> String {
        lazy_static! {
            static ref RE1: Regex = Regex::new(r":.*$").unwrap();
            static ref RE2: Regex = Regex::new(r"\.\d$").unwrap();
        }
        let id2 = RE1.replace_all(id, "").into_owned();
        RE2.replace_all(&id2, "").into_owned()
    }

    fn set_other_types(&mut self, element: &bio::io::gff::Record, id: &str) {
        let id = self.fix_id(id);
        let o = match self.alternate_gene_subclasses.get(element.feature_type()) {
            Some(_) => Some(element.clone()),
            None => None,
        };
        self.other_types
            .entry(element.feature_type().to_string())
            .or_insert(HashMap::new())
            .insert(id, o);
    }

    fn process_gff_element(
        &mut self,
        element: &bio::io::gff::Record,
        orth_ids: &mut HashSet<String>,
    ) {
        lazy_static! {
            static ref VALID_GENE_TYPES: Vec<&'static str> = vec!["gene", "mRNA", "pseudogene"];
        }

        if element.attributes().contains_key("ID") {
            let id = element.attributes()["ID"].clone();
            self.gff.insert(id.clone(), element.clone());
            match element.attributes().get("Parent") {
                Some(parent_id) => self
                    .parent2child
                    .entry(parent_id.to_string())
                    .or_insert(vec![])
                    .push((id.clone(), element.feature_type().to_string())),
                None => {}
            }
        }

        if !VALID_GENE_TYPES.contains(&element.feature_type()) {
            match element.attributes().get("ID") {
                Some(id) => {
                    self.set_other_types(&element, id);
                    match element.attributes().get("Parent") {
                        Some(parent_id) => self.set_other_types(&element, parent_id),
                        None => {}
                    }
                }
                None => {}
            }
        }

        // Orthologs
        match &self.specific_genes_only {
            Some(genes) => match element.attributes().get("ID") {
                Some(id) => {
                    // TODO HACKISH
                    let mut id = id.clone();
                    id.pop();
                    id.pop();
                    if !genes.contains(&id) {
                        return;
                    }
                }
                None => return,
            },
            None => {}
        }
        self.get_orthologs_from_gff_element(&element)
            .iter()
            .for_each(|x| {
                orth_ids.insert(x.1.to_owned());
            });
    }

    // Returns (species,genedb_id)
    fn get_orthologs_from_gff_element(&self, gff: &bio::io::gff::Record) -> Vec<(String, String)> {
        lazy_static! {
            static ref RE_ORTH: Regex = Regex::new(r"^\s*(\S*):(\S+)").unwrap();
        }
        match gff.attributes().get("orthologous_to") {
            Some(orth) => {
                let orth = self.fix_attribute_value(orth);
                let v: Vec<&str> = orth.split(',').collect();
                v.iter()
                    .filter(|o| RE_ORTH.is_match(o))
                    .map(|o| {
                        RE_ORTH
                            .captures_iter(o)
                            .map(|m| return (m[1].to_string(), m[2].to_string()))
                            .next()
                    })
                    .filter(|o| o.is_some())
                    .map(|o| o.unwrap())
                    .collect()
            }
            None => vec![],
        }
    }

    fn fix_attribute_value(&self, s: &str) -> String {
        percent_decode(s.as_bytes())
            .decode_utf8()
            .unwrap()
            .to_string()
    }

    pub fn load_gaf_file_from_url(&mut self, url: &str) -> Result<(), Box<Error>> {
        println!("{}", url);
        let mut res = reqwest::get(url)?;
        let decoder = Decoder::new(&mut res)?;
        let mut reader = gaf::Reader::new(decoder, gaf::GafType::GAF2);
        for element in reader.records() {
            match element {
                Ok(e) => {
                    let id = e.db_object_id().to_string();
                    if !self.gaf.contains_key(&id) {
                        self.gaf.insert(id.clone(), vec![]);
                    }
                    self.gaf.get_mut(&id).unwrap().push(e);
                }
                _ => continue,
            }
        }
        if self.gaf.is_empty() {
            panic!("Can't get GAF data from {}", url);
        }
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

        let qualifiers = vec![];
        let mut new_item = Entity::new_empty_item();
        new_item.set_label(LocaleString::new("en", &(taxon_name + " reference genome")));
        new_item.add_claim(Statement::new_normal(
            Snak::new_item("P279", "Q7307127"),
            qualifiers.clone(),
            self.references(),
        ));
        new_item.add_claim(Statement::new_normal(
            Snak::new_item("P703", species_q.as_str()),
            qualifiers.clone(),
            self.references(),
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

    fn new_time_today(&self) -> Snak {
        let today = Local::now();
        Snak::new_time(
            "P813",
            &format!("{}", today.format("+%Y-%m-%dT00:00:00Z")),
            11,
        )
    }

    fn process_gene(&mut self, genedb_id: String) {
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

        let mut item = Entity::new_empty_item();
        let item_to_diff = match self.get_entity_for_genedb_id(&genedb_id) {
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
        let reference = Reference::new(vec![
            Snak::new_item("P248", "Q5531047"),
            self.new_time_today(),
        ]);
        let ga_quals = vec![
            Snak::new_item("P659", &self.genomic_assembly_q),
            Snak::new_item("P1057", &chr_q),
        ];

        let mut statements_to_create = vec![
            Snak::new_item("P31", gene_type.1),        // Instance of
            Snak::new_item("P703", &self.species_q()), // Found in:Species
            Snak::new_item("P1057", &chr_q),           // Chromosome
            Snak::new_string("P3382", &genedb_id),
        ];

        match gff.strand() {
            Some(strand) => match strand.strand_symbol() {
                "+" => statements_to_create.push(Snak::new_item("P2548", "Q22809680")),
                "-" => statements_to_create.push(Snak::new_item("P2548", "Q22809711")),
                _ => {}
            },
            _ => {}
        }

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

        let protein_entity_ids = self.process_proteins(&genedb_id);
        if protein_entity_ids.len() > 0 {
            statements_to_create.push(Snak::new_item("P279", "Q20747295")); // Subclass of:protein-coding gene

            // Encodes: protein
            for protein_q in &protein_entity_ids {
                statements_to_create.push(Snak::new_item("P688", &protein_q));
            }
        } else if gene_type.0 != "gene" {
            statements_to_create.push(Snak::new_item("P688", gene_type.1));
        } else {
            let mut subclass_found: bool = false;
            for subclass in &self.alternate_gene_subclasses {
                let class_name = subclass.0;
                let _class_q = subclass.1;
                let ct = match self.other_types.get(&class_name.to_owned()) {
                    Some(ct) => ct,
                    None => continue,
                };
                match ct.get(&genedb_id) {
                    Some(_feature_type) => {
                        subclass_found = true;
                        let mut fake_literature: Vec<Literature> = vec![];
                        self.process_product(
                            &gff,
                            &mut item,
                            &mut fake_literature,
                            &genedb_id,
                            &reference,
                        );
                    }
                    None => {}
                };
            }
            if !subclass_found {
                self.log(
                    &genedb_id,
                    &format!("No subclass found for {:?}", &gene_type),
                );
            }
        }

        // Orthologs
        if protein_entity_ids.len() > 0 {
            self.parent2child
                .get(&genedb_id)
                .unwrap_or(&vec![])
                .iter()
                .filter(|o| o.1 == "mRNA")
                .for_each(|o| match self.gff.get(&o.0) {
                    Some(protein) => {
                        self.get_orthologs_from_gff_element(&protein)
                            .iter()
                            .for_each(|(_species, protein_genedb_id)| {
                                match self.orth_genedb2q.get(protein_genedb_id) {
                                    Some(orth_q) => {
                                        match self.orth_genedb2q_taxon.get(protein_genedb_id) {
                                            Some(orth_q_taxon) => {
                                                item.add_claim(Statement::new_normal(
                                                    Snak::new_item("P684", orth_q),
                                                    vec![Snak::new_item("P703", orth_q_taxon)],
                                                    vec![reference.clone()],
                                                ));
                                            }
                                            None => {}
                                        }
                                    }
                                    None => {}
                                }
                            });
                    }
                    None => {}
                });
        }

        // Create simple (single-reference-only) statements
        statements_to_create.iter().for_each(|s| {
            item.add_claim(Statement::new_normal(
                s.to_owned(),
                vec![],
                vec![reference.clone()],
            ))
        });

        // Apply diff
        let my_props = vec![
            "P279",  // Subclass of (mostly to remove protein-coding gene)
            "P703",  // Found in taxon
            "P680",  // Molecular function
            "P681",  // Cell component
            "P682",  // Biological process
            "P684",  // Ortholog
            "P1057", // Chromosome
            "P2548", // Strand orientation
            "P644",  // Genomic start
            "P645",  // Genomic end
        ];

        let mut params = EntityDiffParams::none();
        params.labels = EntityDiffParam::some(&vec!["en"]);
        params.descriptions = EntityDiffParam::some(&vec!["en"]);
        params.aliases = EntityDiffParam::some(&vec!["en"]);
        params.claims.add = EntityDiffParamState::All;
        params.claims.alter = EntityDiffParamState::All;
        params.claims.remove = EntityDiffParamState::some(&my_props);
        params.references.list.push((
            EntityDiffParamState::some(&my_props),
            EntityDiffParamState::except(&vec!["P813"]),
        ));

        let diff = EntityDiff::new(&item_to_diff, &item, &params);
        println!("{}", diff.actions());

        // TODO apply diff/create new (plus, update internal matches and add protein=>gene)
        /*
        let _gene_q = match self.get_entity_id_for_genedb_id(&genedb_id) {
            Some(q) => q,
            None => return,
        };
        */

        //println!("{}", serde_json::to_string_pretty(&item).unwrap());
    }

    fn log(&self, genedb_id: &String, message: &str) {
        println!("!! {}: {}", genedb_id, message);
    }

    fn process_product(
        &self,
        _data: &bio::io::gff::Record,
        _item: &mut Entity,
        _literature: &mut Vec<Literature>,
        _genedb_id: &String,
        _references: &Reference,
    ) {
        // TODO
    }

    fn process_proteins(&mut self, genedb_id: &String) -> Vec<String> {
        let protein_genedb_ids: Vec<String> = self
            .parent2child
            .get(&genedb_id.clone())
            .unwrap_or(&vec![])
            .iter()
            .filter(|child| child.1 == "mRNA")
            .map(|child| child.0.to_owned())
            .collect();
        protein_genedb_ids
            .iter()
            .map(|protein_genedb_id| self.process_protein(&genedb_id, &protein_genedb_id))
            .filter(|entry| entry.is_some())
            .map(|entry| entry.unwrap())
            .collect()
    }

    fn process_protein(
        &mut self,
        _gene_genedb_id: &String,
        protein_genedb_id: &String,
    ) -> Option<String> {
        let gff = match self.gff.get(protein_genedb_id) {
            Some(gff) => gff.clone(),
            None => return None,
        };
        //let mut label = protein_genedb_id.clone();
        let mut _desc = String::from("");

        let mut item = Entity::new_empty_item();
        let _item_to_diff = match self.get_entity_for_genedb_id(&protein_genedb_id) {
            Some(i) => i.clone(),
            None => Entity::new_empty_item(),
        };

        let _reference = Reference::new(vec![
            Snak::new_item("P248", "Q5531047"),
            self.new_time_today(),
        ]);

        let mut literature: Vec<Literature> = vec![];

        match gff.attributes().get("literature") {
            Some(lit) => lit.split(',').for_each(|l| literature.push(l.to_string())),
            None => {}
        }

        item.set_label(LocaleString::new("en", &protein_genedb_id.clone()));

        self.add_go_annotation(&mut item, &gff, &mut literature);

        None
    }

    fn add_go_annotation(
        &mut self,
        item: &mut Entity,
        gff: &bio::io::gff::Record,
        literature: &mut Vec<Literature>,
    ) {
        let protein_genedb_id = gff.attributes()["ID"].clone();
        let gaf = match self.gaf.get(&protein_genedb_id) {
            Some(gaf) => gaf.clone(),
            None => return,
        };
        let mut new_go_claims: HashMap<String, (Snak, Vec<Reference>, Vec<Snak>)> = HashMap::new();
        for ga in gaf {
            let go_term = ga.go_id().to_string();
            let go_q = match self.get_item_for_go_term(&go_term) {
                Some(q) => q,
                None => {
                    println!("No Wikidata item for GO term: '{}'", &go_term);
                    continue;
                }
            };

            let aspect = ga.aspect().to_string();
            let aspect_p = match self.aspects.get(&aspect) {
                Some(p) => p.clone(),
                None => {
                    println!("Unknown aspect: '{}'", &aspect);
                    continue;
                }
            };

            let evidence_code = ga.evidence_code().to_string();
            let evidence_code_q = match self.evidence_codes.get(&evidence_code) {
                Some(q) => q.clone(),
                None => {
                    println!("Unknown evidence code: '{}'", &evidence_code);
                    continue;
                }
            };

            let mut literature_sources: Vec<Snak> = vec![];
            for (k, v) in ga.db_ref().iter() {
                // Literature
                match k.as_str() {
                    "GO_REF" => literature_sources.push(Snak::new_string("P854", &format!("https://github.com/geneontology/go-site/blob/master/metadata/gorefs/goref-{}.md",&v))),
                    "PMID" => match self.get_or_create_paper_item(k, v) {
                        Some(paper_q) => literature_sources.push(Snak::new_item("P248", &paper_q)),
                        None => println!("Can't find/create paper for '{}/{}'", &k, &v),
                    },
                    other => {
                        println!("> {}",other);
                        for w_f in ga.with_from() {
                            let parts : Vec<&str> = w_f.split(':').collect();
                            if parts.len() == 2 && parts[0] == "InterPro" {
                                literature_sources.push(Snak::new_string("P2926", parts[1]));
                            }
                        }
                    }
                }

                // Qualifiers
                let mut qualifiers = vec![Snak::new_item("P459", &evidence_code_q)];
                for qual in ga.qualifier() {
                    if qual == "NOT" {
                        qualifiers.push(Snak::new_item("P6477", "Q186290"));
                    }
                }

                // Qualifiers from with_from
                for w_f in ga.with_from() {
                    let parts: Vec<&str> = w_f.split(':').collect();
                    match self.get_with_from_qualifier(&parts) {
                        Some(snak) => qualifiers.push(snak),
                        None => {}
                    }
                }

                // References for GO terms
                let references: Vec<Reference> = literature_sources
                    .iter()
                    .map(|ls| {
                        Reference::new(vec![
                            ls.clone(),
                            Snak::new_item("P1640", "Q5531047"),
                            self.new_time_today(),
                        ])
                    })
                    .collect();

                let new_claim_key = json!([aspect_p.clone(), go_q.clone(), qualifiers.clone()]);
                let new_claim_key = new_claim_key.as_str().unwrap().to_string();
                if new_go_claims.contains_key(&new_claim_key) {
                    let ngc = new_go_claims.get_mut(&new_claim_key).unwrap();
                    for r in references {
                        ngc.1.push(r);
                    }
                } else {
                    new_go_claims.insert(
                        new_claim_key,
                        (
                            Snak::new_item(aspect_p.clone(), go_q.clone()),
                            references,
                            qualifiers,
                        ),
                    );
                }

                for (k, v) in ga.db_ref().iter() {
                    literature.push(format!("{}|{}", k, v).to_string());
                }

                // Label/aliases
                if !ga.db_object_name().is_empty() {
                    match item.label_in_locale("en") {
                        Some(l) => {
                            if l == protein_genedb_id {
                                item.set_label(LocaleString::new(
                                    "en",
                                    &ga.db_object_name().to_string(),
                                ));
                                item.add_alias(LocaleString::new("en", &protein_genedb_id.clone()));
                            }
                        }
                        None => {}
                    }
                }

                // Aliases
                ga.db_object_synonym().iter().for_each(|syn| {
                    syn.split(',').for_each(|s| {
                        item.add_alias(LocaleString::new("en", &self.fix_alias_name(s)));
                    });
                });
            }
        }

        new_go_claims.iter().for_each(|(_key, claim)| {
            item.add_claim(Statement::new_normal(
                claim.0.to_owned(),
                claim.2.to_owned(),
                claim.1.to_owned(),
            ));
        });
    }

    fn get_with_from_qualifier(&self, parts: &Vec<&str>) -> Option<Snak> {
        if parts.len() != 1 {
            return None;
        }
        match parts[0] {
            "Pfam" => Some(Snak::new_string("P3519", parts[1])),
            "Rfam" => Some(Snak::new_string("P3523", parts[1])),
            "GeneDB" => Some(Snak::new_string("P3382", parts[1])),
            "UniProt" => Some(Snak::new_string("P352", parts[1])),
            "InterPro" => Some(Snak::new_string("P2926", parts[1])),
            "PANTHER" => Some(Snak::new_string(
                "P973",
                &format!(
                    "http://www.pantherdb.org/panther/family.do?clsAccession={}",
                    parts[1]
                ),
            )),
            "CBS" => match parts[1] {
                "TMHMM" => Some(Snak::new_item("P2283", &TMHMM_Q)),
                _ => None,
            },
            _ => None,
        }
    }

    fn get_or_create_paper_item(&mut self, k: &String, v: &String) -> Option<String> {
        if self.paper2q.contains_key(&(k.to_string(), v.to_string())) {
            return Some(
                self.paper2q
                    .get(&(k.to_string(), v.to_string()))
                    .unwrap()
                    .to_string(),
            );
        }
        match k.as_str() {
            "PMID" => {
                let sparql = format!("SELECT ?q {{ ?q wdt:P698 '{}' }}", &v);
                let sparql_result = self.api.sparql_query(&sparql).ok()?;
                let items = self.api.entities_from_sparql_result(&sparql_result, "q");
                match items.len() {
                    1 => {
                        self.paper2q
                            .insert((k.to_string(), v.to_string()), items[0].clone());
                        Some(items[0].clone())
                    }
                    _ => None,
                }
            }
            other => {
                println!("Unknown paper source: '{}'", &other);
                None
            }
        }
    }

    fn get_item_for_go_term(&mut self, go_term: &String) -> Option<String> {
        if self.go_term2q.contains_key(&go_term.clone()) {
            return Some(self.go_term2q.get(&go_term.clone()).unwrap().to_string());
        }
        let sparql = format!("SELECT ?q {{ ?q wdt:P686 '{}' }}", &go_term);
        let sparql_result = self.api.sparql_query(&sparql).ok()?;
        for b in sparql_result["results"]["bindings"].as_array().unwrap() {
            let q = match b["q"]["value"].as_str() {
                Some(s) => self.api.extract_entity_from_uri(s).unwrap(),
                None => continue,
            };
            self.go_term2q.insert(go_term.clone(), q.clone());
            return Some(q);
        }
        None
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
        //let species_q = self.species_q();
        //let _species_i = self.ec.load_entity(&self.api, species_q)?;
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
