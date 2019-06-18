#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate serde_json;
extern crate chrono;
extern crate config;
extern crate libflate;
extern crate mediawiki;
extern crate papers;
extern crate regex;
extern crate reqwest;

//use multimap::MultiMap;
//use std::{thread, time};
use chrono::Local;
use config::{Config, File};
use percent_encoding::percent_decode;
use regex::Regex;
use std::collections::{HashMap, HashSet};
use std::env;
use std::{error::Error, fmt};
use wikibase::entity_container::*;
use wikibase::entity_diff::*;
use wikibase::*;

pub mod gene;
pub mod literature;
pub mod loader;
pub mod protein;

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
pub type Literature = String;

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
    simulate: bool,
    verbose: bool,
    product_term_becomes_label: bool,
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
    xref2prop: HashMap<String, String>,
    aspects: HashMap<String, String>,
    go_term2q: HashMap<String, String>,
    genes2load: Vec<String>,
    pub specific_genes_only: Option<Vec<String>>,
    papers: literature::Papers,
}

impl GeneDBot {
    pub fn new() -> Self {
        let mut builder = reqwest::ClientBuilder::new();
        match std::env::var("http_proxy") {
            Ok(proxy_url) => {
                if !proxy_url.is_empty() {
                    builder = builder.proxy(reqwest::Proxy::http(proxy_url.as_str()).unwrap());
                }
            }
            _ => {}
        }
        match std::env::var("https_proxy") {
            Ok(proxy_url) => {
                if !proxy_url.is_empty() {
                    builder = builder.proxy(reqwest::Proxy::https(proxy_url.as_str()).unwrap());
                }
            }
            _ => {}
        }
        let api =
            mediawiki::api::Api::new_from_builder("https://www.wikidata.org/w/api.php", builder)
                .unwrap();
        Self {
            simulate: false,
            verbose: false,
            product_term_becomes_label: true,
            gff: HashMap::new(),
            gaf: HashMap::new(),
            config: GeneDBotConfig::new_from_json(&json!({})),
            species_key: "".to_string(),
            genomic_assembly_q: "".to_string(),
            ec: EntityContainer::new(),
            api: api.clone(),
            chr2q: HashMap::new(),
            genedb2q: HashMap::new(),
            protein_genedb2q: HashMap::new(),
            evidence_codes_labels: HashMap::new(),
            evidence_codes: HashMap::new(),
            other_types: HashMap::new(),
            orth_genedb2q: HashMap::new(),
            orth_genedb2q_taxon: HashMap::new(),
            parent2child: HashMap::new(),
            go_term2q: HashMap::new(),
            specific_genes_only: None,
            papers: literature::Papers::new(&api),
            genes2load: vec![],
            aspects: vec![("P", "P682"), ("F", "P680"), ("C", "P681")]
                .iter()
                .map(|x| (x.0.to_string(), x.1.to_string()))
                .collect(),
            xref2prop: vec![("UniProtKB", "P352")]
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

    fn log(&self, genedb_id: &String, message: &str) {
        println!("{}: {}", genedb_id, message);
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

    fn references(&self) -> Vec<Reference> {
        vec![Reference::new(vec![
            Snak::new_item("P248", "Q5531047"),
            self.new_time_today(),
        ])]
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

    fn normalize_evidence_label(&self, s: &String) -> String {
        s.trim().to_lowercase()
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

    fn get_or_create_chromosome_entity(&mut self, id: &str) -> Option<String> {
        match self.chr2q.get(id) {
            Some(q) => return Some(q.clone()),
            None => {}
        }

        let mut new_item = Entity::new_empty_item();
        new_item.set_label(LocaleString::new("en", id));
        new_item.add_claim(Statement::new_normal(
            Snak::new_item("P31", "Q37748"),
            vec![],
            self.references(),
        ));
        new_item.add_claim(Statement::new_normal(
            Snak::new_item("P703", self.species_q().as_str()),
            vec![],
            self.references(),
        ));

        let params = EntityDiffParams::all();
        let mut diff = EntityDiff::new(&Entity::new_empty_item(), &new_item, &params);
        diff.set_edit_summary(self.get_edit_summary());
        match self.ec.apply_diff(&mut self.api, &diff) {
            Some(q) => {
                self.chr2q.insert(id.to_string(), q.clone());
                Some(q)
            }
            None => panic!("Could not create chromosome item for '{}'", &id),
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

    fn get_edit_summary(&self) -> Option<String> {
        Some("Syncing to GeneDB (V3)".to_string())
    }

    // Returns (species,genedb_id)
    fn get_orthologs_from_gff_element(&self, gff: &bio::io::gff::Record) -> Vec<(String, String)> {
        lazy_static! {
            static ref RE_ORTH: Regex = Regex::new(r"^\s*(\S*):(\S+)").unwrap();
        }

        match gff.attributes().get_vec("orthologous_to") {
            Some(orths) => {
                let mut ret: Vec<(String, String)> = vec![];
                for orth in orths {
                    let orth = self.fix_attribute_value(orth);
                    let v: Vec<&str> = orth.split("type=orthologous_to,").collect();
                    let mut result: Vec<(String, String)> = v
                        .iter()
                        .filter(|o| RE_ORTH.is_match(o))
                        .map(|o| {
                            RE_ORTH
                                .captures_iter(o)
                                .map(|m| return (m[1].to_string(), m[2].to_string()))
                                .next()
                        })
                        .filter(|o| o.is_some())
                        .map(|o| o.unwrap())
                        .collect();
                    ret.append(&mut result);
                }
                ret
            }
            None => vec![],
        }
    }

    fn fix_attribute_value(&self, s: &str) -> String {
        percent_decode(s.as_bytes())
            .decode_utf8()
            .unwrap()
            .trim_end_matches(';')
            .to_string()
    }

    fn is_product_type(&self, s: &str) -> bool {
        s == "mRNA" || s == "pseudogenic_transcript"
    }

    fn process_orthologs(&mut self, item: &mut Entity, genedb_id: &String, reference: &Reference) {
        let mut had_that: HashSet<String> = HashSet::new();
        self.parent2child
            .get(genedb_id)
            .unwrap_or(&vec![])
            .iter()
            .filter(|o| self.is_product_type(&o.1))
            .for_each(|o| match self.gff.get(&o.0) {
                Some(protein) => {
                    self.get_orthologs_from_gff_element(&protein)
                        .iter()
                        .for_each(|(_species, protein_genedb_id)| {
                            match self.orth_genedb2q.get(protein_genedb_id) {
                                Some(orth_q) => {
                                    if had_that.contains(orth_q) {
                                        return;
                                    }
                                    had_that.insert(orth_q.to_string());
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

    fn is_item(&self, q: &String) -> bool {
        lazy_static! {
            static ref RE1: Regex = Regex::new(r"^Q\d+$").unwrap();
        }
        RE1.is_match(q)
    }

    fn process_product_controlled_curation(
        &mut self,
        gff: &bio::io::gff::Record,
        item: &mut Entity,
        literature: &mut HashSet<Literature>,
        genedb_id: &String,
        _reference: &Reference,
    ) {
        lazy_static! {
            static ref RE1: Regex = Regex::new(r"^(.+?)=(.+)$").unwrap();
        }

        match gff.attributes().get_vec("controlled_curation") {
            Some(parts) => {
                for part in parts {
                    let mut kv: HashMap<String, String> = HashMap::new();
                    let part = self.fix_attribute_value(part);
                    let subparts: Vec<&str> = part.split(';').collect();
                    for subpart in subparts {
                        RE1.captures_iter(&subpart).for_each(|m| {
                            kv.insert(m[1].to_string(), m[2].to_string());
                        });
                    }

                    let mut statement;

                    // Is there a recognizable term?
                    match kv.get("term") {
                        Some(s) => match s.as_str() {
                            "gene deletion phenotype: essential" => {
                                statement = Statement::new_normal(
                                    Snak::new_item("P279", "Q17119234"), // Subclass of: Essential gene
                                    vec![],
                                    vec![],
                                );
                            }
                            "dispensable" => {
                                statement = Statement::new_normal(
                                    Snak::new_item("P279", "Q63092631"), // Subclass of: Dispensable gene
                                    vec![],
                                    vec![],
                                );
                            }
                            other => {
                                self.log(
                                    &genedb_id,
                                    &format!("unknown controlled curation term '{}'", &other),
                                );
                                continue;
                            }
                        },
                        _ => {
                            // No term
                            continue;
                        }
                    }

                    // Valid term confirmed, new statement primed
                    let mut reference = Reference::new(vec![
                        Snak::new_item("P1640", "Q5531047"),
                        self.new_time_today(),
                    ]);

                    // Add evidence code
                    match kv.get("evidence") {
                        Some(evidence_text) => {
                            match self
                                .evidence_codes_labels
                                .get(&self.normalize_evidence_label(evidence_text))
                            {
                                Some(ecq) => {
                                    statement.add_qualifier_snak(Snak::new_item("P459", ecq));
                                }
                                None => {
                                    self.log(
                                        &genedb_id,
                                        &format!("Unrecognized evidence: '{}'", &evidence_text),
                                    );
                                }
                            }
                        }
                        None => {}
                    }

                    // Literature
                    match kv.get("db_xref") {
                        Some(xref) => {
                            let parts: Vec<&str> = xref.split(':').collect();
                            if parts.len() == 2 {
                                let paper_item = self.papers.get_or_create_item(
                                    &parts[0].to_string(),
                                    &parts[1].to_string(),
                                );
                                match paper_item {
                                    Some(q) => {
                                        let mut snaks = reference.snaks().clone();
                                        snaks.push(Snak::new_item("P248", &q));
                                        reference.set_snaks(snaks);
                                        literature.insert(
                                            format!("{}:{}", &parts[0], &parts[1]).to_string(),
                                        );
                                    }
                                    None => {
                                        self.log(
                                            &genedb_id,
                                            &format!(
                                                "Unrecognized literature: '{}:{}'",
                                                &parts[0], &parts[1]
                                            ),
                                        );
                                    }
                                }
                            }
                        }
                        None => {}
                    }

                    // Add the reference to the statement
                    statement.set_references(vec![reference]);

                    // Add statement to item
                    item.add_claim(statement);
                }
            }
            None => {}
        }
    }

    fn process_product(
        &mut self,
        gff: &bio::io::gff::Record,
        item: &mut Entity,
        literature: &mut HashSet<Literature>,
        genedb_id: &String,
        reference: &Reference,
    ) {
        lazy_static! {
            static ref RE1: Regex = Regex::new(r"^(.+?)=(.+)$").unwrap();
            static ref RE2: Regex = Regex::new(r"^term=(.+)$").unwrap();
        }

        self.process_product_controlled_curation(gff, item, literature, genedb_id, reference);

        match gff.attributes().get_vec("Dbxref") {
            Some(xrefs) => {
                xrefs.iter().for_each(|xref| {
                    let parts: Vec<&str> = xref.split(':').collect();
                    if parts.len() == 2 {
                        match self.xref2prop.get(parts[0]) {
                            Some(prop) => item.add_claim(Statement::new_normal(
                                Snak::new_string(prop.clone(), parts[1].to_string()),
                                vec![],
                                vec![reference.clone()],
                            )),
                            None => {}
                        }
                    }
                });
            }
            None => {}
        }

        let mut apk: HashMap<String, String> = HashMap::new();
        let products = match gff.attributes().get_vec("product") {
            Some(products) => products,
            None => return,
        };

        products.iter().for_each(|product| {
            let parts: Vec<&str> = product.split(',').collect();
            for part in parts {
                let part = self.fix_attribute_value(part);
                RE1.captures_iter(&part).for_each(|m| {
                    apk.insert(m[1].to_string(), m[2].to_string());
                });
                RE2.captures_iter(&part).for_each(|m| {
                    let m_fixed = self.fix_alias_name(&m[1]);
                    match item.label_in_locale("en") {
                        Some(label) => {
                            if self.product_term_becomes_label {
                                if label == genedb_id {
                                    item.set_label(LocaleString::new("en", &m_fixed));
                                } else {
                                    item.add_alias(LocaleString::new("en", &m_fixed));
                                }
                            } else {
                                item.add_alias(LocaleString::new("en", &m_fixed));
                            }
                        }
                        None => {
                            item.set_label(LocaleString::new("en", &m_fixed));
                        }
                    }
                });
            }
        });

        self.set_evidence(&apk, item, literature);
    }

    fn set_evidence(
        &mut self,
        apk: &HashMap<String, String>,
        item: &mut Entity,
        literature: &mut HashSet<Literature>,
    ) {
        if !apk.contains_key("evidence") {
            return;
        } // Why is this?

        let reference = Reference::new(vec![
            Snak::new_item("P1640", "Q5531047"),
            self.new_time_today(),
        ]);
        let mut qualifiers = vec![];

        match apk.get("evidence") {
            Some(evidence) => {
                match self
                    .evidence_codes_labels
                    .get(&self.normalize_evidence_label(evidence))
                {
                    Some(ecq) => qualifiers.push(Snak::new_item("P459", ecq)),
                    None => {}
                }
            }
            None => {}
        }
        match apk.get("term") {
            Some(term) => qualifiers.push(Snak::new_string("P1810", term)),
            None => {}
        }
        match apk.get("with") {
            Some(with) => {
                let parts: Vec<&str> = with.split('|').collect();
                match self.get_with_from_qualifier(&parts) {
                    Some(qual) => qualifiers.push(qual),
                    None => {}
                }
            }
            None => {}
        }

        lazy_static! {
            static ref RE1: Regex = Regex::new(r"^(PMID):(.+)$").unwrap();
        }

        let mut lit_q: Option<String> = None;
        match apk.get("db_xref") {
            Some(xref) => {
                RE1.captures_iter(xref).for_each(|m| {
                    let k = m[1].to_string();
                    let v = m[2].to_string();
                    literature.insert(format!("{}:{}", k, v).to_string());
                    lit_q = self.papers.get_or_create_item(&k, &v);
                });
            }
            None => {}
        }

        match lit_q {
            Some(q) => item.add_claim(Statement::new_normal(
                Snak::new_item("P1343", &q),
                qualifiers,
                vec![reference],
            )),
            None => item.add_claim(Statement::new_normal(
                Snak::new_unknown_value("P1343", SnakDataType::SomeValue),
                qualifiers,
                vec![reference],
            )),
        }
    }

    fn process_proteins(&mut self, genedb_id: &String) -> Vec<String> {
        let protein_genedb_ids: Vec<String> = self
            .parent2child
            .get(&genedb_id.clone())
            .unwrap_or(&vec![])
            .iter()
            .filter(|child| self.is_product_type(&child.1))
            .map(|child| child.0.to_owned())
            .collect();
        protein_genedb_ids
            .iter()
            .map(|protein_genedb_id| {
                protein::process(self, &genedb_id, &protein_genedb_id)
                //self.process_protein(&genedb_id, &protein_genedb_id)
            })
            .filter(|entry| entry.is_some())
            .map(|entry| entry.unwrap())
            .collect()
    }

    fn get_with_from_qualifier(&self, parts: &Vec<&str>) -> Option<Snak> {
        if parts.len() != 2 {
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

    fn get_gene_ids_to_process(&self) -> Vec<String> {
        match &self.specific_genes_only {
            Some(genes) => genes.clone(),
            None => {
                let mut ret: Vec<String> =
                    self.genedb2q.iter().map(|(id, _)| id.to_owned()).collect();
                ret.extend(self.genes2load.iter().cloned());
                ret.sort();
                ret.dedup();
                ret
            }
        }
    }

    pub fn run(&mut self) -> Result<(), Box<Error>> {
        self.get_gene_ids_to_process()
            .iter()
            .for_each(|genedb_id| gene::process(self, genedb_id.to_string()));
        Ok(())
    }

    pub fn init(&mut self) -> Result<(), Box<Error>> {
        loader::init(self)
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
    settings.merge(File::with_name("bot.ini")).unwrap();
    let lgname = settings.get_str("user.user").unwrap();
    let lgpass = settings.get_str("user.pass").unwrap();

    let mut bot = GeneDBot::new();
    bot.api.set_user_agent("GeneDBot/3.0");
    bot.api.set_edit_delay(Some(500)); // Half a second between edits
    if args.len() == 3 {
        bot.specific_genes_only = Some(vec![args[2].to_string()]);
    }
    bot.api.login(lgname, lgpass).unwrap();
    bot.load_config_file(species_key)
        .expect("Can't load config file");
    bot.init().unwrap();
    bot.run().unwrap();
}
