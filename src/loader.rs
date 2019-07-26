use crate::genedbot::*;
//use reqwest::header::USER_AGENT;
//use crate::{GeneDBot, Toolbox};
use bio::io::{gaf, gff};
use libflate::gzip::Decoder;
use regex::Regex;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use wikibase::entity_diff::*;
use wikibase::*;

// These work, but are slow for some reason, so tests fail through timeout
pub static TEST_URL_JSON: &str =
    "http://raw.githubusercontent.com/sanger-pathogens/genedbot_rs/master/test_files/dummy.json";
pub static TEST_URL_GFF_GZ: &str =
    "http://raw.githubusercontent.com/sanger-pathogens/genedbot_rs/master/test_files/test.gff.gz";
pub static TEST_URL_GAF_GZ: &str =
    "http://raw.githubusercontent.com/sanger-pathogens/genedbot_rs/master/test_files/test.gaf.gz";
/*
// Alternative
pub static TEST_URL_JSON: &str = "http://magnusmanske.de/genedbot/dummy.json";
pub static TEST_URL_GFF_GZ: &str = "http://magnusmanske.de/genedbot/test.gff.gz";
pub static TEST_URL_GAF_GZ: &str = "http://magnusmanske.de/genedbot/test.gaf.gz";
*/

pub fn init(bot: &mut GeneDBot) -> Result<(), Box<Error>> {
    load_gff_file(bot)?; //.expect(&format!("Can't load GFF file '{}'", gff_url(bot)));
    load_gaf_file(bot)?; //.expect(&format!("Can't load GAF file '{}'", gaf_url(bot)));
    find_genomic_assembly(bot, true)?;
    load_basic_items(bot)?;
    Ok(())
}

pub fn get_text_from_url(url: &str) -> Result<String, reqwest::Error> {
    Ok(reqwest::get(url)?.text()?)
}

pub fn get_json_from_url(url: &str) -> Result<serde_json::Value, reqwest::Error> {
    Ok(reqwest::get(url)?.json()?)
}

pub fn gff_url(bot: &GeneDBot) -> String {
    format!(
        "ftp://ftp.sanger.ac.uk/pub/genedb/releases/latest/{}/{}.gff.gz",
        bot.species_key, bot.species_key
    )
}

pub fn gaf_url(bot: &GeneDBot) -> String {
    format!(
        "ftp://ftp.sanger.ac.uk/pub/genedb/releases/latest/{}/{}.gaf.gz",
        bot.species_key, bot.species_key
    )
}

pub fn load_gff_file_from_url(bot: &mut GeneDBot, url: &str) -> Result<(), Box<Error>> {
    let mut orth_ids: HashSet<String> = HashSet::new();
    /*
    let builder = GeneDBot::get_builder();
    let client = builder.gzip(false).build()?; // timeout(Duration::from_secs(10)).
    let mut res = client.get(url).header(USER_AGENT, "genedbot").send()?;
    let request_builder = client.get(url);
    */
    let mut res = reqwest::get(url)?;
    let decoder = Decoder::new(&mut res)?;
    let mut reader = gff::Reader::new(decoder, gff::GffType::GFF3);
    for element in reader.records() {
        match element {
            Ok(e) => {
                process_gff_element(bot, &e, &mut orth_ids);
            }
            _ => continue,
        }
    }

    if bot.gff.is_empty() {
        return Err(From::from(format!("Can't get GFF data from {}", url)));
    }
    bot.orthologs.load(&bot.api, orth_ids)
}

fn fix_id(id: &str) -> String {
    lazy_static! {
        static ref RE1: Regex = Regex::new(r":.*$").unwrap();
        static ref RE2: Regex = Regex::new(r"\.\d$").unwrap();
    }
    let id2 = RE1.replace_all(id, "").into_owned();
    RE2.replace_all(&id2, "").into_owned()
}

fn set_other_types(bot: &mut GeneDBot, element: &bio::io::gff::Record, id: &str) {
    let id = fix_id(id);
    let o = match bot.alternate_gene_subclasses.get(element.feature_type()) {
        Some(_) => Some(element.clone()),
        None => None,
    };
    bot.other_types
        .entry(element.feature_type().to_string())
        .or_insert(HashMap::new())
        .insert(id, o);
}

fn process_gff_element(
    bot: &mut GeneDBot,
    element: &bio::io::gff::Record,
    orth_ids: &mut HashSet<String>,
) {
    lazy_static! {
        static ref VALID_GENE_TYPES: Vec<&'static str> = vec!["gene", "mRNA", "pseudogene"];
    }

    if element.attributes().contains_key("ID") {
        let id = element.attributes()["ID"].clone();
        bot.gff.insert(id.clone(), element.clone());
        match element.attributes().get("Parent") {
            Some(parent_id) => bot
                .parent2child
                .entry(parent_id.to_string())
                .or_insert(vec![])
                .push((id.clone(), element.feature_type().to_string())),
            None => match element.feature_type().to_string().as_str() {
                "gene" => bot.genes2load.push(id.clone()),
                "pseudogene" => bot.genes2load.push(id.clone()),
                _ => {}
            },
        }
    }

    if !VALID_GENE_TYPES.contains(&element.feature_type()) {
        match element.attributes().get("ID") {
            Some(id) => {
                set_other_types(bot, &element, id);
                match element.attributes().get("Parent") {
                    Some(parent_id) => set_other_types(bot, &element, parent_id),
                    None => {}
                }
            }
            None => {}
        }
    }

    // Orthologs
    match &bot.specific_genes_only {
        Some(genes) => match element.attributes().get("ID") {
            Some(id) => {
                let mut ok = false;

                // TODO HACKISH but only relevant for testing, so...
                let mut id2 = id.clone();
                id2.pop();
                id2.pop();
                if genes.contains(&id2) {
                    ok = true;
                }

                if id.ends_with(":mRNA") {
                    let id2 = &id[..(id.len() - 5)];
                    if genes.contains(&id2.to_string()) {
                        ok = true;
                    }
                }

                if !ok {
                    return;
                }
            }
            None => return,
        },
        None => {}
    }

    bot.orthologs
        .get_from_gff_element(&element)
        .iter()
        .for_each(|x| {
            orth_ids.insert(x.1.to_owned());
        });
}

pub fn load_gaf_file_from_url(bot: &mut GeneDBot, url: &str) -> Result<(), Box<Error>> {
    let mut res = reqwest::get(url)?;
    let decoder = Decoder::new(&mut res)?;
    let mut reader = gaf::Reader::new(decoder, gaf::GafType::GAF2);
    for element in reader.records() {
        match element {
            Ok(e) => {
                let id = e.db_object_id().to_string();
                if !bot.gaf.contains_key(&id) {
                    bot.gaf.insert(id.clone(), vec![]);
                }
                bot.gaf.get_mut(&id).unwrap().push(e);
            }
            _ => continue,
        }
    }
    if bot.gaf.is_empty() {
        return Err(From::from(format!("Can't get GAF data from {}", url)));
    }
    Ok(())
}

pub fn load_gff_file(bot: &mut GeneDBot) -> Result<(), Box<Error>> {
    load_gff_file_from_url(bot, gff_url(bot).as_str())
}

pub fn load_gaf_file(bot: &mut GeneDBot) -> Result<(), Box<Error>> {
    load_gaf_file_from_url(bot, gaf_url(bot).as_str())
}

fn create_genomic_assembly_item(bot: &mut GeneDBot) -> Result<Entity, Box<Error>> {
    let species_q = bot.species_q();
    let species_i = bot.ec.load_entity(&bot.api, species_q.clone())?;
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
        bot.references(),
    ));
    new_item.add_claim(Statement::new_normal(
        Snak::new_item("P703", species_q.as_str()),
        qualifiers.clone(),
        bot.references(),
    ));
    Ok(new_item)
}

fn create_genomic_assembly(bot: &mut GeneDBot) -> Result<String, Box<Error>> {
    let new_item = create_genomic_assembly_item(bot)?;
    let params = EntityDiffParams::all();
    let mut diff = EntityDiff::new(&Entity::new_empty_item(), &new_item, &params);
    diff.set_edit_summary(bot.get_edit_summary());
    match bot.ec.apply_diff(&mut bot.api, &diff) {
        Some(q) => Ok(q),
        None => Err(From::from("Could not create genomic assembly item")),
    }
}

fn find_genomic_assembly(bot: &mut GeneDBot, do_create_if_missing: bool) -> Result<(), Box<Error>> {
    let sparql = format!(
        "SELECT ?q {{ ?q wdt:P279 wd:Q7307127 ; wdt:P703 wd:{} }}",
        bot.species_q()
    );
    let res = bot.api.sparql_query(&sparql)?;
    let candidates = bot.api.entities_from_sparql_result(&res, "q");
    bot.genomic_assembly_q = match candidates.len() {
        0 => {
            if do_create_if_missing {
                create_genomic_assembly(bot)?
            } else {
                return Err(From::from(format!(
                    "Can't find genomic assembly for {}, and not allowed to create one:\n{}",
                    &bot.species_q(),
                    &sparql
                )));
            }
        }
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

fn load_basic_items_chr(bot: &mut GeneDBot) -> Result<(), Box<Error>> {
    let sparql = format!(
        "SELECT ?q {{ ?q wdt:P31 wd:Q37748 ; wdt:P703 wd:{} }}",
        &bot.species_q()
    );
    let res = bot.api.sparql_query(&sparql)?;
    let mut items = bot.api.entities_from_sparql_result(&res, "q");
    items.push(bot.species_q());
    items.push(TMHMM_Q.to_string());
    bot.ec.load_entities(&bot.api, &items)?;

    items
        .iter()
        .for_each(|q| match bot.ec.get_entity(q.to_string()) {
            Some(i) => {
                if i.has_target_entity("P31", "Q37748") {
                    match i.label_in_locale("en") {
                        Some(label) => {
                            bot.chr2q.insert(label.to_string(), q.to_string());
                        }
                        None => {}
                    }
                }
            }
            None => {}
        });
    Ok(())
}

fn sparql_result_to_pairs(
    api: &mut mediawiki::api::Api,
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
            let v1 = api
                .extract_entity_from_uri(v1)
                .unwrap_or(v1.to_string())
                .to_string();
            let v2 = b[k2]["value"].as_str().unwrap();
            let v2 = api
                .extract_entity_from_uri(v2)
                .unwrap_or(v2.to_string())
                .to_string();
            (v1, v2)
        })
        .collect()
}

pub fn load_basic_items_genes(bot: &mut GeneDBot) -> Result<(), Box<Error>> {
    let mut species_list = vec![bot.species_q()];
    match bot.parent_taxon_q() {
        Some(parent_q) => species_list.push(parent_q),
        None => {}
    }
    let species_list = format!(" VALUES ?species {{ wd:{} }} ", species_list.join(" wd:"));

    let alternate_gene_subclasses_values: Vec<String> = bot
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
    let res = bot.api.sparql_query(&sparql)?;
    bot.genedb2q = sparql_result_to_pairs(&mut bot.api, &res["results"]["bindings"], "genedb", "q")
        .into_iter()
        .collect();

    // Proteins
    let sparql = format!("SELECT DISTINCT ?q ?genedb {{ {} . ?q wdt:P31 wd:Q8054 ; wdt:P703 ?species ; wdt:P3382 ?genedb }}",&species_list) ;
    let res = bot.api.sparql_query(&sparql)?;
    bot.protein_genedb2q =
        sparql_result_to_pairs(&mut bot.api, &res["results"]["bindings"], "genedb", "q")
            .into_iter()
            .collect();
    Ok(())
}

fn get_gene_entities_to_process(bot: &GeneDBot) -> Vec<String> {
    match &bot.specific_genes_only {
        Some(_genes) => vec![],
        None => bot
            .genedb2q
            .iter()
            .map(|(_, v)| v.to_owned())
            .chain(bot.protein_genedb2q.iter().map(|(_, v)| v.to_owned()))
            .collect(),
    }
}

pub fn load_basic_items_entities(bot: &mut GeneDBot) -> Result<(), Box<Error>> {
    let mut items_to_load: Vec<String> = get_gene_entities_to_process(bot).to_vec();
    bot.evidence
        .code2q
        .iter()
        .map(|(_, v)| v.to_owned())
        .for_each(|s| items_to_load.push(s));
    if bot.verbose {
        println!("Loading {} items", items_to_load.len());
    }
    bot.ec.load_entities(&bot.api, &items_to_load)?;
    Ok(())
}

pub fn load_basic_items(bot: &mut GeneDBot) -> Result<(), Box<Error>> {
    load_basic_items_chr(bot)?;
    load_basic_items_genes(bot)?;
    bot.evidence.load_from_wikidata(&mut bot.api)?;
    load_basic_items_entities(bot)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    /*
    TODO:
    fn process_gff_element(
    fn load_basic_items_genes(bot: &mut GeneDBot) -> Result<(), Box<Error>>
    fn load_basic_items(bot: &mut GeneDBot) -> Result<(), Box<Error>>
    */

    fn json_url() -> &'static str {
        TEST_URL_JSON
    }

    #[test]
    fn test_get_text_from_url() {
        assert_eq!(
            get_text_from_url(json_url()).unwrap(),
            "{\"test\":1234}".to_string()
        );
    }

    #[test]
    fn test_get_json_from_url() {
        assert_eq!(get_json_from_url(json_url()).unwrap(), json!({"test":1234}));
    }

    #[test]
    fn test_gff_url() {
        let mut bot = GeneDBot::new();
        bot.species_key = "test1234".to_string();
        assert_eq!(
            gff_url(&bot),
            "ftp://ftp.sanger.ac.uk/pub/genedb/releases/latest/test1234/test1234.gff.gz",
        );
        // NOTE: that URL does not exist, but we are only checking the building of it
    }

    #[test]
    fn test_gaf_url() {
        let mut bot = GeneDBot::new();
        bot.species_key = "test1234".to_string();
        assert_eq!(
            gaf_url(&bot),
            "ftp://ftp.sanger.ac.uk/pub/genedb/releases/latest/test1234/test1234.gaf.gz"
        );
        // NOTE: that URL does not exist, but we are only checking the building of it
    }

    #[test]
    fn test_sparql_result_to_pairs() {
        let mut api = mediawiki::api::Api::new("https://www.wikidata.org/w/api.php").unwrap();
        let j = json!({"head":{"vars":["q","genedb"]},"results":{"bindings":[{"q":{"type":"uri","value":"http://www.wikidata.org/entity/Q18968266"},"genedb":{"datatype":"http://www.w3.org/2001/XMLSchema#string","type":"literal","value":"PF3D7_0220200"}}]}});
        let expected = vec![("PF3D7_0220200".to_string(), "Q18968266".to_string())];
        let result = sparql_result_to_pairs(&mut api, &j["results"]["bindings"], "genedb", "q");
        assert_eq!(result, expected,)
    }

    #[test]
    fn test_find_genomic_assembly() {
        let mut bot = GeneDBot::new();
        bot.config.wikidata_id = "Q61779043".to_string();
        find_genomic_assembly(&mut bot, false).unwrap();
        assert_eq!(bot.genomic_assembly_q, "Q61815002");
    }

    #[test]
    fn test_create_genomic_assembly_item() {
        let mut bot = GeneDBot::new();
        bot.config.wikidata_id = "Q61779043".to_string();
        let item = create_genomic_assembly_item(&mut bot).unwrap();
        assert_eq!(
            item.label_in_locale("en"),
            Some("Plasmodium falciparum 3D7 reference genome")
        );
        assert!(item.has_target_entity("P279", "Q7307127"));
        assert!(item.has_target_entity("P703", "Q61779043"));
    }

    #[test]
    fn test_fix_id() {
        assert_eq!(fix_id("abc:def"), "abc");
        assert_eq!(fix_id("abc.9"), "abc");
        assert_eq!(fix_id("abc.123"), "abc.123");
    }

    #[test]
    fn test_get_gene_entities_to_process() {
        let mut bot = GeneDBot::new();
        bot.genedb2q.insert("foo".to_string(), "Q123".to_string());
        bot.protein_genedb2q
            .insert("bar".to_string(), "Q456".to_string());
        assert_eq!(get_gene_entities_to_process(&bot), vec!["Q123", "Q456"]);
        bot.specific_genes_only = Some(vec!["Q123".to_string()]);
        let empty_vec: Vec<String> = vec![];
        assert_eq!(get_gene_entities_to_process(&bot), empty_vec);
    }

    #[test]
    fn test_load_basic_items_entities() {
        let mut bot = GeneDBot::new();
        bot.genedb2q
            .insert("Count Count".to_string(), "Q12345".to_string());
        bot.protein_genedb2q
            .insert("Douglas Adams".to_string(), "Q42".to_string());
        bot.evidence
            .code2q
            .insert("Tim Berners-Lee".to_string(), "Q80".to_string());
        load_basic_items_entities(&mut bot).unwrap();
        assert!(bot.ec.has_entity("Q12345"));
        assert!(bot.ec.has_entity("Q42"));
        assert!(bot.ec.has_entity("Q80"));
        assert!(!bot.ec.has_entity("Q22686"));
    }

    #[test]
    fn test_load_basic_items_chr() {
        let mut bot = GeneDBot::new();
        bot.config.wikidata_id = "Q61779043".to_string();
        load_basic_items_chr(&mut bot).unwrap();
        assert!(bot.ec.has_entity(bot.config.wikidata_id)); // Force-loaded
        assert!(bot.ec.has_entity(TMHMM_Q)); // Force-loaded
        assert!(bot.ec.has_entity("Q61866468")); // Pf3D7_03_v3
        assert!(!bot.ec.has_entity("Q12345")); // Count Count
    }

    #[test]
    fn test_load_gaf_file_from_url() {
        let mut bot = GeneDBot::new();
        load_gaf_file_from_url(&mut bot, TEST_URL_GAF_GZ).unwrap();
        assert!(bot.gaf.contains_key("PF3D7_0100100.1"));
        assert_eq!(
            bot.gaf
                .get("PF3D7_0100100.1")
                .unwrap()
                .get(0)
                .unwrap()
                .go_id(),
            "GO:0009405"
        );
        // Just testing correct loading, not testing GAF parsing any further here
    }

    #[test]
    fn test_load_gff_file_from_url() {
        let mut bot = GeneDBot::new();
        load_gff_file_from_url(&mut bot, TEST_URL_GFF_GZ).unwrap();
        assert!(bot.gff.contains_key("Pfalciparum_REP_25"));
        assert_eq!(*bot.gff.get("Pfalciparum_REP_25").unwrap().start(), 9313);
        // Just testing correct loading, not testing GFF parsing any further here
    }

    #[test]
    fn test_set_other_types() {
        let mut bot = GeneDBot::new();
        bot.alternate_gene_subclasses
            .insert("polypeptide_motif".to_string(), "Q12345".to_string());
        load_gff_file_from_url(&mut bot, TEST_URL_GFF_GZ).unwrap();
        bot.gff
            .clone()
            .iter()
            .for_each(|element| match element.1.attributes().get("ID") {
                Some(id) => set_other_types(&mut bot, &element.1, id),
                None => {}
            });

        let other = bot.other_types.get("CDS").unwrap();
        assert_eq!(other.len(), 2);
        assert_eq!(other.iter().filter(|n| n.1.is_none()).count(), 2);
        assert_eq!(
            other
                .iter()
                .filter(|n| n.0 == "PF3D7_0100100" || n.0 == "PF3D7_0100200")
                .count(),
            2
        );

        let other = bot.other_types.get("polypeptide_motif").unwrap();
        assert_eq!(other.len(), 2);
        assert_eq!(
            other
                .iter()
                .filter(|n| *n.1.to_owned().unwrap().start() == 29568
                    || *n.1.to_owned().unwrap().start() == 39847)
                .count(),
            2
        );
        assert_eq!(
            other
                .iter()
                .filter(|n| n.0 == "chr01" || n.0 == "PF3D7_0100200")
                .count(),
            2
        );
    }

    // Wrappers, untested:

    #[test]
    fn test_load_gaf_file() {
        // Just a wrapper for gaf_url() and load_gaf_file_from_url()
    }

    #[test]
    fn test_load_gff_file() {
        // Just a wrapper for gff_url() and load_gff_file_from_url()
    }

    #[test]
    fn test_init() {
        // just a wrapper around other functions tested individually
    }

    #[test]
    fn create_genomic_assembly() {
        // just a wrapper around create_genomic_assembly_item and diff
    }

}
