use crate::{GeneDBot, TMHMM_Q};
use bio::io::{gaf, gff};
use libflate::gzip::Decoder;
use regex::Regex;
use std::collections::{HashMap, HashSet};
use std::error::Error;
use wikibase::entity_diff::*;
use wikibase::*;

pub fn init(bot: &mut GeneDBot) -> Result<(), Box<Error>> {
    println!("1");
    load_gff_file(bot)?; //.expect(&format!("Can't load GFF file '{}'", gff_url(bot)));
    println!("2");
    load_gaf_file(bot)?; //.expect(&format!("Can't load GAF file '{}'", gaf_url(bot)));
    println!("3");
    find_genomic_assembly(bot)?;
    println!("4");
    load_basic_items(bot)?;
    println!("5");
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
        panic!("Can't get GFF data from {}", url);
    }

    load_orthologs(bot, orth_ids)?;
    Ok(())
}

fn load_orthologs(
    bot: &mut GeneDBot,
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
            let sparql_result = bot.api.sparql_query(&sparql)?;
            for b in sparql_result["results"]["bindings"].as_array().unwrap() {
                let q = match b["q"]["value"].as_str() {
                    Some(s) => bot.api.extract_entity_from_uri(s).unwrap(),
                    None => continue,
                };
                let taxon_q = match b["taxon"]["value"].as_str() {
                    Some(s) => bot.api.extract_entity_from_uri(s).unwrap(),
                    None => continue,
                };
                let genedb = match b["genedb"]["value"].as_str() {
                    Some(s) => s.to_string(),
                    None => continue,
                };
                bot.orth_genedb2q.insert(genedb.to_string(), q);
                bot.orth_genedb2q_taxon.insert(genedb.to_string(), taxon_q);
            }
        }
    } else {
        // Retrieven 'em all and let HashSet sort 'em out...
        let sparql = "SELECT ?q ?genedb ?taxon { ?q wdt:P3382 ?genedb ; wdt:P703 ?taxon }";
        let sparql_result = bot.api.sparql_query(&sparql)?;
        for b in sparql_result["results"]["bindings"].as_array().unwrap() {
            let genedb = match b["genedb"]["value"].as_str() {
                Some(s) => s.to_string(),
                None => continue,
            };
            if !orth_ids.contains(&genedb) {
                continue;
            }
            let q = match b["q"]["value"].as_str() {
                Some(s) => bot.api.extract_entity_from_uri(s).unwrap(),
                None => continue,
            };
            let taxon_q = match b["taxon"]["value"].as_str() {
                Some(s) => bot.api.extract_entity_from_uri(s).unwrap(),
                None => continue,
            };
            bot.orth_genedb2q.insert(genedb.to_string(), q);
            bot.orth_genedb2q_taxon.insert(genedb.to_string(), taxon_q);
        }
    }
    Ok(())
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

    bot.get_orthologs_from_gff_element(&element)
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
        panic!("Can't get GAF data from {}", url);
    }
    Ok(())
}

pub fn load_gff_file(bot: &mut GeneDBot) -> Result<(), Box<Error>> {
    load_gff_file_from_url(bot, gff_url(bot).as_str())
}

pub fn load_gaf_file(bot: &mut GeneDBot) -> Result<(), Box<Error>> {
    load_gaf_file_from_url(bot, gaf_url(bot).as_str())
}

fn create_genomic_assembly(bot: &mut GeneDBot) -> Result<String, Box<Error>> {
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

    let params = EntityDiffParams::all();
    let mut diff = EntityDiff::new(&Entity::new_empty_item(), &new_item, &params);
    diff.set_edit_summary(bot.get_edit_summary());
    match bot.ec.apply_diff(&mut bot.api, &diff) {
        Some(q) => Ok(q),
        None => Err(From::from("Could not create genomic assembly item")),
    }
}

fn find_genomic_assembly(bot: &mut GeneDBot) -> Result<(), Box<Error>> {
    let sparql = format!(
        "SELECT ?q {{ ?q wdt:P279 wd:Q7307127 ; wdt:P703 wd:{} }}",
        bot.species_q()
    );
    let res = bot.api.sparql_query(&sparql)?;
    let candidates = bot.api.entities_from_sparql_result(&res, "q");
    bot.genomic_assembly_q = match candidates.len() {
        0 => create_genomic_assembly(bot)?,
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
    bot: &mut GeneDBot,
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
            let v1 = bot
                .api
                .extract_entity_from_uri(v1)
                .unwrap_or(v1.to_string())
                .to_string();
            let v2 = b[k2]["value"].as_str().unwrap();
            let v2 = bot
                .api
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
    bot.genedb2q = sparql_result_to_pairs(bot, &res["results"]["bindings"], "genedb", "q")
        .into_iter()
        .collect();

    // Proteins
    let sparql = format!("SELECT DISTINCT ?q ?genedb {{ {} . ?q wdt:P31 wd:Q8054 ; wdt:P703 ?species ; wdt:P3382 ?genedb }}",&species_list) ;
    let res = bot.api.sparql_query(&sparql)?;
    bot.protein_genedb2q = sparql_result_to_pairs(bot, &res["results"]["bindings"], "genedb", "q")
        .into_iter()
        .collect();
    Ok(())
}

pub fn load_evidence_codes(bot: &mut GeneDBot) -> Result<(), Box<Error>> {
    let sparql = "SELECT DISTINCT ?q ?qLabel ?qAltLabel { ?q wdt:P31 wd:Q23173209 SERVICE wikibase:label { bd:serviceParam wikibase:language \"en\" } }" ;
    let res = bot.api.sparql_query(&sparql)?;
    res["results"]["bindings"]
        .as_array()
        .unwrap()
        .iter()
        .filter(|b| b["q"]["value"].as_str().is_some())
        .filter(|b| b["qLabel"]["value"].as_str().is_some())
        .for_each(|b| {
            let q = b["q"]["value"].as_str().unwrap();
            let q = bot.api.extract_entity_from_uri(q).unwrap().to_string();
            let label = b["qLabel"]["value"].as_str().unwrap_or("").to_string();
            let alt_label = b["qAltLabel"]["value"].as_str().unwrap_or("").to_string();
            bot.evidence_codes_labels
                .insert(bot.normalize_evidence_label(&alt_label), q.clone());
            bot.evidence_codes.insert(label, q.clone());
        });
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
    bot.evidence_codes
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
    load_evidence_codes(bot)?;
    load_basic_items_entities(bot)?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test1() {
        let bot = GeneDBot::new();
        println!("{:?}", &bot.genes2load);
    }
}
