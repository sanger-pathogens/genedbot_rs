use crate::Toolbox;
use regex::Regex;
use std::collections::{HashMap, HashSet};
use wikibase::*;

const MIN_LARGE_GROUP: usize = 100;

#[derive(Debug, Clone)]
pub struct Orthologs {
    pub genedb2q: HashMap<String, String>,
    pub genedb2taxon_q: HashMap<String, String>,
}

impl Toolbox for Orthologs {}

impl Orthologs {
    pub fn new() -> Self {
        Self {
            genedb2q: HashMap::new(),
            genedb2taxon_q: HashMap::new(),
        }
    }

    // Returns (species,genedb_id)
    pub fn get_from_gff_element(&self, gff: &bio::io::gff::Record) -> Vec<(String, String)> {
        lazy_static! {
            static ref RE_ORTH: Regex =
                Regex::new(r"^\s*(\S*):(\S+)").expect("RE_ORTH does not compile");
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
                        .map(|o| o.expect("get_orthologs_from_gff_element: [1]"))
                        .collect();
                    ret.append(&mut result);
                }
                ret
            }
            None => vec![],
        }
    }

    pub fn process(
        &mut self,
        child: &Vec<(String, String)>,
        gff: &HashMap<String, bio::io::gff::Record>,
        item: &mut Entity,
        reference: &Reference,
    ) {
        let mut had_that: HashSet<String> = HashSet::new();
        child
            .iter()
            .filter(|o| self.is_product_type(&o.1))
            .for_each(|o| match gff.get(&o.0) {
                Some(protein) => {
                    self.get_from_gff_element(&protein).iter().for_each(
                        |(_species, protein_genedb_id)| match self.genedb2q.get(protein_genedb_id) {
                            Some(orth_q) => {
                                if had_that.contains(orth_q) {
                                    return;
                                }
                                had_that.insert(orth_q.to_string());
                                match self.genedb2taxon_q.get(protein_genedb_id) {
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
                        },
                    );
                }
                None => {}
            });
    }

    pub fn load(
        self: &mut Self,
        api: &mediawiki::api::Api,
        orth_ids: HashSet<String>,
    ) -> Result<(), Box<dyn (::std::error::Error)>> {
        if orth_ids.is_empty() {
            return Ok(());
        }

        if orth_ids.len() < MIN_LARGE_GROUP {
            self.load_small_group(api, orth_ids)
        } else {
            self.load_large_group(api, orth_ids)
        }
    }

    fn load_small_group(
        self: &mut Self,
        api: &mediawiki::api::Api,
        mut orth_ids: HashSet<String>,
    ) -> Result<(), Box<dyn (::std::error::Error)>> {
        // Usually for testing
        let orth_ids: Vec<String> = orth_ids.drain().collect();
        for chunk in orth_ids.chunks(100) {
            let sparql = format!("SELECT ?q ?genedb ?taxon {{ VALUES ?genedb {{'{}'}} . ?q wdt:P3382 ?genedb ; wdt:P703 ?taxon }}",chunk.join("' '"));
            let sparql_result = api.sparql_query(&sparql)?;
            for b in sparql_result["results"]["bindings"].as_array().unwrap() {
                let q = match b["q"]["value"].as_str() {
                    Some(s) => api.extract_entity_from_uri(s).unwrap(),
                    None => continue,
                };
                let taxon_q = match b["taxon"]["value"].as_str() {
                    Some(s) => api.extract_entity_from_uri(s).unwrap(),
                    None => continue,
                };
                let genedb = match b["genedb"]["value"].as_str() {
                    Some(s) => s.to_string(),
                    None => continue,
                };
                self.genedb2q.insert(genedb.to_string(), q);
                self.genedb2taxon_q.insert(genedb.to_string(), taxon_q);
            }
        }
        Ok(())
    }

    fn load_large_group(
        self: &mut Self,
        api: &mediawiki::api::Api,
        orth_ids: HashSet<String>,
    ) -> Result<(), Box<dyn (::std::error::Error)>> {
        // Retrieven 'em all and let HashSet sort 'em out...
        let sparql = "SELECT ?q ?genedb ?taxon { ?q wdt:P3382 ?genedb ; wdt:P703 ?taxon }";
        let sparql_result = api.sparql_query(&sparql)?;
        for b in sparql_result["results"]["bindings"].as_array().unwrap() {
            let genedb = match b["genedb"]["value"].as_str() {
                Some(s) => s.to_string(),
                None => continue,
            };
            if !orth_ids.contains(&genedb) {
                continue;
            }
            let q = match b["q"]["value"].as_str() {
                Some(s) => api.extract_entity_from_uri(s).unwrap(),
                None => continue,
            };
            let taxon_q = match b["taxon"]["value"].as_str() {
                Some(s) => api.extract_entity_from_uri(s).unwrap(),
                None => continue,
            };
            self.genedb2q.insert(genedb.to_string(), q);
            self.genedb2taxon_q.insert(genedb.to_string(), taxon_q);
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::loader::{load_gff_file_from_url, TEST_URL_GFF_GZ};

    #[test]
    fn test_new() {
        let o = Orthologs::new();
        assert_eq!(o.genedb2q.len(), 0);
        assert_eq!(o.genedb2taxon_q.len(), 0);
    }

    #[test]
    fn test_get_from_gff_element() {
        let o = Orthologs::new();
        let mut bot = crate::genedbot::GeneDBot::new();
        load_gff_file_from_url(&mut bot, TEST_URL_GFF_GZ).unwrap();
        let gff_element = bot.gff.get("PF3D7_0100200.1").unwrap();
        let result = o.get_from_gff_element(&gff_element);
        assert!(result
            .iter()
            .any(|r| r.0 == "Preichenowi" && r.1 == "PRCDC_0042600"));
    }

    #[test]
    fn test_process() {
        let mut o = Orthologs::new();
        let mut bot = crate::genedbot::GeneDBot::new();
        load_gff_file_from_url(&mut bot, TEST_URL_GFF_GZ).unwrap();
        o.genedb2q
            .insert("PRCDC_0042600".to_string(), "Q123".to_string());
        o.genedb2taxon_q
            .insert("PRCDC_0042600".to_string(), "Q456".to_string());
        let child = vec![("PF3D7_0100200.1".to_string(), "mRNA".to_string())];
        let mut item = Entity::new_empty_item();
        let reference = Reference::new(vec![Snak::new_string("P214", "test")]);
        let _result = o.process(&child, &bot.gff, &mut item, &reference);
        assert!(item.has_target_entity("P684", "Q123"));
        assert_eq!(
            *item.claims().get(0).unwrap().qualifiers().get(0).unwrap(),
            Snak::new_item("P703", "Q456")
        );
        assert_eq!(
            *item
                .claims()
                .get(0)
                .unwrap()
                .references()
                .get(0)
                .unwrap()
                .snaks()
                .get(0)
                .unwrap(),
            Snak::new_string("P214", "test")
        );
    }

    #[test]
    fn test_load_small() {
        let mut o = Orthologs::new();
        let mut api = mediawiki::api::Api::new("https://www.wikidata.org/w/api.php").unwrap();
        let mut orth_ids: HashSet<String> = HashSet::new();
        orth_ids.insert("PF3D7_0102600".to_string());
        o.load(&mut api, orth_ids).unwrap();
        assert_eq!(o.genedb2q.len(), 1);
        assert_eq!(o.genedb2taxon_q.len(), 1);
        assert_eq!(o.genedb2q.get("PF3D7_0102600").unwrap(), "Q18968367");
        assert_eq!(o.genedb2taxon_q.get("PF3D7_0102600").unwrap(), "Q61779043");
    }

    #[test]
    fn test_load_large() {
        let mut o = Orthologs::new();
        let mut api = mediawiki::api::Api::new("https://www.wikidata.org/w/api.php").unwrap();
        let mut orth_ids: HashSet<String> = HashSet::new();
        orth_ids.insert("PF3D7_0102600".to_string());
        for x in 0..MIN_LARGE_GROUP + 5 {
            orth_ids.insert("this does not exist".to_string() + &format!("{}", x));
        }
        o.load(&mut api, orth_ids).unwrap();
        assert_eq!(o.genedb2q.len(), 1);
        assert_eq!(o.genedb2taxon_q.len(), 1);
        assert_eq!(o.genedb2q.get("PF3D7_0102600").unwrap(), "Q18968367");
        assert_eq!(o.genedb2taxon_q.get("PF3D7_0102600").unwrap(), "Q61779043");
    }
}
