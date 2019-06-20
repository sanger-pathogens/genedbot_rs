use std::collections::{HashMap, HashSet};
//use std::{error::Error, fmt};

#[derive(Debug, Clone)]
pub struct Orthologs {
    pub genedb2q: HashMap<String, String>,
    pub genedb2taxon_q: HashMap<String, String>,
}

impl crate::Toolbox for Orthologs {}

impl Orthologs {
    pub fn new() -> Self {
        Self {
            genedb2q: HashMap::new(),
            genedb2taxon_q: HashMap::new(),
        }
    }

    pub fn load(
        self: &mut Self,
        api: &mediawiki::api::Api,
        orth_ids: HashSet<String>,
    ) -> Result<(), Box<::std::error::Error>> {
        if orth_ids.is_empty() {
            return Ok(());
        }

        if orth_ids.len() < 1000 {
            self.load_small_group(api, orth_ids)
        } else {
            self.load_large_group(api, orth_ids)
        }
    }

    fn load_small_group(
        self: &mut Self,
        api: &mediawiki::api::Api,
        mut orth_ids: HashSet<String>,
    ) -> Result<(), Box<::std::error::Error>> {
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
    ) -> Result<(), Box<::std::error::Error>> {
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
