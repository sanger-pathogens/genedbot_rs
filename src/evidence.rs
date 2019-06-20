use std::collections::HashMap;
use std::error::Error;

#[derive(Debug, Clone)]
pub struct Evidence {
    pub label2q: HashMap<String, String>,
    pub code2q: HashMap<String, String>,
}

impl Evidence {
    pub fn new() -> Self {
        Self {
            label2q: HashMap::new(),
            code2q: HashMap::new(),
        }
    }

    pub fn load_from_wikidata(
        self: &mut Self,
        api: &mut mediawiki::api::Api,
    ) -> Result<(), Box<Error>> {
        let sparql = "SELECT DISTINCT ?q ?qLabel ?qAltLabel { ?q wdt:P31 wd:Q23173209 SERVICE wikibase:label { bd:serviceParam wikibase:language 'en' } }" ;
        let res = api.sparql_query(&sparql)?;
        res["results"]["bindings"]
            .as_array()
            .unwrap()
            .iter()
            .filter(|b| b["q"]["value"].as_str().is_some())
            .filter(|b| b["qLabel"]["value"].as_str().is_some())
            .for_each(|b| {
                let q = b["q"]["value"].as_str().unwrap();
                let q = api.extract_entity_from_uri(q).unwrap().to_string();
                let label = b["qLabel"]["value"].as_str().unwrap_or("").to_string();
                let alt_label = b["qAltLabel"]["value"].as_str().unwrap_or("").to_string();
                self.label2q
                    .insert(self.normalize_evidence_label(&alt_label), q.clone());
                self.code2q.insert(label, q.clone());
            });
        Ok(())
    }

    pub fn normalize_evidence_label(&self, s: &String) -> String {
        s.trim().to_lowercase()
    }
}
