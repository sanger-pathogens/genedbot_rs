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
        api: &mut wikibase::mediawiki::api::Api,
    ) -> Result<(), Box<dyn Error>> {
        let sparql = "SELECT DISTINCT ?q ?qLabel ?qAltLabel { ?q wdt:P31 wd:Q23173209 SERVICE wikibase:label { bd:serviceParam wikibase:language 'en' } }" ;
        let res = api.sparql_query(&sparql)?;
        match res["results"]["bindings"].as_array() {
            Some(bindings) => {
                bindings
                    .iter()
                    .filter(|b| b["q"]["value"].as_str().is_some())
                    .filter(|b| b["qLabel"]["value"].as_str().is_some())
                    .for_each(|b| {
                        let q = b["q"]["value"].as_str().unwrap();
                        let q = api.extract_entity_from_uri(q).unwrap().to_string();
                        let label = b["qLabel"]["value"].as_str().unwrap_or("").to_string();
                        let alt_label = b["qAltLabel"]["value"].as_str().unwrap_or("").to_string();
                        self.label2q
                            .insert(self.normalize_label(&alt_label), q.clone());
                        self.code2q.insert(label, q.clone());
                    });
                Ok(())
            }
            None => Err(From::from(format!("load_from_wikidata failed"))),
        }
    }

    pub fn normalize_label(&self, s: &String) -> String {
        s.trim().to_lowercase()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let e = Evidence::new();
        assert!(e.label2q.is_empty());
        assert!(e.code2q.is_empty());
    }

    #[test]
    fn test_normalize_label() {
        let e = Evidence::new();
        assert_eq!(
            e.normalize_label(&" ThIs iS a tEst  ".to_string()),
            "this is a test"
        );
    }

    #[test]
    fn test_load_from_wikidata() {
        let mut api =
            wikibase::mediawiki::api::Api::new("https://www.wikidata.org/w/api.php").unwrap();
        let mut e = Evidence::new();
        e.load_from_wikidata(&mut api).unwrap();
        assert_eq!(
            e.label2q.get("inferred by curator"),
            Some(&"Q23190856".to_string())
        );
        assert_eq!(e.code2q.get("ISA"), Some(&"Q23190738".to_string()));
    }

}
