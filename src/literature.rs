use papers::crossref2wikidata::Crossref2Wikidata;
use papers::orcid2wikidata::Orcid2Wikidata;
use papers::pubmed2wikidata::Pubmed2Wikidata;
use papers::semanticscholar2wikidata::Semanticscholar2Wikidata;
use papers::wikidata_papers::WikidataPapers;
use papers::wikidata_string_cache::WikidataStringCache;
use papers::*;
use std::collections::HashMap;
use std::sync::{Arc, RwLock};

trait Wikibase {
    fn api(self: &mut Self) -> &mut wikibase::mediawiki::api::Api;

    fn search_wikibase(&mut self, query: &String) -> Vec<String> {
        let params: HashMap<_, _> = vec![
            ("action", "query"),
            ("list", "search"),
            ("srnamespace", "0"),
            ("srsearch", &query.as_str()),
        ]
        .into_iter()
        .map(|(x, y)| (x.to_string(), y.to_string()))
        .collect();
        let res = self.api().get_query_api_json(&params).unwrap();
        match res["query"]["search"].as_array() {
            Some(items) => items
                .iter()
                .map(|item| item["title"].as_str().unwrap().to_string())
                .collect(),
            None => vec![],
        }
    }
}

#[derive(Debug, Clone)]
pub struct Papers {
    paper2q: HashMap<(String, String), String>,
    api: wikibase::mediawiki::api::Api,
}

impl Wikibase for Papers {
    fn api(self: &mut Self) -> &mut wikibase::mediawiki::api::Api {
        &mut self.api
    }
}

impl Papers {
    pub fn new(api: &wikibase::mediawiki::api::Api) -> Self {
        Papers {
            paper2q: HashMap::new(),
            api: api.clone(),
        }
    }

    pub fn get_or_create_item(&mut self, k: &String, v: &String) -> Option<String> {
        if v == "PMID" {
            return None;
        }
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
                /*
                let sparql = format!("SELECT ?q {{ ?q wdt:P698 '{}' }}", &v);
                let sparql_result = self.api.sparql_query(&sparql).ok()?;
                let items = self.api.entities_from_sparql_result(&sparql_result, "q");
                */
                let items = self.search_wikibase(&format!("haswbstatement:P698={}", &v));
                match items.len() {
                    0 => {
                        let mut ids = vec![GenericWorkIdentifier::new_prop(PROP_PMID, v)];
                        match self.create_paper_item(&mut ids) {
                            Some(q) => {
                                println!(
                                    "CREATED NEW PAPER ITEM https://www.wikidata.org/wiki/{}",
                                    &q
                                );
                                self.paper2q
                                    .insert((k.to_string(), v.to_string()), q.clone());
                                Some(q.clone())
                            }
                            None => {
                                println!("FAILED TO CREATE WIKIDATA ITEM FOR PMID:{}", &v);
                                None
                            }
                        }
                    }
                    1 => {
                        self.paper2q
                            .insert((k.to_string(), v.to_string()), items[0].clone());
                        Some(items[0].clone())
                    }
                    _ => {
                        // Multiple, use first one
                        self.paper2q
                            .insert((k.to_string(), v.to_string()), items[0].clone());
                        Some(items[0].clone())
                    }
                }
            }
            other => {
                println!("Unknown paper source: '{}'", &other);
                None
            }
        }
    }

    fn create_paper_item(&mut self, ids: &mut Vec<GenericWorkIdentifier>) -> Option<String> {
        println!("Creating paper item for {:?}", &ids);
        let mw_api = Arc::new(RwLock::new(self.api.clone()));
        let cache = Arc::new(WikidataStringCache::new(mw_api.clone()));
        let mut wdp = WikidataPapers::new(cache);
        wdp.add_adapter(Box::new(Pubmed2Wikidata::new()));
        wdp.add_adapter(Box::new(Crossref2Wikidata::new()));
        wdp.add_adapter(Box::new(Semanticscholar2Wikidata::new()));
        wdp.add_adapter(Box::new(Orcid2Wikidata::new()));
        let ids = wdp.update_from_paper_ids(&ids);
        match wdp.create_or_update_item_from_ids(mw_api, &ids) {
            Some(edit_result) => Some(edit_result.q),
            None => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new() {
        let api = wikibase::mediawiki::api::Api::new("https://www.wikidata.org/w/api.php").unwrap();
        let papers = Papers::new(&api);
        assert!(papers.paper2q.is_empty());
    }

    #[test]
    fn test_api() {
        let api = wikibase::mediawiki::api::Api::new("https://www.wikidata.org/w/api.php").unwrap();
        let mut papers = Papers::new(&api);
        assert_eq!(
            papers
                .api()
                .get_site_info_string("general", "sitename")
                .unwrap(),
            "Wikidata"
        );
    }

    #[test]
    fn test_search_wikibase() {
        let api = wikibase::mediawiki::api::Api::new("https://www.wikidata.org/w/api.php").unwrap();
        let mut papers = Papers::new(&api);
        let result =
            papers.search_wikibase(&"\"Charles Darwin\" haswbstatement:P31=Q5".to_string());
        assert!(result.iter().any(|s| s == "Q1035"));
        assert!(!result.iter().any(|s| s == "Q1064071"));
    }

    #[test]
    fn test_get_or_create_item() {
        let api = wikibase::mediawiki::api::Api::new("https://www.wikidata.org/w/api.php").unwrap();
        let mut papers = Papers::new(&api);
        assert_eq!(
            papers.get_or_create_item(&"PMID".to_string(), &"27998271".to_string()),
            Some("Q28030910".to_string())
        );
        assert_eq!(
            papers.get_or_create_item(&"PMID".to_string(), &"0".to_string()),
            None
        );
    }

    #[test]
    fn test_create_paper_item() {
        let api = wikibase::mediawiki::api::Api::new("https://www.wikidata.org/w/api.php").unwrap();
        let mut papers = Papers::new(&api);
        assert_eq!(papers.create_paper_item(&mut vec![]), None)
    }
}
