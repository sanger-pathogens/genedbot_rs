use crate::{GeneDBot, Literature};
use regex::Regex;
use std::collections::{HashMap, HashSet};
use wikibase::entity_diff::*;
use wikibase::*;

pub fn process(
    bot: &mut GeneDBot,
    gene_genedb_id: &String,
    protein_genedb_id: &String,
) -> Option<String> {
    let gff = match bot.gff.get(protein_genedb_id) {
        Some(gff) => gff.clone(),
        None => return None,
    };

    let mut item = Entity::new_empty_item();
    let item_to_diff = match bot.get_entity_for_genedb_id(&protein_genedb_id) {
        Some(i) => i.clone(),
        None => Entity::new_empty_item(),
    };

    let reference = Reference::new(vec![
        Snak::new_item("P248", "Q5531047"),
        bot.new_time_today(),
    ]);

    let mut literature: HashSet<Literature> = HashSet::new();

    match gff.attributes().get("literature") {
        Some(lit) => lit.split(',').for_each(|l| {
            literature.insert(l.to_string());
        }),
        None => {}
    }

    item.set_label(LocaleString::new("en", &protein_genedb_id.clone()));

    let mut statements_to_create = vec![
        Snak::new_item("P31", "Q8054"),           // Instance of:Protein
        Snak::new_item("P703", &bot.species_q()), // Found in:Species
        Snak::new_string("P3382", &protein_genedb_id),
    ];

    match gff.feature_type() {
        "mRNA" => {
            statements_to_create.push(Snak::new_item("P279", "Q8054"));
        }
        "pseudogenic_transcript" => {
            statements_to_create.push(Snak::new_item("P279", "Q64698614"));
        }
        other => {
            bot.log(
                &protein_genedb_id,
                &format!("Protein has unknown type {}", other),
            );
        }
    }

    // Encoded by:gene
    match bot.get_entity_for_genedb_id(&gene_genedb_id) {
        Some(q) => statements_to_create.push(Snak::new_item("P702", q.id())),
        None => {}
    }

    // Create simple (single-reference-only) statements
    statements_to_create.iter().for_each(|s| {
        item.add_claim(Statement::new_normal(
            s.to_owned(),
            vec![],
            vec![reference.clone()],
        ))
    });

    bot.process_product(
        &gff,
        &mut item,
        &mut literature,
        &protein_genedb_id,
        &reference,
    );

    add_go_annotation(bot, &mut item, &gff, &mut literature);

    // Apply diff
    let my_props = vec![
        //"P1343", // Decribed by source; CAREFUL!
        "P703",  // Found in taxon
        "P1057", // Chromosome
        "P2548", // Strand orientation
        "P644",  // Genomic start
        "P645",  // Genomic end
        "P680", "P681", "P682",
    ];

    let mut params = EntityDiffParams::none();
    params.labels = EntityDiffParam::some(&vec!["en"]);
    params.descriptions = EntityDiffParam::some(&vec!["en"]);
    params.aliases = EntityDiffParam::some(&vec!["en"]);
    params.claims.add = EntityDiffParamState::All;
    params.claims.alter = EntityDiffParamState::All;
    params.claims.remove = EntityDiffParamState::some(&my_props);
    params.qualifiers = EntityDiffParamSub::all();
    params.references.list.push((
        EntityDiffParamState::some(&my_props),
        EntityDiffParamState::except(&vec!["P813"]),
    ));

    let mut diff = EntityDiff::new(&item_to_diff, &item, &params);
    diff.set_edit_summary(bot.get_edit_summary());
    protein_edit_validator(&mut diff, &item_to_diff);
    if !diff.is_empty() {
        if bot.verbose {
            println!(
                "\nPROTEIN {}/{:?}:\n{}",
                &protein_genedb_id,
                diff.edit_target(),
                serde_json::to_string(&diff.actions()).unwrap()
            );
        }
        if !bot.simulate {
            match bot.ec.apply_diff(&mut bot.api, &diff) {
                Some(q) => {
                    //thread::sleep(time::Duration::from_millis(500));
                    bot.genedb2q.insert(protein_genedb_id.to_string(), q);
                }
                None => bot.log(&protein_genedb_id, "Applying diff returned nothing"),
            }
        }
    }

    // TODO add main subject to literature

    bot.get_entity_id_for_genedb_id(&protein_genedb_id)
}

/// Blocks removal of claims [GO terms] that have a reference a non-GeneDB curator
fn protein_edit_validator(diff: &mut EntityDiff, original_item: &Entity) {
    let actions = diff.actions_mut();
    if actions["claims"].is_null() {
        return;
    }
    let claim_actions = match actions["claims"].as_array_mut() {
        Some(c) => c,
        None => return,
    };
    let props = ["P680", "P681", "P682"];

    // Closure should return false to remove the action
    claim_actions.retain(|action| {
        // Removals only
        if !action["remove"].is_string() {
            return true;
        }
        // Get the ID of the claim to be removed
        let id = match action["id"].as_str() {
            Some(id) => id,
            None => return true,
        };
        // Find that claim in the original item
        let claim = match original_item.claim_with_id(id) {
            Some(claim) => claim,
            None => return true,
        };
        // Is it a GO term property?
        if !props.contains(&claim.main_snak().property()) {
            return true;
        }
        let references = claim.references();
        // No reference? => remove
        if references.is_empty() {
            return false;
        }
        // More than one reference? => do not remove
        if references.len() > 1 {
            return true;
        }
        // Any P1640 (=curator) snaks that are NOT Q5531047 (GeneDB)?
        references
            .get(0)
            .unwrap()
            .snaks()
            .iter()
            .filter(|snak| snak.property() == "P1640")
            .filter(|snak| match snak.data_value() {
                Some(dv) => match dv.value() {
                    Value::Entity(value) => value.id() != "Q5531047",
                    _ => false,
                },
                None => false,
            })
            .count()
            > 0
    });
}

/// Adds GO annotation to the protein item
fn add_go_annotation(
    bot: &mut GeneDBot,
    item: &mut Entity,
    gff: &bio::io::gff::Record,
    literature: &mut HashSet<Literature>,
) {
    lazy_static! {
        static ref RE_DATE: Regex = Regex::new(r"^(\d{4})(\d{2})(\d{2})$").unwrap();
    }

    let protein_genedb_id = gff.attributes()["ID"].clone();
    let gaf = match bot.gaf.get(&protein_genedb_id) {
        Some(gaf) => gaf.clone(),
        None => return,
    };
    let mut new_go_claims: HashMap<String, (Snak, Vec<Reference>, Vec<Snak>)> = HashMap::new();
    for ga in gaf {
        let go_term = ga.go_id().to_string();
        let go_q = match bot.get_item_for_go_term(&go_term) {
            Some(q) => q,
            None => {
                bot.log(
                    &protein_genedb_id,
                    &format!("No Wikidata item for GO term '{}'", &go_term),
                );
                continue;
            }
        };

        let aspect = ga.aspect().to_string();
        let aspect_p = match bot.aspects.get(&aspect) {
            Some(p) => p.clone(),
            None => {
                bot.log(
                    &protein_genedb_id,
                    &format!("Unknown aspect '{}' for GO term '{}'", &aspect, &go_term),
                );
                continue;
            }
        };

        let evidence_code = ga.evidence_code().to_string();
        let evidence_code_q = match bot.evidence_codes.get(&evidence_code) {
            Some(q) => q.clone(),
            None => {
                bot.log(
                    &protein_genedb_id,
                    &format!(
                        "Unknown evidence code '{}' for GO term '{}'",
                        &evidence_code, &go_term
                    ),
                );
                continue;
            }
        };

        // Literature
        let mut literature_sources: Vec<Snak> = vec![];
        for (k, values) in ga.db_ref().iter_all() {
            for v in values {
                match k.as_str() {
                        "GO_REF" => {
                            literature_sources.push(Snak::new_string("P854", &format!("https://github.com/geneontology/go-site/blob/master/metadata/gorefs/goref-{}.md",&v)))
                        }
                        "PMID" => {
                            match bot.papers.get_or_create_item(k, v) {
                                Some(paper_q) => literature_sources.push(Snak::new_item("P248", &paper_q)),
                                None => bot.log(&protein_genedb_id,&format!("Can't find item for PMID '{}'",&v)),
                            }
                        }
                        other => {
                            bot.log(&protein_genedb_id,&format!("Unknown db_ref literature key '{}'",other));
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

                // Date
                match RE_DATE.captures(ga.date()) {
                    Some(caps) => {
                        let time = "+".to_string()
                            + caps.get(1).unwrap().as_str()
                            + "-"
                            + caps.get(2).unwrap().as_str()
                            + "-"
                            + caps.get(3).unwrap().as_str()
                            + "T00:00:00Z";
                        qualifiers.push(Snak::new_time("P585", &time, 11));
                    }
                    None => {}
                }

                // Qualifiers from with_from
                for w_f in ga.with_from() {
                    let parts: Vec<&str> = w_f.split(':').collect();
                    match bot.get_with_from_qualifier(&parts) {
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
                            bot.new_time_today(),
                        ])
                    })
                    .collect();

                let new_claim_key = json!([aspect_p.clone(), go_q.clone(), qualifiers.clone()]);
                let new_claim_key = serde_json::to_string(&new_claim_key).unwrap();

                if !new_go_claims.contains_key(&new_claim_key) {
                    new_go_claims.insert(
                        new_claim_key.clone(),
                        (
                            Snak::new_item(aspect_p.clone(), go_q.clone()),
                            vec![],
                            qualifiers,
                        ),
                    );
                }

                let ngc = new_go_claims.get_mut(&new_claim_key).unwrap(); // This was just added
                for r in references {
                    ngc.1.push(r);
                }
                deduplicate_references(&mut ngc.1);

                for (k, values) in ga.db_ref().iter_all() {
                    for v in values {
                        literature.insert(format!("{}:{}", k, v).to_string());
                    }
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
                        item.add_alias(LocaleString::new("en", &bot.fix_alias_name(s)));
                    });
                });
            }
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

fn deduplicate_references(references: &mut Vec<Reference>) {
    let mut ref_keys: HashSet<String> = HashSet::new();
    references.retain(|r| {
        let ref_key = serde_json::to_string(&r).unwrap();
        if ref_keys.contains(&ref_key) {
            false
        } else {
            ref_keys.insert(ref_key);
            true
        }
    });
}
