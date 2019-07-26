use crate::{GeneDBot, Literature, Toolbox};
use std::collections::HashSet;
use wikibase::entity_diff::*;
use wikibase::*;

pub fn process(bot: &mut GeneDBot, genedb_id: String) {
    let gff = match bot.gff.get(&genedb_id) {
        Some(gff) => gff.clone(),
        None => return,
    };
    let gene_type = match gff.feature_type() {
        "gene" => ("gene", "Q7187"),
        "pseudogene" => ("pseudogene", "Q277338"),
        other => {
            bot.log(&genedb_id, &format!("Gene has unknown type {}", other));
            return;
        }
    };

    let mut item = Entity::new_empty_item();
    let item_to_diff = match bot.get_entity_for_genedb_id(&genedb_id) {
        Some(i) => i.clone(),
        None => Entity::new_empty_item(),
    };
    let chr_q = match bot.get_or_create_chromosome_entity(gff.seqname()) {
        Some(chr_q) => chr_q.to_owned(),
        None => {
            bot.log(
                &genedb_id,
                &format!("Could not create chromosome item for '{}'", gff.seqname()),
            );
            return;
        }
    };

    // Labels and aliases
    match gff.attributes().get("Name") {
        Some(name) => {
            item.set_label(LocaleString::new("en", name));
            //item.add_alias(LocaleString::new("en", &bot.fix_alias_name(&genedb_id)));
        }
        None => item.set_label(LocaleString::new("en", &genedb_id)),
    };
    vec!["previous_systematic_id", "synonym", "alias"]
        .iter()
        .for_each(|key| {
            match gff.attributes().get(&key.to_string()) {
                Some(ids) => ids.split(',').for_each(|id| {
                    item.add_alias(LocaleString::new("en", &bot.fix_alias_name(id)))
                }),
                None => {}
            };
        });

    // Statements
    let reference = Reference::new(vec![
        Snak::new_item("P248", "Q5531047"),
        bot.new_time_today(),
    ]);
    let ga_quals = vec![
        Snak::new_item("P659", &bot.genomic_assembly_q),
        Snak::new_item("P1057", &chr_q),
    ];

    let mut statements_to_create = vec![
        Snak::new_item("P31", gene_type.1),       // Instance of
        Snak::new_item("P703", &bot.species_q()), // Found in:Species
        Snak::new_item("P1057", &chr_q),          // Chromosome
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

    let protein_entity_ids = bot.process_proteins(&genedb_id);
    if protein_entity_ids.len() > 0 {
        if gene_type.0 == "gene" {
            // Genes only, no pseudogene
            statements_to_create.push(Snak::new_item("P279", "Q20747295")); // Subclass of:protein-coding gene
        } else if gene_type.0 == "pseudogene" {
            statements_to_create.push(Snak::new_item("P279", "Q277338")); // Subclass of:pseudogene
        }

        // Encodes: protein
        for protein_q in &protein_entity_ids {
            statements_to_create.push(Snak::new_item("P688", &protein_q));
        }
    } else if gene_type.0 != "gene" {
        statements_to_create.push(Snak::new_item("P279", gene_type.1));
    } else {
        let mut subclass_found: bool = false;
        for subclass in bot.alternate_gene_subclasses.clone() {
            let class_name = subclass.0;
            let class_q = subclass.1;
            let ct = match bot.other_types.get(&class_name.to_owned()) {
                Some(ct) => ct,
                None => continue,
            };
            match ct.get(&genedb_id) {
                Some(gff_tmp_opt) => match gff_tmp_opt.clone() {
                    Some(gff_tmp) => {
                        statements_to_create.push(Snak::new_item("P279", &class_q));
                        subclass_found = true;
                        let mut fake_literature: HashSet<Literature> = HashSet::new();
                        bot.process_product(
                            &gff_tmp,
                            &mut item,
                            &mut fake_literature,
                            &genedb_id,
                            &reference,
                        );
                    }
                    None => {}
                },
                None => {}
            };
        }
        if !subclass_found {
            bot.log(
                &genedb_id,
                &format!("No subclass found for {:?}", &gene_type),
            );
        }
    }

    // Orthologs
    if !protein_entity_ids.is_empty() {
        match bot.parent2child.get(&genedb_id) {
            Some(child) => {
                bot.orthologs
                    .process(&child, &bot.gff, &mut item, &reference);
            }
            None => {}
        }
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
        "P688",  // Encodes
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

    let mut diff = EntityDiff::new(&item_to_diff, &item, &params);
    diff.set_edit_summary(bot.get_edit_summary());
    if !diff.is_empty() {
        if bot.verbose {
            println!(
                "\nGENE {}/{:?}:\n{}",
                &genedb_id,
                diff.edit_target(),
                serde_json::to_string(&diff.actions()).unwrap()
            );
        }
        if !bot.simulate {
            match bot.ec.apply_diff(&mut bot.api, &diff) {
                Some(q) => {
                    //thread::sleep(time::Duration::from_millis(500));
                    bot.genedb2q.insert(genedb_id.to_string(), q);
                }
                None => bot.log(&genedb_id, "Applying diff returned nothing"),
            }
        }
    }

    match bot.get_entity_id_for_genedb_id(&genedb_id) {
        Some(gene_q) => {
            for protein_q in protein_entity_ids {
                link_protein_to_gene(bot, &protein_q, &gene_q);
            }
        }
        None => {}
    }
}

fn link_protein_to_gene(bot: &mut GeneDBot, protein_q: &String, gene_q: &String) {
    if !bot.is_item(gene_q) || !bot.is_item(protein_q) {
        return;
    }
    match bot
        .ec
        .load_entities(&bot.api, &vec![gene_q.clone(), protein_q.clone()])
    {
        Ok(_) => {}
        _ => return,
    }
    let gene_i = match bot.ec.get_entity(gene_q.as_str()) {
        Some(i) => i.clone(),
        None => return,
    };
    let protein_i = match bot.ec.get_entity(protein_q.as_str()) {
        Some(i) => i.clone(),
        None => return,
    };
    link_items(bot, "P688", &gene_i, protein_q.to_string());
    link_items(bot, "P702", &protein_i, gene_q.to_string());
}

fn link_items(bot: &mut GeneDBot, property: &str, item: &Entity, target_q: String) {
    if item.has_target_entity(property, &target_q) {
        return;
    }
    let mut new_item = item.clone();
    new_item.add_claim(Statement::new_normal(
        Snak::new_item(property, &target_q),
        vec![],
        bot.references(),
    ));
    let params = EntityDiffParams::all();
    let mut diff = EntityDiff::new(&item, &new_item, &params);
    diff.set_edit_summary(bot.get_edit_summary());
    /*
    println!(
        "{} ={}=> {} : {}",
        item.id(),
        &property,
        &target_q,
        serde_json::to_string(&diff.actions()).unwrap()
    );
    */
    if !bot.simulate {
        // Run, but ignore result
        match bot.ec.apply_diff(&mut bot.api, &diff) {
            Some(_) => {}
            None => {}
        }
    }
}

#[cfg(test)]
mod tests {
    //use super::*;
}
