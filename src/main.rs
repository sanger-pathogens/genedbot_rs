use std::env;
#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate serde_json;
extern crate chrono;
extern crate clap;
extern crate config;
extern crate libflate;
extern crate papers;
extern crate regex;
extern crate reqwest;

use crate::genedbot::*;
use clap::{App, Arg};
use config::{Config, File};
use std::error::Error;

pub mod evidence;
pub mod gene;
pub mod genedbot;
pub mod literature;
pub mod loader;
pub mod orthologs;
pub mod protein;

fn run_bot_for_species_and_gene(
    species_key: &String,
    genes: &Option<Vec<String>>,
    lgname: &str,
    lgpass: &str,
) -> Result<(), Box<dyn Error>> {
    let mut bot = GeneDBot::new();
    //bot.simulate = true;
    //bot.verbose = true;
    bot.api().write().unwrap().set_user_agent("GeneDBot/3.0");
    bot.api().write().unwrap().set_edit_delay(Some(500)); // Half a second between edits
    bot.specific_genes_only = genes.to_owned();
    bot.api().write().unwrap().login(lgname, lgpass)?;
    bot.load_config_file(species_key)?;
    bot.init()?;
    bot.run()?;
    Ok(())
}

fn main() {
    let matches = App::new("GeneDBot")
        .version("0.1")
        .author("Magnus Manske <mm6@sanger.ac.uk>")
        .about("Updates Wikidata from CHADO GFF/GAF files")
        .arg(
            Arg::with_name("config")
                .short("c")
                .long("config")
                .value_name("FILE")
                .required(false)
                .help("Sets a custom config file")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("gene")
                .short("g")
                .long("gene")
                .required(false)
                .help("Only process specific gene(s)")
                .takes_value(true),
        )
        .arg(
            Arg::with_name("SPECIES_KEY")
                .help("Species key, or 'all'")
                .required(true)
                .index(1),
        )
        .get_matches();

    let ini_file = matches.value_of("config").unwrap_or("bot.ini");
    let mut settings = Config::default();
    settings.merge(File::with_name(ini_file)).unwrap();
    let lgname = settings.get_str("user.user").unwrap();
    let lgpass = settings.get_str("user.pass").unwrap();

    // Use proxy variable in config, if set
    match settings.get_str("user.proxy") {
        Ok(proxy) => env::set_var("http_proxy", proxy),
        _ => {}
    }

    let species_key = matches.value_of("SPECIES_KEY").unwrap();
    if species_key == "all" {
        let config: serde_json::Value = reqwest::get(SPECIES_CONFIG_FILE)
            .expect(format!("Config file '{}' can not be loaded", &SPECIES_CONFIG_FILE).as_str())
            .json()
            .expect(
                format!(
                    "Failed to parse config file '{}' into JSON",
                    &SPECIES_CONFIG_FILE,
                )
                .as_str(),
            );
        for (group, species_list) in config.as_object().unwrap() {
            println!("{}", &group);
            for species in species_list.as_array().unwrap() {
                let conf = GeneDBotConfig::new_from_json(species);
                let species_key = conf.abbreviation;
                println!("> {}", &species_key);
                match run_bot_for_species_and_gene(&species_key, &None, &lgname, &lgpass) {
                    Ok(_) => {}
                    Err(e) => println!("RUN FAILED: {:?}", e),
                }
            }
        }
    } else {
        let gene: Option<Vec<String>> = match matches.value_of("gene") {
            Some(gene) => Some(vec![gene.to_string()]),
            None => None,
        };
        run_bot_for_species_and_gene(&species_key.to_string(), &gene, &lgname, &lgpass).unwrap();
    }
}

#[cfg(test)]
mod tests {
    //use super::*;
}
