#[macro_use]
extern crate lazy_static;
#[macro_use]
extern crate serde_json;
extern crate chrono;
extern crate config;
extern crate libflate;
extern crate mediawiki;
extern crate papers;
extern crate regex;
extern crate reqwest;

use crate::genedbot::*;

use config::{Config, File};

use std::env;
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
) -> Result<(), Box<Error>> {
    let mut bot = GeneDBot::new();
    //bot.simulate = true;
    //bot.verbose = true;
    bot.api_mut().set_user_agent("GeneDBot/3.0");
    bot.api_mut().set_edit_delay(Some(500)); // Half a second between edits
    bot.specific_genes_only = genes.to_owned();
    bot.api_mut().login(lgname, lgpass)?;
    bot.load_config_file(species_key)?;
    bot.init()?;
    bot.run()?;
    Ok(())
}

fn main() {
    let args: Vec<_> = env::args().collect();
    if args.len() < 2 {
        println!("Argument (species key) required\n");
        return;
    }

    let mut settings = Config::default();
    settings.merge(File::with_name("bot.ini")).unwrap();
    let lgname = settings.get_str("user.user").unwrap();
    let lgpass = settings.get_str("user.pass").unwrap();

    let species_key = &args[1];
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
        let gene = if args.len() == 3 {
            Some(vec![args[2].to_string()])
        } else {
            None
        };
        run_bot_for_species_and_gene(&species_key, &gene, &lgname, &lgpass).unwrap();
    }
}

#[cfg(test)]
mod tests {
    //use super::*;
}
