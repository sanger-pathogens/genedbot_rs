[![Build Status](https://travis-ci.org/sanger-pathogens/genedbot.svg)](https://travis-ci.org/sanger-pathogens/genedbot) [![codecov](https://codecov.io/gh/sanger-pathogens/genedbot/branch/master/graph/badge.svg)](https://codecov.io/gh/sanger-pathogens/genedbot)

This is Rust code to update Wikidata items from GeneDB GFF and GAF files.

# Installation
1. [Install Rust](https://rustup.rs/) (it's a one-liner)
2. `git clone https://github.com/sanger-pathogens/genedbot.git`
3. `cd genedbot`
4. `cargo build --release`
5. Create a `bot.ini` file, like so:
```
[user]
user = GeneDBot
pass = INSERT_PASSWORD_HERE
proxy = http://wwwcache.sanger.ac.uk:3128
```

# Usage:
1. As pathpipe@pathpipe-farm4, in `~/genedbot` directory
2. `./target/release/genedbot SPECIES_CODE` (`SPECIES_CODE` eg Pfalciparum), _or_
3. `./target/release/genedbot all` to run all species sequentially (use `run_all.sh` to start this via `bsub` on farm4), _or_
4. `./target/release/genedbot --help` for options

# Update
1. As pathpipe@pathpipe-farm4, in `~/genedbot` directory
2. `git pull ; rustup update ; cargo update ; cargo build --release`
