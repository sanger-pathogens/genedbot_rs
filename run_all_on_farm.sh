#!/bin/bash
mailto="pathdevg@sanger.ac.uk"
jobname="genedbot"
genedbot_base="/nfs/users/nfs_p/pathdb/genedbot"
errfile="${genedbot_base}.e"
outfile="${genedbot_base}.o"
bsub -q long -eo ${errfile} -oo ${outfile} -J ${jobname} -M3500 -R "select[mem>3500] rusage[mem=3500]" "http_proxy=http://wwwcache.sanger.ac.uk:3128 RUST_BACKTRACE=1 genedbot --config /nfs/users/nfs_p/pathdb/genedbot.ini all"
bsub -o /dev/null -e /dev/null -w "ended("${jobname}")" "if [ -s $errfile ]; then cat $errfile | mailx -s \"GeneDBot error\" ${mailto}; fi"
