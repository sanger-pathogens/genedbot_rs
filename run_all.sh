#!/bin/bash
\rm /nfs/users/nfs_p/pathpipe/genedbot/genedbot.o
bsub -q long -J genedbot -M 8000 -R "select[mem>8000] rusage[mem=8000]" -o /nfs/users/nfs_p/pathpipe/genedbot/genedbot.o  "cd /nfs/users/nfs_p/pathpipe/genedbot ; /nfs/users/nfs_p/pathpipe/genedbot/target/release/genedbot all"
