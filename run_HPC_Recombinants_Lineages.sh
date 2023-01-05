#!/bin/bash

set -euo pipefail

printf "Running snakemake...\n"

#snakemake --forceall --dag | dot -Tpdf > dag.pdf

snakemake \
      -s Snakefile_Recombinants_Lineages \
      --cluster-config cluster.json \
      --cluster "qsub -V -l nodes={cluster.nodes}:ppn={cluster.ppn} -q {cluster.name} -l walltime={cluster.walltime} -e logs -o logs" \
      --jobs 20 all \
      --rerun-incomplete \
      --keep-going \
      --reason \
      --latency-wait 60 
