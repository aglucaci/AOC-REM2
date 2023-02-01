"""
2022 - Analysis of Orthologous Collections (AOC).
@Author: Alexander G. Lucaci
"""

#----------------------------------------------------------------------------
# Imports 
#----------------------------------------------------------------------------

import itertools
import os
import sys
import csv
import json
from pathlib import Path
from snakemake.utils import min_version
min_version("6.3")
from Bio import Entrez
from ete3 import NCBITaxa
import pandas as pd
from ete3 import Tree
import glob

#----------------------------------------------------------------------------
# Configuration
#----------------------------------------------------------------------------

configfile: 'config.yml'

print("# Loaded config yaml file")

with open("cluster.json", "r") as fh_cluster:
  cluster = json.load(fh_cluster)
  fh_cluster.close()
#end with

print("# Loaded cluster json file")

Nucleotide_file = config["Nucleotide"]
Protein_file = config["Protein"]
Label = config["Label"]

print("# Using nucleotide data from:", Nucleotide_file)
print("# Using protein data from:", Protein_file)
print("# Using the analysis label:", Label)

HumanRef = config["HumanRef"]

# Set output directory
BASEDIR = os.getcwd()

print("# We are operating out of base directory:", BASEDIR)
OUTDIR = os.path.join(BASEDIR, "results", Label)
OUTDIR_RESULTS = os.path.join(BASEDIR, "results")

print("# We will create and store results in:", OUTDIR)

# Create output dir.
os.makedirs(OUTDIR_RESULTS, exist_ok = True)
print("# Directory '% s' created" % OUTDIR_RESULTS)

os.makedirs(OUTDIR, exist_ok = True)
print("# Directory '% s' created" % OUTDIR)

#Path(OUTDIR_RESULTS).mkdir(parents=True, exist_ok=True)
#Path(OUTDIR).mkdir(parents=True, exist_ok=True)

PPN = cluster["__default__"]["ppn"]

#PPN_MSA = int(PPN) - 7

# Batch files
PREMSA = os.path.join(BASEDIR, config["PREMSA"])
POSTMSA = os.path.join(BASEDIR, config["POSTMSA"])

# Hard-coded HyPhy settings
HYPHY = config["HYPHY"]
HYPHYMPI = config["HYPHYMPI"]
#RES = os.path.join(BASEDIR, config["RES"])

CODON_OUTPUT = os.path.join(OUTDIR, Label)

# Clustering
TN93_T = config["TN93_Threshold"]

CSV = os.path.join(BASEDIR, config["CSV"])

#----------------------------------------------------------------------------
# Rule all
#----------------------------------------------------------------------------

rule all:
    input:
        CODON_OUTPUT,
        os.path.join(OUTDIR, Label + "_protein.fas"),
        os.path.join(OUTDIR, Label + "_nuc.fas"),
        os.path.join(OUTDIR, Label + "_protein.aln"),
        os.path.join(OUTDIR, Label + "_codons.fasta"),
        os.path.join(OUTDIR, Label + "_codons_duplicates.json"),
        os.path.join(OUTDIR, Label + "_codons.SA.fasta.dst"),
        os.path.join(OUTDIR, Label + "_codons.SA.fasta"),
        os.path.join(OUTDIR, Label + "_codons.SA.fasta.treefile"),
        os.path.join(OUTDIR, Label + "_codons.SA.fasta.BUSTEDS.json"),
        os.path.join(OUTDIR, Label + "_codons.SA.fasta.BUSTED.json"),
        os.path.join(OUTDIR, Label + "_codons.SA.fasta.BUSTEDS+MH.json"),
        os.path.join(OUTDIR, Label + "_codons.SA.fasta.BUSTED+MH.json"),
        #os.path.join(OUTDIR, Label + "_codons.SA.fasta.FITMG94.json"),
        os.path.join(OUTDIR, Label + "_codons.SA.fasta.BGM.json"),
        os.path.join(OUTDIR, Label + "_codons.SA.fasta.FEL.json"),
        os.path.join(OUTDIR, Label + "_codons.SA.fasta.FUBAR.json"),
        os.path.join(OUTDIR, Label + "_codons.SA.fasta.FMM.json"),
        os.path.join(OUTDIR, Label + "_codons.SA.fasta.MEME.json"),
        os.path.join(OUTDIR, Label + "_codons.SA.fasta.SLAC.json"),
        os.path.join(OUTDIR, Label + "_codons.SA.fasta.aBSRELS.json"),
        os.path.join(OUTDIR, Label + "_codons.SA.fasta.aBSRELS+MH.json"),
        #os.path.join(OUTDIR, Label + "_codons.SA.fasta.treefile.labelled"),
        os.path.join(OUTDIR, Label + "_codons.SA.fasta.GARD.json"),
        #os.path.join(OUTDIR, Label + "_codons.SA.fasta.RELAX.json"),
        #os.path.join(OUTDIR, Label + "_codons.SA.fasta.CFEL.json") 
#end rule all

print("# Moving on to processing rules")

#----------------------------------------------------------------------------
# Rules
#----------------------------------------------------------------------------

rule get_codons:
    output:
        codons = CODON_OUTPUT, 
    params:
        Nuc = Nucleotide_file,
        Prot = Protein_file,
        Out = CODON_OUTPUT,
        Logfile = CODON_OUTPUT + ".log"
    script:
        "scripts/codons.py"
#end rule

#----------------------------------------------------------------------------
# Alignment
#----------------------------------------------------------------------------

rule pre_msa:
    input: 
        codons = rules.get_codons.output.codons
    output: 
        protein_fas = os.path.join(OUTDIR, Label + "_protein.fas"),
        nucleotide_fas = os.path.join(OUTDIR, Label + "_nuc.fas")
    shell: 
        "mpirun -np {PPN} {HYPHYMPI} {PREMSA} --input {input.codons} --reference {HumanRef}"
        #"mpirun -np {PPN} {HYPHYMPI} {PREMSA} --input {input.codons}"
#end rule 

rule mafft:
    input:
        protein = rules.pre_msa.output.protein_fas
    output:
        protein_aln = os.path.join(OUTDIR, Label + "_protein.aln")
    shell:
        "mafft --auto {input.protein} > {output.protein_aln}"
#end rule

rule post_msa:
    input: 
        protein_aln = rules.mafft.output.protein_aln,
        nucleotide_seqs = rules.pre_msa.output.nucleotide_fas  
    output: 
        codons_fas = os.path.join(OUTDIR, Label + "_codons.fasta"),
        duplicates_json = os.path.join(OUTDIR, Label + "_codons_duplicates.json")
    shell: 
        "mpirun -np {PPN} {HYPHYMPI} {POSTMSA} --protein-msa {input.protein_aln} --nucleotide-sequences {input.nucleotide_seqs} --output {output.codons_fas} --duplicates {output.duplicates_json}"
#end rule 

#----------------------------------------------------------------------------
# Remove ambiguous codons, a source of noise.
#----------------------------------------------------------------------------

rule strike_ambigs:
   input:
       in_msa = rules.post_msa.output.codons_fas
   output:
       out_strike_ambigs = os.path.join(OUTDIR, Label + "_codons.SA.fasta")  
   shell:
       #"hyphy scripts/strike-ambigs.bf --alignment {input.in_msa} --output {output.out_strike_ambigs}"
       "{HYPHY} scripts/strike-ambigs.bf --alignment {input.in_msa} --output {output.out_strike_ambigs}"
#end rule 

#----------------------------------------------------------------------------
# IQ-TREE for ML tree inference
#----------------------------------------------------------------------------

rule iqtree: # Unrooted
    input:
        codons_fas = rules.strike_ambigs.output.out_strike_ambigs
    output:
        tree = os.path.join(OUTDIR, Label + "_codons.SA.fasta.treefile")  
    shell:
        "iqtree -s {input.codons_fas} -T AUTO"
    #end shell
#end rule iqtree

#----------------------------------------------------------------------------
# TN93, genetic distance calculation
#----------------------------------------------------------------------------

rule tn93:
    input:
       input = rules.strike_ambigs.output.out_strike_ambigs
    output:
       output = os.path.join(OUTDIR, Label + "_codons.SA.fasta.dst")
    shell:
       "tn93 -t 1 -o {output.output} {input.input}"
    #end shell
#end rule

#----------------------------------------------------------------------------
# Downsample for GARD
#----------------------------------------------------------------------------

rule tn93_cluster:
    input:
        input = rules.strike_ambigs.output.out_strike_ambigs
    output:
        output = os.path.join(OUTDIR, Label + "_codons.SA.fasta.cluster.json")
    shell:
        "tn93-cluster -f -o {output.output} -t {TN93_T} {input.input}" 
#end rule

rule cluster_to_fasta:
   input: 
       input = rules.tn93_cluster.output.output
   output:
       output = os.path.join(OUTDIR, Label + "_codons.SA.fasta.cluster.fasta")
   shell:
       "python scripts/cluster_to_fasta.py -i {input.input} -o {output.output}"
#end rule 

#----------------------------------------------------------------------------
# Recombination detection
#----------------------------------------------------------------------------

rule recombination_original:
    input: 
        input = rules.strike_ambigs.output.out_strike_ambigs 
    output: 
        output =  os.path.join(OUTDIR, Label + "_codons.SA.fasta.GARD.json") 
    shell: 
        "mpirun -np {PPN} {HYPHYMPI} ENV='TOLERATE_NUMERICAL_ERRORS=1;' GARD --alignment {input.input} --rv GDD --output {output.output}"
#end rule


rule recombination:
    input: 
        input = rules.cluster_to_fasta.output.output 
    output: 
        output = os.path.join(OUTDIR, Label + "_codons.SA.fasta.cluster.fasta.GARD.json")
    shell: 
        "mpirun --use-hwthread-cpus -np {PPN} {HYPHYMPI} GARD --alignment {input.input} --rv GDD --output {output.output}"
#end rule


#rule recombination_clean:
#    input: 
#        input =  rules.strike_ambigs.output.out_strike_ambigs 
#    output: 
#        output = os.path.join(OUTDIR, Label + "_codons.SA.fasta.GARD.json")
#    shell: 
#        "mpirun --use-hwthread-cpus -np {PPN} {HYPHYMPI} GARD --alignment {input.input} --rv GDD --output {output.output}"
#end rule

#----------------------------------------------------------------------------
# Split out GARD partitions
#----------------------------------------------------------------------------

rule gard_parse:
    input:
        input = rules.recombination.output.output
    params:
        genelabel = Label,
        base_dir = BASEDIR
    output:
        output = os.path.join(OUTDIR, Label + ".1.codon.fas")
    script:
        "scripts/GARD_Parser.py"
#end rule

#----------------------------------------------------------------------------
# Selection analyses
#----------------------------------------------------------------------------

rule BUSTEDS:
    input: 
        codon_aln = rules.strike_ambigs.output.out_strike_ambigs,
        tree = rules.iqtree.output.tree
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.SA.fasta.BUSTEDS.json")
    shell: 
        "mpirun -np {PPN} {HYPHYMPI} BUSTED --alignment {input.codon_aln} --tree {input.tree} --output {output.results} --srv Yes --starting-points 10"
#end rule

rule BUSTED:
    input: 
        codon_aln = rules.strike_ambigs.output.out_strike_ambigs,
        tree = rules.iqtree.output.tree
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.SA.fasta.BUSTED.json")
    shell: 
        "mpirun -np {PPN} {HYPHYMPI} BUSTED --alignment {input.codon_aln} --tree {input.tree} --output {output.results} --srv Yes --starting-points 10"
#end rule

rule BUSTEDSMH:
    input: 
        codon_aln = rules.strike_ambigs.output.out_strike_ambigs,
        tree = rules.iqtree.output.tree
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.SA.fasta.BUSTEDS+MH.json")
    shell: 
        "mpirun -np {PPN} {HYPHYMPI} BUSTED --alignment {input.codon_aln} --tree {input.tree} --output {output.results} --srv Yes --starting-points 10 --multiple-hits Double+Triple"
#end rule

rule BUSTEDMH:
    input: 
        codon_aln = rules.strike_ambigs.output.out_strike_ambigs,
        tree = rules.iqtree.output.tree
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.SA.fasta.BUSTED+MH.json")
    shell: 
        "mpirun -np {PPN} {HYPHYMPI} BUSTED --alignment {input.codon_aln} --tree {input.tree} --output {output.results} --srv Yes --starting-points 10 --multiple-hits Double+Triple"
#end rule

rule BGM:
    input: 
        codon_aln = rules.strike_ambigs.output.out_strike_ambigs,
        tree = rules.iqtree.output.tree
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.SA.fasta.BGM.json")
    shell: 
        #"mpirun -np {PPN} {HYPHYMPI} BGM --alignment {input.codon_aln} --tree {input.tree} --output {output.results}"
        "{HYPHY} BGM --alignment {input.codon_aln} --tree {input.tree} --output {output.results}"
#end rule 

rule SLAC:
    input: 
        codon_aln = rules.strike_ambigs.output.out_strike_ambigs,
        tree = rules.iqtree.output.tree
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.SA.fasta.SLAC.json")
    shell: 
        "mpirun -np {PPN} {HYPHYMPI} SLAC --alignment {input.codon_aln} --tree {input.tree} --output {output.results}"
#end rule 

rule ABSRELS:
    input: 
        codon_aln = rules.strike_ambigs.output.out_strike_ambigs,
        tree = rules.iqtree.output.tree
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.SA.fasta.aBSRELS.json")
    shell: 
        "mpirun -np {PPN} {HYPHYMPI} ABSREL --alignment {input.codon_aln} --tree {input.tree} --output {output.results} --srv Yes"
#end rule

rule ABSRELSMH:
    input: 
        codon_aln = rules.strike_ambigs.output.out_strike_ambigs,
        tree = rules.iqtree.output.tree
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.SA.fasta.aBSRELS+MH.json")
    shell: 
        "mpirun -np {PPN} {HYPHYMPI} ABSREL --alignment {input.codon_aln} --tree {input.tree} --output {output.results} --srv Yes --multiple-hits Double+Triple"
#end rule

rule FITMG94:
    input: 
        codon_aln = rules.strike_ambigs.output.out_strike_ambigs,
        tree = rules.iqtree.output.tree
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.SA.fasta.FITMG94.json")
    shell: 
        "{HYPHY} {FITMG94} --alignment {input.codon_aln} --tree {input.tree} --output {output.results} --rooted No --lrt Yes --type global --frequencies CF3x4"
#end rule 


# --- Site methods ---

rule FMM:
    input: 
        codon_aln = rules.strike_ambigs.output.out_strike_ambigs,
        tree = rules.iqtree.output.tree
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.SA.fasta.FMM.json")
    shell: 
        "mpirun -np {PPN} {HYPHYMPI} FMM --alignment {input.codon_aln} --tree {input.tree} --output {output.results} --triple-islands Yes"
#end rule


rule MEME:
    input: 
        codon_aln = rules.strike_ambigs.output.out_strike_ambigs,
        tree = rules.iqtree.output.tree
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.SA.fasta.MEME.json")
    shell: 
        "mpirun -np {PPN} {HYPHYMPI} MEME --alignment {input.codon_aln} --tree {input.tree} --output {output.results}"
#end rule 


rule FEL:
    input: 
        codon_aln = rules.strike_ambigs.output.out_strike_ambigs,
        tree = rules.iqtree.output.tree 
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.SA.fasta.FEL.json")
    shell: 
        "mpirun -np {PPN} {HYPHYMPI} FEL --alignment {input.codon_aln} --tree {input.tree} --output {output.results} --ci Yes"
#end rule 

rule FUBAR:
    input: 
        codon_aln = rules.strike_ambigs.output.out_strike_ambigs,
        tree = rules.iqtree.output.tree
    output: 
        results = os.path.join(OUTDIR, Label + "_codons.SA.fasta.FUBAR.json")
    shell: 
        "mpirun -np {PPN} {HYPHYMPI} FUBAR --alignment {input.codon_aln} --tree {input.tree} --output {output.results}"
#end rule 

#----------------------------------------------------------------------------
# Lineages
#----------------------------------------------------------------------------
rule GatherLineages:
    input:
        out_d  = OUTDIR,
        csv_f  = CSV,
        tree_f = rules.iqtree.output.tree
    output:
        output = os.path.join(OUTDIR, Label + "_codons.SA.fasta.treefile.log")  
    shell:
        "python scripts/LineageAnnotation_Pipeline.py {input.out_d} {input.csv_f} {input.tree_f}"
#end rule

CLADE_FILES = [x for x in glob.glob(os.path.join(OUTDIR, "*.clade"))]
print("# We have", len(CLADE_FILES), "clade files")
print(CLADE_FILES) 

rule AssignLineages:
    input:
        tree = rules.iqtree.output.tree,
        log  = rules.GatherLineages.output.output
    output:
        output = os.path.join(OUTDIR, Label + "_codons.SA.fasta.treefile.labelled")
    run:
        first_time = True
        for clade_file in CLADE_FILES:
            print(clade_file, input[0])
            label      = os.path.basename(clade_file).split(".")[0]

            if first_time == True:
            	cmd = " ".join([HYPHY, 
                                os.path.join(BASEDIR, "scripts", "label-tree.bf"),
                            	"--tree", input[0],
                            	"--list", clade_file,
                            	"--output", output[0],
                            	"--label", label])
                first_time = False
            else:
                cmd = " ".join([HYPHY, 
                                os.path.join(BASEDIR, "scripts", "label-tree.bf"),
                            	"--tree", output[0],
                            	"--list", clade_file,
                            	"--output", output[0],
                            	"--label", label])
            #end if
            print(cmd)
            os.system(cmd)
        #end for
    #end run
#end rule

rule RELAX:
    input:
        treefile = rules.AssignLineages.output.output,
        fasta    = rules.strike_ambigs.output.out_strike_ambigs
    output:
        output = os.path.join(OUTDIR, Label + "_codons.SA.fasta.RELAX.json")
    shell:
        "{HYPHY} RELAX --alignment {input.fasta} --tree {input.treefile} --output {output.output} --reference-group Primates --models All --mode 'Group mode' --starting-points 10 --srv Yes"
#end rule

rule CFEL:
    input:  
        treefile = rules.AssignLineages.output.output,
        fasta    = rules.strike_ambigs.output.out_strike_ambigs
    output: 
        output = os.path.join(OUTDIR, Label + "_codons.SA.fasta.CFEL.json")
    shell:
        "{HYPHY} contrast-fel --alignment {input.fasta} --tree {input.treefile} --output {output.output} --branch-set Primates"
#end file



#----------------------------------------------------------------------------
# End of file
#----------------------------------------------------------------------------

