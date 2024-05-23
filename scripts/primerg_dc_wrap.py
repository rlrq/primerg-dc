## PRIMERg
# Deep-sequencing primer design (with primerg) for multiple targets on a continuous genome

## Source code
## Do not edit unnecessarily

## Gets all gRNA positions, group cleavage sites within max_primary_PCR_amplicon_size of each other.
## Extracts genomic sequence containing cleavages sites per group with max_primary_PCR_amplicon_size flanks on either side of the most upstream and most downstream cleavage site in group.
## Each extracted genomic sequence is passed to primerg_contiguous.py for primer generation.
## Output is merged into a single dataframe and written to file.

output_pd = False

primerg_global_args = {
    'pam', 'gRNA_len', 'cleavage_pos',
    'primer_opt_size', 'primer_min_size', 'primer_max_size',
    'primer_opt_tm', 'primer_min_tm', 'primer_max_tm',
    'primer_min_GC', 'primer_max_GC',
    'primer_dna_conc', 'primer_dNTP_conc', 'primer_salt_divalent',
    'primer_pair_max_diff_TM', 'primer_num_return',
    'primer_num_return_for_unique',
    'evalue', 'word_size', 'max_target_seqs',
    'min_primer_len', 'total_mismatch', 'three_prime_match',
    'primary_PCR_amplicon_size', 'secondary_PCR_amplicon_size',
    'num_primer_pair_per_cleavage_pos', 'increment',
    'max_primary_PCR_amplicon_size',
    'min_acceptable_off_target_primary_amplicon_size',
    'min_acceptable_off_target_secondary_amplicon_size',
    'on_target_ranges', 'verbose',
    'output_excel', 'output_tsv'
}

## Import packages
import sys
import os
import pandas as pd
import tempfile
import csv
import importlib
import logging
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from .index_fasta import IndexedFasta

import argparse

parser = argparse.ArgumentParser(prog = "PRIMERg-dc", description = "Discontiguous PRIMERg")
parser.add_argument('-d', "--dir", "--directory", dest = "directory")
parser.add_argument('-t', "--template", dest = "template_DNA_fasta")
parser.add_argument('-p', "--param", dest = "param_file")
parser.add_argument('-g', "--grna", dest = "gRNA_fasta")
parser.add_argument('-b', "--bg", "--background", dest = "background_DNA_fastas",
                    nargs = '+', action = "extend", default = [])
args = parser.parse_args()

## check if param file was passed to programme
## if no parameter file given, import primerg_parameters from current directory
if args.param_file is None:
    from primerg_parameters import *
## if parameter file is passed to programme
else:
    param_file = os.path.abspath(args.param_file)
    ## try to parse as excel
    try:
        df = pd.read_excel(param_file)
        tmp_tsv = tempfile.NamedTemporaryFile(dir = os.getcwd()).name
        tmp_py = tempfile.NamedTemporaryFile(dir = os.getcwd(), suffix = ".py").name
        # print(tmp_tsv, tmp_py)
        df.to_csv(tmp_tsv, index = False, header = True, sep = '\t', quoting = csv.QUOTE_NONE)
        with open(tmp_tsv, 'r') as f:
            raw_contents = [line[:-1].split('\t') for line in f.readlines()]
        with open(tmp_py, 'w+') as f:
            param_entries = [line for line in raw_contents if line[1] != '' and line[0] != "Parameter"]
            f.write('\n'.join([(line[0] + '=' + ("None" if line[1] == '' else line[1]))
                               for line in param_entries]))
        ## import module by name in variable
        module_path = os.path.dirname(tmp_py)
        sys.path.append(module_path)
        mdl = importlib.import_module(os.path.splitext(os.path.basename(tmp_py))[0])
        if "__all__" in mdl.__dict__:
            names = mdl.__dict__["__all__"]
        else:
            names = [x for x in mdl.__dict__ if not x.startswith("_")]
        globals().update({k: getattr(mdl, k) for k in names})
        ## remove tmp file
        os.remove(tmp_tsv)
        os.remove(tmp_py)
    ## try to directly import file
    except ValueError:
        param_file_dir = os.path.dirname(param_file)
        ## import module by name in variable
        sys.path.append(param_file_dir)
        mdl = importlib.import_module(os.path.splitext(os.path.basename(param_file))[0])
        if "__all__" in mdl.__dict__:
            names = mdl.__dict__["__all__"]
        else:
            names = [x for x in mdl.__dict__ if not x.startswith("_")]
        globals().update({k: getattr(mdl, k) for k in names})

## terminate if no output files
if output_excel is None and output_tsv is None:
    print("At least one of the following output file options is required:")
    print(" - Use 'output_excel' to write output to an Excel sheet")
    print(" - Use 'output_tsv' to write output to a tab-separated text file")
    exit(1)

## parse paths & potentially clashing inputs
if "directory" not in globals() or directory is None:
    directory = os.getcwd() if args.directory is None else os.path.abspath(args.directory)
else:
    directory = os.path.abspath(directory)

if "gRNA_fasta" not in globals() or gRNA_fasta is None:
    gRNA_fasta = None if args.gRNA_fasta is None else os.path.abspath(args.gRNA_fasta)
else:
    gRNA_fasta = os.path.abspath(gRNA_fasta)

if "template_DNA_fasta" not in globals() or template_DNA_fasta is None:
    template_DNA_fasta = None if args.template_DNA_fasta is None else os.path.abspath(args.template_DNA_fasta)
else:
    template_DNA_fasta = os.path.abspath(template_DNA_fasta)

new_background_DNA_fastas = [os.path.abspath(fname) for fname in args.background_DNA_fastas]
if "background_DNA_fastas" not in globals() or background_DNA_fastas is []:
    background_DNA_fastas = new_background_DNA_fastas
else:
    if isinstance(background_DNA_fastas, str):
        background_DNA_fastas = [os.path.abspath(fname) for fname in background_DNA_fastas.split(',')]
    background_DNA_fastas.extend(new_background_DNA_fastas)

## parse other args
if on_target_ranges is None: on_target_ranges = {}

## exit if required inputs are not given
if gRNA_fasta is None:
    print("Fasta file of gRNA sequences is required. Use '--grna <path to file>' at the command line or update gRNA_fasta in the parameters file.")
    exit(1)

if template_DNA_fasta is None:
    print("Fasta file of template is required. Use '--template <path to file>' at the command line or update template_DNA_fasta in the parameters file.")
    exit(1)


## import other primerg modules & update globals
from . import primerg_header as ph
from . import primerg_contiguous as pc

ph.addglobals({varname: val for varname, val in globals().items()
               if varname in primerg_global_args})
pc.addglobals({varname: val for varname, val in globals().items()
               if varname in primerg_global_args})

# test_args = primerg_global_args.union({"gRNA_fasta", "background_DNA_fastas", "template_DNA_fasta", "directory"})
# print({varname: val for varname, val in globals().items()
#        if varname in test_args})

# exit(0)

## use verbose as flag to prevent users from accessing logging levels < INFO
logging_level = logging.INFO if ph.verbose else logging.WARNING

## devs to modify this as necessary
logging.getLogger().setLevel(level = logging_level)

## Reading parameters
os.chdir(directory)
#template read as string
genomic_templates = IndexedFasta(template_DNA_fasta)
ifasta = genomic_templates

# Obtain list of gRNAs and cleavage site positions
gRNA_list = ph.gRNA_listing(gRNA_fasta)
logging.info(''.join(map(str, ["# gRNA: ", len(gRNA_list)])))
gRNA_pam_pos_FR = ph.gRNA_pam_pos_in_IndexedFasta(gRNA_list, pam, ifasta)
num_gRNA_pam_pos = sum(len(subject_dat['F']) + len(subject_dat['R'])
                       for subject, subject_dat in gRNA_pam_pos_FR.items())
if num_gRNA_pam_pos == 0:
    print("No gRNA binding sites detected. Please verify the contents of gRNA_fasta and template_DNA_fasta.")
    exit(1)
logging.info(''.join(map(str, ["# cleavage sites: ", num_gRNA_pam_pos])))
gRNA_seq_cleavage_pos = ph.gRNA_pam_pos_to_seq_and_cleavage_pos(gRNA_pam_pos_FR, ifasta,
                                                                upper = True, as_string = True)

## group cleavage positions
gRNA_cleavage_pos_grouped = {
    subject_id: ph.group_seq_cleavage_pos_within_distance(FR["F"], FR["R"],
                                                          distance = max_primary_PCR_amplicon_size)
    for subject_id, FR in gRNA_seq_cleavage_pos.items()}
logging.info(''.join(map(str, ["gRNA_cleavage_pos_grouped: ", gRNA_cleavage_pos_grouped])))


## print status info
num_cleavage_pos_total = sum((len(cleavage_seq_pos['F'])+len(cleavage_seq_pos['R']))
                             for cleavage_seq_pos_groups in gRNA_cleavage_pos_grouped.values()
                             for cleavage_seq_pos in cleavage_seq_pos_groups)
print(f"--# cleavage positions total: {num_cleavage_pos_total}--")

## get sequences for groups & execute primerg_contiguous
cleavage_pos_processed = 0
df_all = []
for subject_id, cleavage_seq_pos_groups in gRNA_cleavage_pos_grouped.items():
    print(f"-----{subject_id}-----")
    for cleavage_seq_pos in cleavage_seq_pos_groups:
        ## get template DNA sequence for this group of cleavage positions and write to file
        cleavage_positions_FR = [pos for seq, pos in (cleavage_seq_pos['F'] + cleavage_seq_pos['R'])]
        merged_range_start = max(0, min(cleavage_positions_FR) - max_primary_PCR_amplicon_size)
        merged_range_end = max(cleavage_positions_FR) + max_primary_PCR_amplicon_size
        print(f"-----group range: {merged_range_start}-{merged_range_start}-----")
        print(f"-----# cleavage positions in range: {len(cleavage_positions_FR)}-----")
        sequence = ifasta[subject_id][merged_range_start:merged_range_end].upper()
        fasta_template = "primerg_group_template_DNA.fasta"
        SeqIO.write(
            [SeqRecord(Seq(sequence),
                       id = f"group_template|{subject_id}|{merged_range_start+1}-{merged_range_end}",
                       description = f"group_template|{subject_id}|{merged_range_start+1}-{merged_range_end}")],
            fasta_template, "fasta")
        ## this should be in format [(<gRNA + PAM seq>, <cleavage pos>), (<gRNA + PAM seq>, <cleavage pos>)]
        gRNA_pos_adj = [(seq, pos - merged_range_start)
                        for seq, pos in (cleavage_seq_pos['F'] + cleavage_seq_pos['R'])]
        ## execute primerg_contiguous
        df = pc.primerg(gRNA_pos_adj, fasta_template, background_DNA_fastas + [template_DNA_fasta],
                        suppress_output_message = True, print_df = False)
        ## insert column with subject_id as sequence_name
        df.insert(0, "sequence_name", subject_id)
        ## adjust positions from local (in fasta_template) to global (in template_DNA_fasta)
        for i in range(len(df)):
            cleavage_pos = df.loc[i, "cleavage_pos"]
            df.loc[i, "cleavage_pos"] = (cleavage_pos
                                         if not isinstance(cleavage_pos, int)
                                         else cleavage_pos + merged_range_start)
        ## update df_final
        df_all.append(df)
        ## print status update
        cleavage_pos_processed += len(cleavage_positions_FR)
        print(f"--# cleavage positions processed: {cleavage_pos_processed}--")

## merge dataframes
df_final = pd.concat(df_all, ignore_index = True)
## sort dataframe in ascending order
df_final.sort_values(by=["sequence_name", "cleavage_pos"], inplace = True)
logging.debug(str(df_final))

####  WRITE  ####
# Export to excel
if output_excel:
    try:
        df_final.to_excel(output_excel, index = False)
        print(f"Complete. Excel sheet exported as {directory}/{output_excel}")
    except PermissionError:
        print(f"Please close {output_excel} and try again")
## export to tsv
if output_tsv:
    try:
        df_final.to_csv(output_tsv, index = False, sep = '\t')
        print(f"Complete. TSV file exported as {directory}/{output_tsv}")
    except PermissionError:
        print(f"Please close {output_tsv} and try again")

## if output pandas df
if output_pd:
    df_final

exit(0)
