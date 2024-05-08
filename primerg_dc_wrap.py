## PRIMERg
# Deep-sequencing primer design (with primerg) for multiple targets on a continuous genome

## Source code
## Do not edit unnecessarily

## Gets all gRNA positions, group cleavage sites within max_primary_amplicon_size of each other.
## Extracts genomic sequence containing cleavages sites per group with max_primary_amplicon_size flanks on either side of the most upstream and most downstream cleavage site in group.
## Each extracted genomic sequence is passed to primerg_contiguous.py for primer generation.
## Output is merged into a single dataframe and written to file.

## Import packages
from index_fasta import IndexedFasta
from primerg_contiguous import *

## terminate if no output files
if output_excel is None and output_tsv is None:
    print("At least one of the following output file options is required:")
    print(" - Use 'output_excel' to write output to an Excel sheet")
    print(" - Use 'output_tsv' to write output to a tab-separated text file")
    exit(1)

## use verbose as flag to prevent users from accessing logging levels < INFO
logging_level = logging.INFO if verbose else logging.WARNING

## devs to modify this as necessary
logging.getLogger().setLevel(level = logging_level)

## Reading parameters
os.chdir(directory)
#template read as string
genomic_templates = IndexedFasta(genomic_DNA_fasta)
ifasta = genomic_templates

# Obtain list of gRNAs and cleavage site positions
gRNA_list = gRNA_listing(gRNA_fasta)
logging.info(''.join(map(str, ["# gRNA: ", len(gRNA_list)])))
gRNA_pam_pos_FR = gRNA_pam_pos_in_IndexedFasta(gRNA_list, pam, ifasta)
num_gRNA_pam_pos = sum(len(subject_dat['F']) + len(subject_dat['R'])
                       for subject, subject_dat in gRNA_pam_pos_FR.items())
if num_gRNA_pam_pos == 0:
    print("No gRNA binding sites detected. Please verify the contents of gRNA_fasta and genomic_DNA_fasta.")
    exit(1)
logging.info(''.join(map(str, ["# cleavage sites: ", num_gRNA_pam_pos])))
gRNA_seq_cleavage_pos = gRNA_pam_pos_to_seq_and_cleavage_pos(gRNA_pam_pos_FR, ifasta, upper = True, as_string = True)

## group cleavage positions
gRNA_cleavage_pos_grouped = {
    subject_id: group_seq_cleavage_pos_within_distance(FR["F"], FR["R"],
                                                       distance = max_primary_amplicon_size)
    for subject_id, FR in gRNA_seq_cleavage_pos.items()}
logging.info(''.join(map(str, ["gRNA_cleavage_pos_grouped: ", gRNA_cleavage_pos_grouped])))


## get sequences for groups & execute primerg_contiguous
df_all = []
for subject_id, cleavage_seq_pos_groups in gRNA_cleavage_pos_grouped.items():
    print(f"-----{subject_id}-----")
    for cleavage_seq_pos in cleavage_seq_pos_groups:
        ## get template DNA sequence for this group of cleavage positions and write to file
        cleavage_positions_FR = [pos for seq, pos in (cleavage_seq_pos['F'] + cleavage_seq_pos['R'])]
        merged_range_start = max(0, min(cleavage_positions_FR) - max_primary_amplicon_size)
        merged_range_end = max(cleavage_positions_FR) + max_primary_amplicon_size
        sequence = ifasta[subject_id][merged_range_start:merged_range_end].upper()
        fasta_template = "primerg_group_template_DNA.fasta"
        SeqIO.write([SeqRecord(Seq(sequence),
                               id = f"group_template|{subject_id}|{merged_range_start+1}-{merged_range_end}",
                               description = f"group_template|{subject_id}|{merged_range_start+1}-{merged_range_end}")],
                    fasta_template, "fasta")
        ## this should be in format [(<gRNA + PAM seq>, <cleavage pos>), (<gRNA + PAM seq>, <cleavage pos>)]
        gRNA_pos_adj = [(seq, pos - merged_range_start) for seq, pos in (cleavage_seq_pos['F'] + cleavage_seq_pos['R'])]
        ## execute primerg_contiguous
        df = primerg(gRNA_pos_adj, fasta_template, suppress_output_message = True, print_df = False)
        ## insert column with subject_id as sequence_name
        df.insert(0, "sequence_name", subject_id)
        ## adjust positions from local (in fasta_template) to global (in genomic_DNA_fasta)
        for i in range(len(df)):
            cleavage_pos = df.loc[i, "cleavage_pos"]
            df.loc[i, "cleavage_pos"] = cleavage_pos if not isinstance(cleavage_pos, int) else cleavage_pos + merged_range_start
        ## update df_final
        df_all.append(df)

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

exit(0)
