## PRIMERg
# Deep-sequencing primer design (with primerg) for multiple targets on a continuous genome

##  Input: (Quick start) [1] gRNA list as fasta, [2] continuous genome seq as fasta and [3] directory to fasta files  
##  Output: Excel sheet of primers if output_excel is True. Pandas dataframe if output_excel is False.

#  Required input
directory = '/mnt/chaelab/rachelle/scripts/primerg/primerg-dc/test2' ## set to None if specifying at command line
gRNA_fasta = "test_gRNA.fasta"                   #list of gRNAs in 5'-3' direction (PAM is not required here) (set to None if specifying at command line)
template_DNA_fasta = "test_genomic.fasta"          # fasta file of genome in which gRNA is meant to be identified (set to None if specifying at command line)
background_DNA_fastas = []                         # list of fasta file of background genomes against which specificity is to be checked (optional if template_DNA_fasta is the only background)
# gRNA_fasta = "test_gRNA_ADR1-NRG1.fasta"                   #list of gRNAs in 5'-3' direction (PAM is not required here)
# genomic_DNA_fasta = "/mnt/chaelab/shared/genomes/TAIR10/fasta/a_thaliana_all.fasta"          #full fasta file of genome against which to check specificity



#  Optional input

#  CRISPR-Cas parameters
pam = "NGN"                           #e.g. "NGG", let pam be '' if pamless, follows IUPAC DNA (e.g. R = A/G)
gRNA_len = 20
cleavage_pos = 3                      #nth base from 3' end of gRNA where Cas9 cleaves genomic template,
                                      #default = 3 where NNNNNNNNNNNNNNNNN|NNN
#  Primer3
primer_opt_size = 20                  #recommend 20 bp
primer_min_size = 18                  #recommend 18 bp
primer_max_size = 30                  #recommend 30 bp
primer_opt_tm = 65                    #recommend 65
primer_min_tm = 60                    #recommend 60
primer_max_tm = 70                    #recommend 70
primer_min_GC = 30                    #recommend 40%; 30% for leniency ###check###
primer_max_GC = 60                    #recommend 60%
primer_dna_conc = 50                  #in nm
primer_dNTP_conc = 0.5                #in mM
primer_pair_max_diff_TM = 5           #recommend 5
primer_salt_divalent = 1.5            #mM of divalent cations e.g. Mg2+
primer_num_return = 500                #maximum number of primer pairs to return per primer3 run
primer_num_return_for_unique = 20      #maximum number of primer partners to return per primer3 run for unique primers

#  Pseudo-Primer-BLAST
evalue = 1000
word_size = 7
max_target_seqs = 50000
min_primer_len = 10                    # The minimum length of alignment between primer and template to consider a primer as unspecific; default 10 bp is minimum for annealing prior to extension
total_mismatch = 6                     # The minimum number of mismatches between primer and template to consider a primer as unspecific; default is 6 (each unaligned position counts as a mismatch)
three_prime_match = 3                  # At least X matches within the last 5 bp at the 3' end; default is 3 (counts from the end of the primer; each unaligned position at the 3' end counts as a mismatch)
valid_amplicon_size = 500              # Size of amplicon produced by unspecific primers to be considered valid amplicon, default is 500 ***************** 500  (defunct. I'm not sure it was even used by the original primerg)
NGS_amplicon_size = 750                # Length of amplicon of EACH primer in paired-end sequencing; default is set at 150 for iSeq **********  (defunct. see secondary_PCR_amplicon_size)

#  Primerg
output_excel = "sample_output.xlsx"           #if an excel name is specified (e.g. filename.xlsx), output will be an excel sheet.
output_tsv = "sample_output.tsv"              #if an tsv name is specified (e.g. filename.tsv), output will be an tsv file. (this is IN ADDITION to output_excel)
## If both output_excel and output_tsv are False, output will be a pandas dataframe.

#  New inputs
primary_PCR_amplicon_size = [1200, 2000] # First PCR product size (size of amplicon as template for second PCR)
secondary_PCR_amplicon_size = [250, 280] # Second PCR product size (size of amplicon for sequencing)
num_primer_pair_per_cleavage_pos = 10    # Generate x number of primer pairs per cleavage site
increment = 300                          # Increment of allowed primary amplicon size when no primer can be generated
max_primary_amplicon_size = 3500         # Maximum primary amplicon size
min_acceptable_off_target_primary_amplicon_size = float("Inf")    # Length of minimum acceptable off-target amplicon size. Definition from https://manual.geneious.com/en/latest/13-Primers.html: 'Off-target primer pairs that result in amplicon sizes larger than specified value are considered as primer candidates due to decreasing efficiency of PCR as the amplicon size increases.' Set to float("Inf") to discard off-target primer pairs regardless of off-target amplicon size.
min_acceptable_off_target_secondary_amplicon_size = float("Inf")    # Length of minimum acceptable off-target amplicon size. Definition from https://manual.geneious.com/en/latest/13-Primers.html: 'Off-target primer pairs that result in amplicon sizes larger than specified value are considered as primer candidates due to decreasing efficiency of PCR as the amplicon size increases.' Set to float("Inf") to discard off-target primer pairs regardless of off-target amplicon size.
on_target_ranges = {}               # dict of {<subject_id>: [(<start>, <end>), (<start>, <end>)]}; ignore any off-target amplicons within these ranges (both primers of a pair must be within)

#######################################################################################################################


## Miscellaneous printing options
verbose = True          ## controls verbosity of printout to screen; set to False for less verbose execution
