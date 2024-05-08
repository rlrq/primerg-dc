#### FUNCTIONS AND WHATNOT ####

## Source code
#  Do not edit unnecessarily

#  Required input
directory = None
gRNA_fasta = None                 #list of gRNAs in 5'-3' direction (PAM is not required here)
genomic_DNA_fasta = None          #full fasta file of genome against which to check specificity
# template_DNA_fasta = None         #fasta file of a single sequence trimmed around one or more cleavage sites


#  Optional input

#  CRISPR-Cas parameters
pam = None                          #e.g. "NGG", let pam be '' if pamless, follows IUPAC DNA (e.g. R = A/G)
gRNA_len = None
cleavage_pos = None                    #nth base from 3' end of gRNA where Cas9 cleaves genomic template,
                                      #default = 3 where NNNNNNNNNNNNNNNNN|NNN
#  Primer3
primer_opt_size = None                 #recommend 20 bp
primer_min_size = None                 #recommend 18 bp
primer_max_size = None                  #recommend 30 bp
primer_opt_tm = None                    #recommend 65
primer_min_tm = None                    #recommend 60
primer_max_tm = None                   #recommend 70
primer_min_GC = None                    #recommend 40%; 30% for leniency ###check###
primer_max_GC = None                   #recommend 60%
primer_dna_conc = None                 #in nm
primer_dNTP_conc = None                #in mM
primer_pair_max_diff_TM = None           #recommend 5
primer_salt_divalent = None            #mM of divalent cations e.g. Mg2+
primer_num_return = None                #maximum number of primers to return per primer3 run 

#  Pseudo-Primer-BLAST
evalue = 1000
word_size = 7
max_target_seqs = 50000
min_primer_len = None                   # The minimum length of alignment between primer and template to consider a primer as unspecific; default 10 bp is minimum for annealing prior to extension
total_mismatch = None                    # The minimum number of mismatches between primer and template to consider a primer as unspecific; default is 6 (each unaligned position counts as a mismatch)
three_prime_match = None                  # At least X matches within the last 5 bp at the 3' end; default is 3 (counts from the end of the primer; each unaligned position at the 3' end counts as a mismatch)
valid_amplicon_size = None             # Size of amplicon produced by unspecific primers to be considered valid amplicon, default is 500 ***************** 500  (defunct. I'm not sure it was even used by the original primerg)
NGS_amplicon_size = None               # Length of amplicon of EACH primer in paired-end sequencing; default is set at 150 for iSeq **********  (defunct. see secondary_PCR_amplicon_size)

#  Primerg
output_excel = None           #if an excel name is specified (e.g. filename.xlsx), output will be an excel sheet.
output_tsv = None              #if an tsv name is specified (e.g. filename.tsv), output will be an tsv file. (this is IN ADDITION to output_excel)
## If both output_excel and output_tsv are False, output will be a pandas dataframe.

#  New inputs
## [rlrq: idk why primary_PCR_amplicon_size and secondary_PCR_amplicon_size are lists and not just integer values, but I won't change them in case it breaks the programme somehow]
primary_PCR_amplicon_size = None # First PCR product size (size of amplicon as template for second PCR)
secondary_PCR_amplicon_size = None# Second PCR product size (size of amplicon for sequencing)
num_primer_pair_per_cleavage_pos = None     # Generate x number of primer pairs per cleavage site
increment = None                         # Increment of allowed primary amplicon size when no primer can be generated
max_primary_amplicon_size = None      # Maximum primary amplicon size
min_acceptable_off_target_primary_amplicon_size = None    # Length of minimum acceptable off-target amplicon size. Definition from https://manual.geneious.com/en/latest/13-Primers.html: 'Off-target primer pairs that result in amplicon sizes larger than specified value are considered as primer candidates due to decreasing efficiency of PCR as the amplicon size increases.' Set to float("Inf") to discard off-target primer pairs regardless of off-target amplicon size.
min_acceptable_off_target_secondary_amplicon_size = None    # Length of minimum acceptable off-target amplicon size. Definition from https://manual.geneious.com/en/latest/13-Primers.html: 'Off-target primer pairs that result in amplicon sizes larger than specified value are considered as primer candidates due to decreasing efficiency of PCR as the amplicon size increases.' Set to float("Inf") to discard off-target primer pairs regardless of off-target amplicon size.
on_target_ranges = None               # dict of {<subject_id>: [(<start>, <end>), (<start>, <end>)]}; ignore any off-target amplicons within these ranges (both primers of a pair must be within)

## Import packages
from primerg_parameters import *
from primerg_classes import *
from sys import exit
import os
import primer3
import tempfile
import pandas as pd
import regex as re
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
# from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast import NCBIXML
from Bio import SearchIO
from Bio.Data import IUPACData

import copy
import subprocess

## Miscellaneous constants
IUPAC_dict = {amb_letter: (disamb if len(disamb) == 1 else f"[{disamb}]")
              for amb_letter, disamb in IUPACData.ambiguous_dna_values.items()}

# ## Reading parameters
# os.chdir(directory)
# genomic_template = str(SeqIO.read(template_DNA_fasta, "fasta").seq).upper()    #template read as string

## Read gRNA fasta file into a list
def gRNA_listing(gRNA_fasta):
    # Time: O(n) where n is length of gRNA list
    # Space: O(n) where n is length of gRNA list
    gRNA_list = []
    try:
        gRNA_fasta_list = SeqIO.parse(gRNA_fasta, "fasta")
        for record in gRNA_fasta_list:
            gRNA_list += [str(record.seq).upper()]            
    except ValueError:    #to catch fasta files with single gRNA (use SeqIO.read instead of SeqIO.parse)
        gRNA_fasta_list = SeqIO.read(gRNA_fasta, "fasta")
        gRNA_list += [str(gRNA_fasta_list.seq.upper())] 
    return gRNA_list
              

## Find all gRNA + pam sequences
def gRNA_pam_listing(gRNA_list, pam, genomic_template):
    # Time and space depends on Regex's re.findall
    reg_pam = ''
    F = []
    R = []
    for letter in pam:
        reg_pam += IUPAC_dict[letter]
    
    gRNA_pam_list = [x + reg_pam for x in gRNA_list]
    
    for gRNA_pam in gRNA_pam_list:
        F += re.findall(gRNA_pam, genomic_template)
        R += re.findall(gRNA_pam, str(Seq(genomic_template).reverse_complement()))
        
    return [list(set(F)), list(set(R))]

## Primer design

#  Create error class
class SequenceError(Exception):
    pass

#  Find gRNA in genomic DNA: returns list of positions of each gRNA's cleavage site; only one position per gRNA will be returned
def gRNA_finder(gRNA_pam_list, genomic_template):
    # Time: Depends on Regex's re.finditer
    # Space: O(n) where n is len(gRNA_pam_list)
    
    temp = []
    lst_grand = []
    unpacked_lst = []
    
    #note: data structure of gRNA_pam_list is [ [forward gRNA + pam], [reverse complement gRNA + pam]]
    
    #find position of forward gRNA cleavage site
    for gRNA_pam in gRNA_pam_list[0]:
        #find overlapping matches for forward gRNA AND keep list of Cas9 cleavage site position
        temp  += [ (m.start() + len(gRNA_pam) - len(pam) - cleavage_pos)
                   for m in re.finditer('(?=' + gRNA_pam + ')', genomic_template) ] 
        lst_grand += [(gRNA_pam, temp)]
        temp = []
        
    #find position of reverse complement gRNA cleavage site
    for gRNA_pam in gRNA_pam_list[1]:
        #find overlapping matches with reverse comp gRNA AND keep list Cas9 cleavage site position
        temp += [ (m.start() + len(pam) + cleavage_pos)
                  for m in re.finditer('(?=' + str(Seq(gRNA_pam).reverse_complement()) + ')', genomic_template) ] 
        lst_grand += [(gRNA_pam, temp)]
        temp = []
    
    #structure of lst_grand is e.g. [('GGCGTTGACAGATGAGGGGCAGG', [403, 501]), ('AATGCTGGATTTTCTGCCTGTGG', [643])]
    for i in lst_grand:
        for j in i[1]:
            unpacked_lst += [tuple((i[0], j))]
    
    #reset list to save space
    lst_grand = []
    
    return unpacked_lst

# Find tuple of primers F and R for ONE cleavage site
## returns a CollapsedPrimer3 obj
def primer_design(pos, template, amplicon_size, primer_task = "generic"):
    
    #find targeted regions where primers will be designed
    half_len = round(amplicon_size[1]/2) #half the length of amplicon
    region = template[pos - half_len: pos + half_len]        
    
    #design primers for current position
    try:
        primer3_params = [
            {
                'SEQUENCE_TEMPLATE': region
            },
            {
                'PRIMER_OPT_SIZE': primer_opt_size,
                'PRIMER_MIN_SIZE': primer_min_size,
                'PRIMER_MAX_SIZE': primer_max_size,
                'PRIMER_OPT_TM': primer_opt_tm,
                'PRIMER_MIN_TM': primer_min_tm,
                'PRIMER_MAX_TM': primer_max_tm,
                'PRIMER_MIN_GC': primer_min_GC,
                'PRIMER_MAX_GC': primer_max_GC,

                'PRIMER_DNA_CONC': primer_dna_conc,
                'PRIMER_PRODUCT_SIZE_RANGE': amplicon_size,
                'PRIMER_DNTP_CONC': primer_dNTP_conc,
                'PRIMER_PAIR_MAX_DIFF_TM': primer_pair_max_diff_TM,
                'PRIMER_SALT_DIVALENT': primer_salt_divalent,
                'PRIMER_NUM_RETURN': primer_num_return,
                
                'PRIMER_TASK': primer_task
                
            }
        ]
        design_primer = primer3.bindings.design_primers(*primer3_params)
        
    except(OSError) as e:
        print("Error: PCR product size parameter demand primers to be designed beyond genomic template. Include longer genomic template sequence or 'decrease primer_product_size_range'.")
        raise e
    
    ## retrieve all primers, collapse to unique sequences
    primer3_nr = CollapsedPrimer3(design_primer)
    primer3_nr.params = primer3_params
    primer3_nr.target = half_len
    
    if primer3_nr.raw["PRIMER_PAIR_NUM_RETURNED"] < primer_num_return:
        print(f"User-defined number of forward primers to return = {primer_num_return} is not met. Explanation from Primer3:")
        print("   Pair: " + design_primer['PRIMER_PAIR_EXPLAIN'])
        print("   Left: " + design_primer['PRIMER_LEFT_EXPLAIN'])
        print("   Right: " + design_primer['PRIMER_RIGHT_EXPLAIN'])
    
    # Return 
    return primer3_nr

## returns dict of {<subject_id>: (<start>, <end>)} [range is 0-indexed, end-exclusive]
## representing position of exact match(es)
def template_in_background(template_fasta, background_fasta):
    fout = "template-in-background.blastn.tsv"
    ## BLAST set-up
    custom_fields = ["qseqid", "sseqid", "qstart", "qend", "sstart", "send", "qlen", "gaps", "mismatch"]
    args = list(map(str, [
        "blastn",
        "-query", template_fasta,
        "-subject", background_fasta,
        "-word_size", secondary_PCR_amplicon_size[0],
        "-perc_identity", 100,                     ## exact matches only
        "-ungapped",                               ## ungapped because we're looking for exact matches
        "-dust", "no",
        "-soft_masking", "false",
        "-out", fout,
        "-outfmt", "6 " + ' '.join(custom_fields)
    ]))
    ## execute
    subprocess.run(args, check = True)
    ## filter
    qresult = SearchIO.parse(fout, "blast-tab", fields = custom_fields)
    to_mask = {}
    for query in qresult:
        for hit in query.hits:
            for hsp in hit.hsps:
                ## filter for alignments that span the full length of the query sequence
                ## (pident and gapped have been accounted for in args, but filter again here for insurance)
                if (hsp.query_start == 0 and hsp.query_end == query.seq_len
                    and hsp.gap_num == 0 and hsp.mismatch_num == 0):
                    ## no need to check if hit_start < hit_end; hit_start and hit_end are sorted when using blast-tab
                    to_mask[hit.id] = to_mask.get(hit.id, set()).union({(hsp.hit_start, hsp.hit_end)})
    os.remove(fout)
    return to_mask

## returns dict of {<subject_id>: (<start>, <end>)} [range is 0-indexed, end-exclusive]
## representing position of exact match(es)
def sequence_in_fasta(query_fasta, subject_fasta):
    fout = "query-in-fasta.blastn.tsv"
    ## get length of shortest query sequence
    min_query_len = min(len(record.seq) for record in SeqIO.parse(query_fasta, "fasta"))
    ## BLAST set-up
    custom_fields = ["qseqid", "sseqid", "qstart", "qend", "sstart", "send", "qlen", "gaps", "mismatch"]
    args = list(map(str, [
        "blastn",
        "-query", query_fasta,
        "-subject", subject_fasta,
        "-word_size", min_query_len,
        "-perc_identity", 100,                     ## exact matches only
        "-ungapped",                               ## ungapped because we're looking for exact matches
        "-dust", "no",
        "-soft_masking", "false",
        "-out", fout,
        "-outfmt", "6 " + ' '.join(custom_fields)
    ]))
    ## execute
    subprocess.run(args, check = True)
    ## filter
    qresult = SearchIO.parse(fout, "blast-tab", fields = custom_fields)
    to_mask = {}
    for query in qresult:
        for hit in query.hits:
            for hsp in hit.hsps:
                ## filter for alignments that span the full length of the query sequence
                ## (pident and gapped have been accounted for in args, but filter again here for insurance)
                if (hsp.query_start == 0 and hsp.query_end == query.seq_len
                    and hsp.gap_num == 0 and hsp.mismatch_num == 0):
                    ## no need to check if hit_start < hit_end; hit_start and hit_end are sorted when using blast-tab
                    to_mask[hit.id] = to_mask.get(hit.id, set()).union({(hsp.hit_start, hsp.hit_end)})
    os.remove(fout)
    return to_mask

def num_exact_in_fasta(seq, subject_fasta):
    fp, query_fasta = tempfile.mkstemp()
    seq_records = [SeqRecord(Seq(seq), id = "query", description = "query")]
    SeqIO.write(seq_records, query_fasta, "fasta")
    exact_matches = sequence_in_fasta(query_fasta, subject_fasta)
    os.remove(query_fasta)
    num_exact_matches = sum(len(coords) for coords in exact_matches.values())
    return num_exact_matches

## blasts a sequence against a fasta file and checks if it only appears once (match across full length, no mismatches, no indels)
def unique_in_fasta(seq, fasta):
    num_exact_matches = num_exact_in_fasta(seq, fasta)
    if num_exact_matches == 0: return None
    return num_exact_in_fasta(seq, fasta) == 1

# Filter primer_list (primer_list = ((list of F primers), (list of R primers)) by pseudo_primer_blast 
# Returns 2 elements:
# [1] list of uniquely binding primers (str of sequences)
# [2] set of tuples of disallowed primer pairs (str of sequences)
## previously: pseudo_primer_blast
## previously: invalid_primer_pairs
## exclude should be a CollapsedSeqs object. All sequences in 'exclude' will be excluded from BLAST.
## merge should be a PrimerOffTargetChecker object. All bindings in 'merge' will be added.
def make_off_target_checker(collapsed_primer3, template_DNA_fasta, background_DNA_fasta,
                            min_acceptable_off_target_amplicon_size,
                            exclude_collapsed_seqs = None, merge_primers_off_target_checker = None,
                            intersect_primers_off_target_checker = None,
                            prefix = None):
    ## write unique sequences to file in prep for blast
    uniq_primers_fasta = "primerg_specificity_check_primers.fasta"
    collapsed_primer3.write_fasta_uniq(uniq_primers_fasta, prefix = prefix, exclude = exclude_collapsed_seqs)
    
    ## get exact position of template_DNA_fasta in background_DNA_fasta (so we can add to on_target_ranges)
    ## this is a dict {<chrom>: (<start>, <end>)} that we can use to extend on_target_ranges
    to_mask = template_in_background(template_DNA_fasta, background_DNA_fasta)
    
    ## instantiate output object PrimersOffTargetChecker to track chrom & position of hits
    ## allows filtering to allow primer pairs with amplicon >= min_acceptable_off_target_amplicon_size
    binding_primers = PrimersOffTargetChecker(uniq_primers_fasta)
    _num_primers_for_screening = len(binding_primers.primers_map)
    binding_primers.forgive_coords = {subject: set().union(to_mask.get(subject, set()),
                                                           on_target_ranges.get(subject, set()))
                                      for subject in set(list(to_mask.keys()) + list(on_target_ranges.keys()))}
    binding_primers.template_coords = to_mask
    
    ## merge
    if merge_primers_off_target_checker: binding_primers.merge(merge_primers_off_target_checker)
    ## intersect
    if intersect_primers_off_target_checker:
        shared_seqs = [seq for seq in intersect_primers_off_target_checker.seqs()
                       if collapsed_primer3.has_seq(seq)]
        if shared_seqs:
            binding_primers.merge(intersect_primers_off_target_checker.subset(seqs = set(shared_seqs)))
    
    ## return binding_primers immediately if no new primers for screening
    if _num_primers_for_screening == 0:
        print("No new primer sequences. Skipping BLAST.")
        return binding_primers

    # BLAST set-up (NcbiblastnCommandline is/will be deprecated, use subprocess instead
    ## we'll use outfmt 6 instead with custom fields to reduce file size
    fout_blast = "primer-blast.tsv"
    custom_fields = ["qseqid", "sseqid", "qstart", "qend", "sstart", "send", "qlen", "qframe", "sframe", "btop"]
    args = list(map(str, [
        "blastn",
        "-query", uniq_primers_fasta,
        "-subject", background_DNA_fasta,
        "-evalue", evalue,                #high e-value to allow more chance hits so we can evaluate all possible binding sites
        "-word_size", word_size,          #small seed is required as primer seq is short; allows more possible hits; default=1000
        "-task", "blastn-short",          #BLASTN program optimized for sequences shorter than 50 bases; default=7
        "-dust", "no",                    #No masking of low complexity sequence
        "-soft_masking", "false",         #No masking of low complexity sequence
        # "-strand", "plus",                #Only plus strand is used as query ## this doesn't do what wy thought it did (+ of SUBJECT)
        "-max_target_seqs", max_target_seqs,           #Allow more hits to be shown; default=50000
        "-out", fout_blast,               #store BLAST result in tab file
        "-outfmt", "6 " + ' '.join(custom_fields)      #format of output as tsv means outfmt = 6
    ]))
    
    # Run BLAST
    print("--Executing BLASTN--", end='\r')
    subprocess.run(args, check = True)
    print('', end = '\r')
    
    # Read BLAST results
    blast_records = BlastResult(fout_blast, "blast-tab", fields = custom_fields)
    
    ## Iterate through BLAST results and add non-specific primers to PrimersOffTargetChecker object
    ## Pseudo-primer-blast checks if primer binds to site if they meet the requirements below
    ## (parameters adjustable as variables above):
    
    ## [1] Length of alignment must exceed 10
    ## [2] Total number of mismatches does NOT exceed 6
    ## [3] Last 5 bp of 3' end of aligned sequences has at least 3 matches
    
    for i, hsp in enumerate(blast_records):
        sbtop = hsp.sbtop
        if (hsp.query.seq_len >= min_primer_len
            and hsp.num_mismatch() + hsp.num_unaligned() <= total_mismatch):
            # ## TO COMMENT OUT
            # if (binding_primers.primers_map[hsp.query.id]
            #     in {'GCCTCCACACGCTCTTTCGT', 'ACGCCTCCACACGCTCTTTCGT', 'ACGCCTCCACACGCTCTTT',
            #         'ACGCCTCCACACGCTCTTTC', 'TCGCATGGCCTTTGTTCGGT'}):
            #     print("highlight:", binding_primers.primers_map[hsp.query.id],
            #           hsp.hsp.hit_frame, hsp.hsp.query_frame, hsp.btop, hsp.sbtop)
            if hsp.sbtop[-5:].count('.') >= three_prime_match:
                binding_primers.add_primer(hsp.query.id, hsp.subject.id,
                                           (1 if hsp.hsp.hit_frame == hsp.hsp.query_frame else -1),
                                           hsp.hsp.hit_start, hsp.hsp.hit_end)
    
    # ## read uniq_primers_fasta for mapping primer ids to sequences
    # uniq_primers_dict = {record.id: str(record.seq) for record in SeqIO.parse(uniq_primers_fasta)}
    # ## set of tuples of disallowed primer pairs (str of primer id)
    # ## note that on_target_ranges MUST be provided otherwise all primers will be disqualified
    # disallowed_pairs = binding_primers.invalid_primer_pairs(max_off_target_amplicon_size = max_off_target_amplicon_size,
    #                                                         forgive_coords = merged_mask_ranges)
    # disallowed_pairs_seq = set((uniq_primers_dict[plus_primer], uniq_primers_dict[minus_primer])
    #                            for plus_primer, minus_primer in disallowed_pairs)
    # ## get uniquely-binding primers
    # uniq_binding_primers = binding_primers.uniquely_binding_primers()
    # uniq_binding_seq = [uniq_primers_dict[primer_id] for primer_id in uniq_binding_primers]
    # ## return disallowed primer pair sequences
    return binding_primers


# Generate sets of primers given primer_list ((list of F primers), (list of R primers))
# N is the number of primer sets desired
def generate_primer_pairs(primer_list, N, template, amplicon_size):
    
    original_N = N # keeps track of starting N value
    n = N          # keeps track of alternative list length

    # Generate F and R primer lists where positions in each list is corresponding
    # e.g. lst[0][0] corresponds to list[1][0]
    F_lst = []     
    R_lst = []
    alt_F_lst = [] # alternate lists keep specific + unspecific primer pairs
    alt_R_lst = []
    
    local_template = template[int(primer_list[0]-(0.5*amplicon_size[1])):
                              int(primer_list[0]+(0.5*amplicon_size[1])) ]
    
    for i in range(1,3):
        if i == 1: # Means we are looking at forward (left) primer
            if N <= 0:
                break
            
            for primer in primer_list[i]:
                if N <= 0:
                    break
                
                try:

                    design_primer = primer3.bindings.design_primers(
                    {
                        'SEQUENCE_TEMPLATE': local_template
                    },
                    {
                        'PRIMER_OPT_SIZE': primer_opt_size,
                        'PRIMER_MIN_SIZE': primer_min_size,
                        'PRIMER_MAX_SIZE': primer_max_size,
                        'PRIMER_OPT_TM': primer_opt_tm,
                        'PRIMER_MIN_TM': primer_min_tm,
                        'PRIMER_MAX_TM': primer_max_tm,
                        'PRIMER_MIN_GC': primer_min_GC,
                        'PRIMER_MAX_GC': primer_max_GC,
        
                        'PRIMER_DNA_CONC': primer_dna_conc,
                        'PRIMER_PRODUCT_SIZE_RANGE': amplicon_size,
                        'PRIMER_DNTP_CONC': primer_dNTP_conc,
                        'PRIMER_PAIR_MAX_DIFF_TM': primer_pair_max_diff_TM,
                        'PRIMER_SALT_DIVALENT': primer_salt_divalent,
                        'PRIMER_NUM_RETURN': primer_num_return,
                        
                        'SEQUENCE_PRIMER': primer
                        
                    })
                
                except OSError as e:
                    print("Error: PCR product size parameter demand primers to be designed beyond genomic template. Please provide longer genomic template sequence or decrease 'primer_product_size_range'.")
                    raise e
                
                except IndexError as e:
                    print("Error: local_template cannot be generated because cleavage site position +/- 0.5*desired_amplicon_size is beyond user-supplied genomic template. Please provide longer genomic template sequence.")
                    raise e
                    
                except: # List is empty
                    pass
        
                
                for j in range(len(design_primer)):
                    if N <= 0:
                        break
                    try:
                        if design_primer['PRIMER_RIGHT_' + str(j) + '_SEQUENCE'] in primer_list[2]:
                            F_lst += [primer]
                            R_lst += [design_primer['PRIMER_RIGHT_' + str(j) + '_SEQUENCE']]
                            N -= 1
                            break # ensures we only use the same primer only one time
                        else:
                            if n > 0:
                                alt_F_lst += [primer]
                                alt_R_lst += [design_primer['PRIMER_RIGHT_' + str(j) + '_SEQUENCE']] 
                                n -= 1
                                break
           
                    except KeyError:
                        break

        elif i == 2: # Means we are looking at reverse (right) primer
            if N <= 0:
                break
            
            for primer in primer_list[i]:
                if N <= 0:
                    break
                
                try:
                    design_primer = primer3.bindings.design_primers(
                    {
                        'SEQUENCE_TEMPLATE': local_template
                    },
                    {
                        'PRIMER_OPT_SIZE': primer_opt_size,
                        'PRIMER_MIN_SIZE': primer_min_size,
                        'PRIMER_MAX_SIZE': primer_max_size,
                        'PRIMER_OPT_TM': primer_opt_tm,
                        'PRIMER_MIN_TM': primer_min_tm,
                        'PRIMER_MAX_TM': primer_max_tm,
                        'PRIMER_MIN_GC': primer_min_GC,
                        'PRIMER_MAX_GC': primer_max_GC,
        
                        'PRIMER_DNA_CONC': primer_dna_conc,
                        'PRIMER_PRODUCT_SIZE_RANGE': amplicon_size,
                        'PRIMER_DNTP_CONC': primer_dNTP_conc,
                        'PRIMER_PAIR_MAX_DIFF_TM': primer_pair_max_diff_TM,
                        'PRIMER_SALT_DIVALENT': primer_salt_divalent,
                        'PRIMER_NUM_RETURN': primer_num_return,
                        
                        'SEQUENCE_PRIMER_REVCOMP': primer
                        
                    })
                
                except(OSError):
                    print("Error: PCR product size parameter demand primers to be designed beyond genomic template. Include longer genomic template sequence or 'decrease primer_product_size_range'.")
                    exit()
                    
                except: # List is empty
                    pass
                
                for j in range(len(design_primer)):
                    if N <= 0:
                        break
                    try:
                        if design_primer['PRIMER_LEFT_' + str(j) + '_SEQUENCE'] in primer_list[1]:
                            F_lst += [design_primer['PRIMER_LEFT_' + str(j) + '_SEQUENCE']]
                            R_lst += [primer]
                            N -= 1
                            break # ensures we only use the same primer only one time
                        else:
                            if n > 0:
                                alt_F_lst += [design_primer['PRIMER_LEFT_' + str(j) + '_SEQUENCE']]
                                alt_R_lst += [primer]
                                n -= 1
                                break # ensures we only use the same primer only one time
           
                    except KeyError:
                        break                
                
                
    if N > 0: # when specific F cannot correspond to specific R
        try:
            F_lst += alt_F_lst[0: original_N - (original_N-N)]
            R_lst += alt_R_lst[0: original_N - (original_N-N)]
        except: # in case alt_F or alt_R is empty, which is unlikely
            pass
     
    # Returns [pos, [forward primers], [reverse primers]], where position of pairing primers are corresponding
    # e.g. list[0][0] corresponds to list[1][0]
    return [primer_list[0], F_lst, R_lst]

## previously: Generate sets of primers given primer_list ((list of F primers), (list of R primers))
## current: primers_nr is an object of class CollapsedPrimer3
## note that this function ADDS new primers, doesn't overwrite existing primers in primers_nr, unless overwrite=True
def primer_design_for_uniquely_binding_primers(pos, primers_nr, template, amplicon_size,
                                               unique_in_template = False, overwrite = False):
    
    ## uniquely binding primers that we can try to find new pairs for
    if unique_in_template:
        uniquely_binding_primers = primers_nr.off_target_checker.uniquely_binding_primers_in_template(seq = True)
    else:
        uniquely_binding_primers = primers_nr.off_target_checker.uniquely_binding_primers(seq = True)
    if not uniquely_binding_primers: return primers_nr
    
    local_template = template[int(pos-(0.5*amplicon_size[1])):
                              int(pos+(0.5*amplicon_size[1])) ]
    
    # print("# primers premerge left:", len(primers_nr.primers_left_unique()),
    #       "; # primers premerge right:", len(primers_nr.primers_right_unique()),
    #       "; # primers premerge total:",
    #       len(primers_nr.primers_left_unique()) + len(primers_nr.primers_right_unique()))
    
    ## new object for merging CollapsedPrimer3 objects
    if overwrite:
        merged_primers_nr = CollapsedPrimer3({"PRIMER_LEFT": [], "PRIMER_RIGHT": []})
        merged_primers_nr.params = primers_nr.params
        merged_primers_nr.target = primers_nr.target
    else:
        merged_primers_nr = copy.deepcopy(primers_nr)
        
    ## uniq_primer_i is an integer prepended to primers
    ## unique for each execution of primer3.bindings.design_primers to facilitate merging primers into merged_primers_nr
    uniq_primer_i = 0
    for primer_seq in uniquely_binding_primers:
        ## if a single primer sequence can be found on both + and - strands, we design_primers for both situations
        primer_strands = primers_nr.get_strand_by_seq(primer_seq, return_type = 'c')
        for primer_strand in primer_strands:
            try:
                design_primer = primer3.bindings.design_primers(
                {
                    'SEQUENCE_TEMPLATE': local_template
                },
                {
                    'PRIMER_OPT_SIZE': primer_opt_size,
                    'PRIMER_MIN_SIZE': primer_min_size,
                    'PRIMER_MAX_SIZE': primer_max_size,
                    'PRIMER_OPT_TM': primer_opt_tm,
                    'PRIMER_MIN_TM': primer_min_tm,
                    'PRIMER_MAX_TM': primer_max_tm,
                    'PRIMER_MIN_GC': primer_min_GC,
                    'PRIMER_MAX_GC': primer_max_GC,

                    'PRIMER_DNA_CONC': primer_dna_conc,
                    'PRIMER_PRODUCT_SIZE_RANGE': amplicon_size,
                    'PRIMER_DNTP_CONC': primer_dNTP_conc,
                    'PRIMER_PAIR_MAX_DIFF_TM': primer_pair_max_diff_TM,
                    'PRIMER_SALT_DIVALENT': primer_salt_divalent,
                    'PRIMER_NUM_RETURN': primer_num_return_for_unique,

                    ('SEQUENCE_PRIMER' if primer_strand == '+' else 'SEQUENCE_PRIMER_REVCOMP'): primer_seq

                })
                
            except OSError as e:
                print("Error: PCR product size parameter demand primers to be designed beyond genomic template. Please provide longer genomic template sequence or decrease 'primer_product_size_range'.")
                print("Primer sequence:", primer_seq)
                print("Primer strand:", primer_strand)
                raise e
            except IndexError as e:
                print("Error: local_template cannot be generated because cleavage site position +/- 0.5*desired_amplicon_size is beyond user-supplied genomic template. Please provide longer genomic template sequence.")
                raise e
            
            ## retrieve all primers, collapse to unique sequences
            primer3_nr = CollapsedPrimer3(design_primer)
            ## merge with merged_primers_nr
            for primer_id in primer3_nr.primers_left:
                primer = primer3_nr.primers[primer_id]
                primer.seqid = f"{uniq_primer_i}|{primer_id}"
                merged_primers_nr.add_primer(primer, 1)
            for primer_id in primer3_nr.primers_right:
                primer = primer3_nr.primers[primer_id]
                primer.seqid = f"{uniq_primer_i}|{primer_id}"
                merged_primers_nr.add_primer(primer, -1)
            ## update index for uniq_primer_i
            uniq_primer_i += 1
            
    return merged_primers_nr

# Generates PCR amplicon made by F and R primer (without primer sequence within the returned amplicon)
def generate_amplicon(pos, F, R, template, amplicon_size):
    local_template = template[int(pos-(0.5*amplicon_size[1])): int(pos+(0.5*amplicon_size[1])) ]
    pos_start = str(local_template).find(str(F)) + len(F)
    pos_end = str(local_template).find(str(Seq(R).reverse_complement())) 
    return local_template[pos_start: pos_end]

## collapse primer pairs to unique sequence combinations AND sort in order of increasing penalty
## primer ids will be arbitrary selected for sequences associated with >1 seqids
## outputs list of PrimerPair objects
def sort_unique_primer_pairs(primer_pairs):
    unique_primer_pairs = {}
    ## iterate through all pairs and get unique sequence combos
    for primer_pair in primer_pairs:
        ## get sequence combo
        primer_pair_seq = primer_pair.as_seq_tuple()
        ## ignore of sequence combo already seen
        if primer_pair_seq in unique_primer_pairs: continue
        ## add if not seen before
        else: unique_primer_pairs[primer_pair_seq] = primer_pair
    ## sort
    sorted_primer_pairs = sorted(unique_primer_pairs.values(), key = lambda primer_pair:primer_pair.penalty())
    return sorted_primer_pairs
    
#######################
##   FUNCTIONS FOR   ##
##  PRIMERG_DC_WRAP  ##
#######################

def _gRNA_pam_pos_sliding(gRNA_pam_list, subject, window = 100000, overlap = 100):
    """
    Searches for regex pattern in memory saving way using sliding
    within ``window`` and overlap ``overlap``.
    Outputs dict of 2 elements: 'F' indexed pertains to gRNA-PAM in fwd strand, 'R' to rvs.
    Each element is a list of tuples of [(start1, end1), (start2, end2), ...].
    Positions pertaining to rvs strand matches are relative to fwd strand (i.e. slice fwd then reverse complement to get gRNA-PAM).
    """
    # Time and space depends on Regex's re.findall
    F, R = [], []
    for i in range(0, len(subject), window - overlap):
        fwd = subject[i:i+window]
        rvs = str(fwd.reverse_complement()).upper()
        fwd = str(fwd).upper()
        window_truesize = len(fwd)
        for gRNA_pam in gRNA_pam_list:
            # F.extend(list(re.findall(gRNA_pam, fwd, overlapped = True)))
            # R.extend(list(re.findall(gRNA_pam, rvs, overlapped = True)))
            F.extend([(i+x.start(), i+x.end()) for x in re.finditer(gRNA_pam, fwd, overlapped = True)])
            # R.extend([len(subject)-i-window_truesize+x.start()
            #           for x in re.finditer(gRNA_pam, rvs, overlapped = True)]) ## relative to rvs strand
            R.extend([(i+window_truesize-x.end(), i+window_truesize-x.start())
                       for x in re.finditer(gRNA_pam, rvs, overlapped = True)]) ## relative to fwd strand
    return {"F": sorted(set(F)), "R": sorted(set(R))}

def gRNA_pam_pos_in_IndexedFasta(gRNA_list, pam, ifasta, window = 100000, overlap = 100):
    """
    Searches for regex pattern in memory saving way using sliding
    within ``window`` and overlap ``overlap``.
    Outputs dict of N elements of dict of 2 elements which are lists: {<chrom>: {'F': [], 'R': []}}
    where 'F' indexed pertains to gRNA-PAM in fwd strand, 'R' to rvs.
    Each list has the format [(start1, end1), (start2, end2), ...].
    Positions pertaining to rvs strand matches are relative to fwd strand (i.e. slice fwd then reverse complement to get gRNA-PAM).
    """
    # Time and space depends on Regex's re.findall
    reg_pam = ''
    ## make PAM pattern
    for letter in pam:
        reg_pam += IUPAC_dict[letter]
    ## make gRNA+PAM pattern
    gRNA_pam_list = [x + reg_pam for x in gRNA_list]
    ## find gRNA+PAM
    output = {}
    for seq in ifasta:
        fr = _gRNA_pam_pos_sliding(gRNA_pam_list, seq, window = window, overlap = overlap)
        output[seq.name] = fr
    return output

# ## group positions if they are within distance from each other
# ## F and R are iterables of integers (position)
# def group_cleavage_sites_within_distance(F, R, distance = 1000):
#     sorted_F = sorted(F)
#     sorted_R = sorted(R)
#     output = []
#     last_pos = min(F+R)
#     group = {"F": [], "R": []}
#     ## while both sorted_F and sorted_R are not empty
#     while sorted_F and sorted_R:
#         if sorted_F[0] <= sorted_R[0]:
#             ## if too far apart, close current group and start new group
#             if (sorted_F[0] - last_pos) >= distance:
#                 output.append(group)
#                 group = {"F": [], "R": []}
#             ## add current F pos to new group
#             last_pos = sorted_F.pop(0)
#             group["F"].append(last_pos)
#         else:
#             ## if too far apart, close current group and start new group
#             if (sorted_R[0] - last_pos) >= distance:
#                 output.append(group)
#                 group = {"F": [], "R": []}
#             ## add current R pos to new group
#             last_pos = sorted_R.pop(0)
#             group["R"].append(last_pos)
#     ## if only one of sorted_F and sorted_R is not empty
#     if sorted_F or sorted_R:
#         leftover_strand = "F" if sorted_F else "R"
#         leftover_pos = sorted_F if sorted_F else sorted_R
#         while leftover_pos:
#             ## if too far apart, close current group and start new group
#             if (leftover_pos[0] - last_pos) >= distance:
#                 output.append(group)
#                 group = {"F": [], "R": []}
#             ## add current R pos to new group
#             last_pos = leftover_pos.pop(0)
#             group[leftover_strand].append(last_pos)
#     output.append(group)
#     return output

## group positions if they are within distance from each other
## F and R are iterables of (<gRNA+PAM seq>, <int (position)>)
def group_seq_cleavage_pos_within_distance(F, R, distance = 1000):
    if not F and not R:
        return []
    get_pos = lambda x:x[1]
    sorted_F = sorted(F, key = get_pos)
    sorted_R = sorted(R, key = get_pos)
    output = []
    last_pos = get_pos(min(F+R, key = get_pos))
    group = {"F": [], "R": []}
    ## while both sorted_F and sorted_R are not empty
    while sorted_F and sorted_R:
        pos_F, pos_R = get_pos(sorted_F[0]), get_pos(sorted_R[0])
        if pos_F <= pos_R:
            ## if too far apart, close current group and start new group
            if (pos_F - last_pos) >= distance:
                output.append(group)
                group = {"F": [], "R": []}
            ## add current F pos to new group
            last_pos = pos_F
            group["F"].append(sorted_F.pop(0))
        else:
            ## if too far apart, close current group and start new group
            if (pos_R - last_pos) >= distance:
                output.append(group)
                group = {"F": [], "R": []}
            ## add current R pos to new group
            last_pos = pos_R
            group["R"].append(sorted_R.pop(0))
    ## if only one of sorted_F and sorted_R is not empty
    if sorted_F or sorted_R:
        leftover_strand = "F" if sorted_F else "R"
        leftover_pos = sorted_F if sorted_F else sorted_R
        while leftover_pos:
            ## if too far apart, close current group and start new group
            if (get_pos(leftover_pos[0]) - last_pos) >= distance:
                output.append(group)
                group = {"F": [], "R": []}
            ## add current R pos to new group
            last_pos = get_pos(leftover_pos[0])
            group[leftover_strand].append(leftover_pos.pop(0))
    output.append(group)
    return output

## find cleavage positions from output of gRNA_pam_pos_in_IndexedFasta
def gRNA_pam_pos_to_cleavage_pos(gRNA_pam_pos_FR):
    """
    Takes output of gRNA_pam_pos_in_IndexedFasta.
    Outputs dict of 2 elements: 'F' indexed pertains to gRNA-PAM in fwd strand, 'R' to rvs.
    Each element is a dictionary of lists of cleavage positions of {<chr1>: [pos1, pos2, ...], <chr2>: [pos1, pos2, ...], ...}.
    Positions pertaining to rvs strand matches are relative to fwd strand (i.e. slice fwd then reverse complement to get gRNA-PAM).
    cleavage_pos assumed to be 1-indexed.
    """
    F, R = {}, {}
    pam_len = len(pam)
    ## fwd
    for seq_name, ranges in gRNA_pam_pos_FR['F'].items():
        F[seq_name] = [end - pam_len - cleavage_pos for start, end in ranges] ## no -1 cuz end is exclusive
    ## rvs
    for seq_name, ranges in gRNA_pam_pos_FR['R'].items():
        R[seq_name] = [start + pam_len + cleavage_pos - 1 for start, end in ranges] ## -1 because start is inclusive
    return {'F': F, 'R': R}

## find cleavage positions from output of gRNA_pam_pos_in_IndexedFasta
def gRNA_pam_pos_to_seq_and_cleavage_pos(gRNA_pam_pos_FR, ifasta, upper = True, as_string = True):
    """
    Takes output of gRNA_pam_pos_in_IndexedFasta.
    Outputs dict of N elements corresponding to the number of sequences in ifasta.
    Each element is a dict of 2 element: 'F' indexed pertains to gRNA-PAM in fwd strand, 'R' to rvs.
    Each element of 'F'/'R' dict is a lists of gRNA+PAM seq and cleavage positions of
      [(<gRNA+PAM 1>, pos1), (<gRNA+PAM 2>, pos2), ...].
    Positions pertaining to rvs strand matches are relative to fwd strand (i.e. slice fwd then reverse complement to get gRNA-PAM).
    cleavage_pos assumed to be 1-indexed.
    """
    pam_len = len(pam)
    output = {}
    for seq_name, seq_dat in gRNA_pam_pos_FR.items():
        output[seq_name] = {}
        for strand, ranges in seq_dat.items():
            if strand in {1, '+', 'F', "plus", "Plus"}:
                new_seq_pos = [(ifasta[seq_name][start:end], end - pam_len - cleavage_pos)
                               for start, end in ranges]
            else:
                new_seq_pos = [(ifasta[seq_name][start:end].reverse_complement(), start + pam_len + cleavage_pos)
                               for start, end in ranges]
            if upper:
                new_seq_pos = [(grna_pam_seq.upper(), pos) for grna_pam_seq, pos in new_seq_pos]
            if as_string:
                new_seq_pos = [(str(grna_pam_seq), pos) for grna_pam_seq, pos in new_seq_pos]
            output[seq_name][strand] = new_seq_pos
    return output
