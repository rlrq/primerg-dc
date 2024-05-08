## PRIMERg
# Deep-sequencing primer design (with primerg) for multiple targets on a continuous genome

## Import packages
from primerg_header import *

# Test
#Input:
#x = "aatttcagaaatctgtttttttttctctctctatcttcttccgcgtaatgattgaaaacctttttttcaatctatgtgtaattttattaaaataagaattaataataatctttttatttcgtttttttgagatctatgtttcatgaattttttgcataa"
#F = "tcttcttccgcgtaa"
#R = "caaaaaattcatga"
#generate_amplicon(F, R, x)

#Output: 'tgattgaaaacctttttttcaatctatgtgtaattttattaaaataagaattaataataatctttttatttcgtttttttgagatctatgtt'


### Execution

## there should only be ONE sequence in template_DNA_fasta file
def primerg(gRNA_pos, template_DNA_fasta, suppress_output_message = False, print_df = True):
    ## Set up dataframe
    df = pd.DataFrame(columns = ['cleavage_pos', 'gRNA_pam', 
                                 'F1', 'R1', 
                                 'F2', 'R2',
                                 'unique_F2', 'unique_R2'])
    
    ## Ensure excel file is not opened
    if output_excel:
        try: 
            df.to_excel(output_excel, index = False)
        except PermissionError as e:
            print(f"Please close {output_excel} and try again")
            raise e
    
    ## Ensure tsv file is not opened
    if output_tsv:
        try: 
            df.to_csv(output_tsv, index = False)
        except PermissionError as e:
            print(f"Please close {output_tsv} and try again")
            raise e
    
    # Initialize empty variables
    amplicon = None
    filtered_primer_list = None
    genomic_template = str(SeqIO.read(template_DNA_fasta, "fasta").seq).upper() ## there should only be ONE sequence in this file
    
    # # Obtain list of gRNAs and cleavage site positions
    # gRNA_pos = gRNA_finder(gRNA_pam_list, genomic_template)
    
    index = 0
    start_index = 0
    original_value = primary_PCR_amplicon_size[1]

    for pos in gRNA_pos:

        ## TODO: comment out after we've solved the issue of primers that are non-specific in template
        ## there should be no specific F/R1 primers for this pos
        # if pos[1] not in {6384, 10236}:
        # if pos[1] != 10236:
        # if pos[1] != 6384:
        # if pos[1] != 20099:
        # if pos[1] != 6214:
        #     continue

        start_index = index
        # Update meta-data
        df2 = pd.DataFrame([ [pos[1], pos[0], '', '', '', '', '', ''] ,],
                           columns = ['cleavage_pos', 'gRNA_pam', 'F1', 'R1', 'F2', 'R2', 'unique_F2', 'unique_R2'])

        for i in range(num_primer_pair_per_cleavage_pos):
            df = pd.concat([df, df2], ignore_index = True)

        print(f"--Processing gRNA + PAM {df.loc[index, 'gRNA_pam']} at position {df.loc[index, 'cleavage_pos']}")

        ###################
        ##  PRIMARY PCR  ##
        ###################
        
        print("[[ Primary amplicon ]]")
        
        # Update primary PCR primer (first PCR)
        ## previously: output stored in variable 'primer_list', with structure [pos, <left primer seqs>, <right primer seqs>]
        ## current: output is CollapsedPrimer3 obj
        print("Generating primary primer list")
        primers_nr = primer_design(pos[1], genomic_template, primary_PCR_amplicon_size)

        ## new output of pseudo_primer_blast (now invalid_primer_pairs) is a set of disallowed primer pairs
        ## previously: output stored in variable 'filtered_primer_list', with structure [pos, <left primer seqs>, <right primer seqs>]
        ## current: output is object of class PrimersOffTargetChecker
        ## NOTE: previous output was VALID primers. current output is object of class PrimersOffTargetChecker
        print("Filtering primary primer list")
        primers_nr.off_target_checker = make_off_target_checker(
            primers_nr, template_DNA_fasta, genomic_DNA_fasta,
            min_acceptable_off_target_primary_amplicon_size)

        ## if no valid primers, increase primary amplicon size until there ARE valid primers or until max amplicon size is reached
        i = 1
        while (
                ## and no uniquely binding primers
                primers_nr.off_target_checker.uniquely_binding_primers() == []
                ## there are no valid primer pairs
                and not primers_nr.has_primer_pairs(
                    check_valid = True, max_amplicon_size = max_primary_amplicon_size,
                    min_acceptable_off_target_amplicon_size = min_acceptable_off_target_primary_amplicon_size,
                    primer3_args = primers_nr.params, full = False)
        ):

            # Break if max amplicon size is reached
            if primary_PCR_amplicon_size[1] > max_primary_amplicon_size:
                break

            primary_PCR_amplicon_size[1] += increment
            print(f"No primary primer found. Increasing allowed primary amplicon size to {primary_PCR_amplicon_size}")

            # Generate new primary primer list with increased allowed amplicon size
            print("Re-generating primary primer list")
            new_primers_nr = primer_design(pos[1], genomic_template, primary_PCR_amplicon_size)

            print("Filtering primary primer list")
            new_primers_nr.off_target_checker = make_off_target_checker(
                new_primers_nr, template_DNA_fasta, genomic_DNA_fasta,
                min_acceptable_off_target_primary_amplicon_size,
                exclude_collapsed_seqs = primers_nr,
                intersect_primers_off_target_checker = primers_nr.off_target_checker,
                prefix = f"{i}|")
            
            primers_nr = new_primers_nr
            i += 1

        print("Re-generating primary primer partners for uniquely-binding primers")
        ## note that this ADDS new primers, doesn't overwrite existing primers in primers_nr
        primers_nr_primary = primer_design_for_uniquely_binding_primers(
            pos[1], primers_nr, genomic_template, primary_PCR_amplicon_size, overwrite = False)

        print("Filtering final primary primer list")
        primers_nr_primary.off_target_checker = make_off_target_checker(
            primers_nr_primary, template_DNA_fasta, genomic_DNA_fasta,
            min_acceptable_off_target_primary_amplicon_size,
            exclude_collapsed_seqs = primers_nr,
            merge_primers_off_target_checker = primers_nr.off_target_checker,
            prefix = "uniq|")

        print("Generating primary primer pairs")
        primer_pairs_primary = primers_nr_primary.primer_pairs(
            check_valid = True, max_amplicon_size = max_primary_amplicon_size,
            min_acceptable_off_target_amplicon_size = min_acceptable_off_target_primary_amplicon_size,
            primer3_args = primers_nr.params, full = False, specific_in_template = 1)
        ## collapse to unique primer sequence combinations and sort by penalty
        primer_pairs_primary = sort_unique_primer_pairs(primer_pairs_primary)

        ### Add primers to df ###
        logging.info(''.join(map(str, ["  # valid primary primer pairs: ", len(primer_pairs_primary)])))
        for i in range(num_primer_pair_per_cleavage_pos):
            if i < len(primer_pairs_primary):
                primer_pair_seq = primer_pairs_primary[i].as_seq_tuple()
                df.loc[index, 'F1'] = primer_pair_seq[0]
                df.loc[index, 'R1'] = primer_pair_seq[1]
            else:
                df.loc[index, 'F1'] = "NA"
                df.loc[index, 'R1'] = "NA"
            index += 1

        ## reset allowed amplicon size
        primary_PCR_amplicon_size[1] = original_value

        #####################
        ##  SECONDARY PCR  ##
        #####################
        
        print("[[ Secondary amplicon ]]")

        # Update secondary PCR primer (second PCR)
        print("Processing amplicon for secondary PCR")
        max_amplicon2_len = int(secondary_PCR_amplicon_size[1]/2)
        local_genomic_template = genomic_template[int(pos[1]-max_amplicon2_len) : int(pos[1] + max_amplicon2_len)]
        local_template_DNA_fasta = "primerg_local_template.fasta"
        SeqIO.write([SeqRecord(Seq(local_genomic_template), id = "local_template", description = "local_template")],
                    local_template_DNA_fasta, "fasta")
        local_pos = int(0.5*len(local_genomic_template))

        print("Generating secondary primer list")
        primers_nr = primer_design(local_pos, local_genomic_template, secondary_PCR_amplicon_size)

        print("Screening secondary primer list against genomic template to get globally specific primers")
        # Filter by genomic template to get globally specific primers
        primers_nr.off_target_checker = make_off_target_checker(
            primers_nr, local_template_DNA_fasta, genomic_DNA_fasta,
            min_acceptable_off_target_secondary_amplicon_size)

        print("Generating globally specific secondary primer pairs")
        primer_pairs_secondary = primers_nr.primer_pairs(
            check_valid = True, max_amplicon_size = secondary_PCR_amplicon_size[1], include_pos = local_pos,
            min_acceptable_off_target_amplicon_size = min_acceptable_off_target_secondary_amplicon_size,
            full = False, specific_in_template = 1)
        ## collapse to unique primer sequence combinations and sort by penalty
        primer_pairs_secondary = sort_unique_primer_pairs(primer_pairs_secondary)

        if len(primer_pairs_secondary) == 0:
            print("Globally specific primers are non-exstient")
            print("Screening secondary primer list against first PCR amplicon to get locally specific primers")
            # Filter by genomic template to get locally specific primers
            primers_nr = primer_design_for_uniquely_binding_primers(
                local_pos, primers_nr, local_genomic_template, secondary_PCR_amplicon_size,
                unique_in_template = True, overwrite = True)

            print("Generating locally specific secondary primer pairs")
            primer_pairs_secondary = primers_nr.primer_pairs(
                check_valid = True, max_amplicon_size = secondary_PCR_amplicon_size[1], include_pos = local_pos,
                min_acceptable_off_target_amplicon_size = min_acceptable_off_target_secondary_amplicon_size,
                full = False, specific_in_template = 1)
            ## collapse to unique primer sequence combinations and sort by penalty
            primer_pairs_secondary = sort_unique_primer_pairs(primer_pairs_secondary)

        for i in range(num_primer_pair_per_cleavage_pos):
            try:
                primer_pair_seq = primer_pairs_secondary[i].as_seq_tuple()
                df.loc[start_index, 'F2'] = primer_pair_seq[0]
                df.loc[start_index, 'R2'] = primer_pair_seq[1]
                amplicon = generate_amplicon(
                    local_pos, df.loc[start_index, 'F2'], df.loc[start_index, 'R2'],
                    local_genomic_template, secondary_PCR_amplicon_size)
                ampliconF = amplicon[len(df.loc[start_index, 'F2']):secondary_PCR_amplicon_size[1]-len(df.loc[start_index, 'F2'])]
                ampliconR = amplicon[-(secondary_PCR_amplicon_size[1]) : -len(df.loc[start_index, 'R2']):]

                countF = num_exact_in_fasta(ampliconF, genomic_DNA_fasta)
                countR = num_exact_in_fasta(ampliconR, genomic_DNA_fasta)

                if countF == 1:
                    df.loc[start_index, 'unique_F2'] = 1
                elif countF > 1:
                    df.loc[start_index, 'unique_F2'] = 0
                else:
                    df.loc[start_index, 'unique_F2'] = "Error: Amplicon F count in genome is negative"

                if countR == 1:
                    df.loc[start_index, 'unique_R2'] = 1
                elif countR > 1:
                    df.loc[start_index, 'unique_R2'] = 0
                else:
                    df.loc[start_index, 'unique_R2'] = "Error: Amplicon R count in genome is negative"

                start_index += 1

            except IndexError: # when no secondary primer is specific to even second PCR amplicon
                df.loc[start_index, 'F2'] = "NA"
                df.loc[start_index, 'R2'] = "NA"
                start_index += 1            

        if print_df:
            print("Dataframe updated:")    
            print(df)
    
    return df
