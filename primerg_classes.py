#########################################
##  CLASSES & MISCELLANEOUS FUNCTIONS  ##
#########################################

import copy
import itertools
import logging
import primer3
import re

from Bio import SeqIO
from Bio import SearchIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

logging.basicConfig(format="%(message)s")

#####################
##  MISCELLANEOUS  ##
##    FUNCTIONS    ##
#####################

def get_count_dict(iterable):
    output = {}
    for e in iterable:
        output[e] = output.get(e, 0) + 1
    return output

def any_in_set(s):
    for e in s:
        return e

###############
##  CLASSES  ##
###############

class Primer():
    ## stats is the dictionary associatd with a given primer in the output of primer3.bindings.design_primers
    ## e.g. the dict that is <primer3.bindings.design_primers output>["PRIMER_LEFT"][0]
    def __init__(self, seqid, stats):
        self.seqid = seqid
        self.stats = stats
        self.seq = str(stats["SEQUENCE"]).upper()
    def get_stat(self, stat_name):
        return self.stats.get(stat_name, None)

## class that tracks all seqids associated with a given sequence
class CollapsedSeq():
    def __init__(self, seq, *seqids):
        self.seq = str(seq).upper()
        self.ids = set(seqids)
    def add_id(self, seqid):
        self.ids.update({seqid})
    def expand_seqs(self):
        ## returns dict of {<seqid>: <seq>} for all associated seqids
        return {seqid: self.seq for seqid in self.ids}
    def random_id(self):
        for seqid in self.ids:
            return seqid

## class that tracks multiple CollapsedSeq objects
class CollapsedSeqs():
    def __init__(self):
        self.collapsed_seqs = {} ## dict of {<seq>: <CollapsedSeq obj>}
    def get_collapsed_seq_by_id(self, seqid):
        for collapsed_seq in self.collapsed_seqs.values():
            if seqid in collapsed_seq.ids:
                return collapsed_seq
        return None
    def get_seq_by_id(self, seqid):
        collapsed_seq = self.get_collapsed_seq_by_id(seqid)
        if collapsed_seq is not None: return collapsed_seq.seq
        return None
    def iter_seqs(self):
        for seq in self.collapsed_seqs:
            yield seq
        return
    def has_seq(self, seq):
        return str(seq).upper() in self.collapsed_seqs
    def get_collapsed_seq_by_seq(self, seq):
        return self.collapsed_seqs.get(str(seq).upper(), None)
    def add_collapsed_seq(self, seq):
        if not self.has_seq(seq):
            self.collapsed_seqs[str(seq).upper()] = CollapsedSeq(seq)
        return
    def add_seq(self, seq, seqid):
        self.add_collapsed_seq(seq)
        collapsed_seq = self.get_collapsed_seq_by_seq(seq)
        collapsed_seq.add_id(seqid)
        return
    ## returns dict of {<seqid>: <seq>} for [1] all sequences if seqs=None or [2] only sequences in seq if seqs!=None
    def expand_seqs(self, *seqs):
        output = {}
        seqs = set(str(seq).upper() for seq in seq)
        for collapsed_seq in self.collapsed_seqs:
            if not seqs or collapsed_seq.seq in seqs:
                output.update(collapsed_seq.expand_seqs())
        return output
    ## writes seqs to fasta file
    ## [1] all sequences if seqs=None or [2] only sequences in seq if seqs!=None
    def write_fasta(self, fout, *seqs, exclude = None, prefix = None):
        if prefix is None: prefix = ''
        seqs = self.expand_seqs(*seqs)
        seq_records = [SeqRecord(Seq(seq), id = seqid, description = seqid) for seqid, seq in seqs.items()]
        SeqIO.write(seq_records, fout, "fasta")
        return
    ## writes unique sequences to fasta file with arbitrary sequence ids
    def write_fasta_uniq(self, fout, prefix = None, exclude = None):
        if prefix is None: prefix = ''
        if exclude is not None:
            exclude_seqs = set(collapsed_seq for collapsed_seq in exclude.iter_seqs())
        else:
            exclude_seqs = set()
        seq_records = [SeqRecord(Seq(seq), id = f"collapsed_{i}", description = f"collapsed_{i}")
                       for i, seq in enumerate(self.iter_seqs())
                       if seq not in exclude_seqs]
        logging.info(''.join(map(str, ["  # exclude: ", len(exclude_seqs), "; # include: ", len(seq_records)])))
        logging.info(''.join(map(str, ["  # unique seqs to write: ", len(seq_records)])))
        SeqIO.write(seq_records, fout, "fasta")
        return

class CollapsedPrimers(CollapsedSeqs):
    def __init__(self):
        super().__init__()
        self.primers = {} ## dict of {<primer id>: <Primer obj>}
        self.primers_left = set()  ## additional variable to track ids of left and right primers separately
        self.primers_right = set() ## additional variable to track ids of left and right primers separately
    def add_primer(self, primer, strand): ## primer is a Primer obj
        self.add_seq(primer.seq, primer.seqid)
        self.primers[primer.seqid] = primer
        if strand in (1, '+', "Plus"):
            self.primers_left.update({primer.seqid})
        elif strand in (-1, '-', "Minus"):
            self.primers_right.update({primer.seqid})
        return
    ## returns set of primer ids for left primers
    ## for primers with multiple ids associated with a single sequence, an arbitrary id will be selected
    def primers_left_unique(self):
        output = []
        for seq, collapsed_seq in self.collapsed_seqs.items():
            # print(seq, collapsed_seq.ids)
            for seqid in collapsed_seq.ids:
                if seqid in self.primers_left:
                    output.append(seqid)
                    break
        return output
    ## returns set of primer ids for right primers
    ## for primers with multiple ids associated with a single sequence, an arbitrary id will be selected
    def primers_right_unique(self):
        output = []
        for seq, collapsed_seq in self.collapsed_seqs.items():
            for seqid in collapsed_seq.ids:
                if seqid in self.primers_right:
                    output.append(seqid)
                    break
        return output
    ## return_type options for plus/minus strand indicators:
    ## [1] 'i' (int) OR 'n' (numeric) --> 1/-1 integers
    ## [2] 'c' (char) --> '+'/'-' characters
    ## [3] 'w' (word) --> "plus"/"minus" strings (lower case)
    ## [4] 'W' (Word) --> "Plus"/"Minus" strings (first letter capitalised)
    def get_strand_by_seqid(self, seqid, return_type = 'c'):
        if return_type == 'i' or return_type == 'n': output = (1, -1)
        elif return_type == 'c': output = ('+', '-')
        elif return_type == 'w': output = ("plus", "minus")
        elif return_type == 'W': output = ("Plus", "Minus")
        else: raise Exception("Invalid input for 'return_type'. Options: 'i', 'n', 'c', 'w', 'W'")
        return (output[0] if seqid in self.primers_left else output[1] if seqid in self.primers_right else None)
    ## this returns set of strands (we're using an iterable to allow for primers that bind to both + and - strands)
    def get_strand_by_seq(self, seq, return_type = 'c'):
        collapsed_seq = self.get_collapsed_seq_by_seq(seq)
        seq_ids = collapsed_seq.ids
        return set(self.get_strand_by_seqid(primer_id, return_type = return_type) for primer_id in seq_ids)

## assumes primers are designed on same molecule/chrom/template
class PrimerPairPrimer3():
    def __init__(self, primer_A, primer_B):
        self.primer_A = primer_A
        self.primer_B = primer_B
        ## cache
        self._valid = None
        self._primer_plus = None
        self._primer_minus = None
    ## these two methods return lists
    def primers_plus(self):
        return [primer for primer in [self.primer_A, self.primer_B] if primer.get_stat("STRAND") in (1, '+', "Plus")]
    def primers_minus(self):
        return [primer for primer in [self.primer_A, self.primer_B] if primer.get_stat("STRAND") in (-1, '-', "Minus")]
    ## assumes only one primer on plus strand. if multiple, returns first primer in self.primers_plus output
    def primer_plus(self):
        if self._primer_plus is not None: return self._primer_plus
        primers = self.primers_plus()
        output = None if not primers else primers[0]
        self._primer_plus = output
        return output
    ## assumes only one primer on minus strand. if multiple, returns first primer in self.primers_plus output
    def primer_minus(self):
        if self._primer_minus is not None: return self._primer_minus
        primers = self.primers_minus()
        output = None if not primers else primers[0]
        self._primer_minus = output
        return output
    ## self.valid() returns True if and only if:
    ## [1] primers are on opposite strands
    ## [2] extension is towards each other
    ## [3] amplicon size is within desired range
    ## [4] optional: position to include is between primers (but not in primers)
    def valid(self, min_amplicon_size = 0, max_amplicon_size = float("Inf"), include_pos = None):
        if self._valid is not None: return self._valid
        ## check there is one on each strand
        strand_plus = self.primers_plus()
        strand_minus = self.primers_minus()
        if not len(strand_plus) == len(strand_minus) == 1:
            self._valid = False
            return False
        ## check if primers are an appropriate distance apart
        primer_plus_start = self.primer_plus().get_stat("START")
        primer_minus_end = self.primer_minus().get_stat("END")
        primer_dist = primer_minus_end - primer_plus_start
        if not min_amplicon_size <= primer_dist <= max_amplicon_size:
            self._valid = False
            return False
        ## check if desired position is within amplicon
        if include_pos is not None:
            output = primer_plus_start <= include_pos < primer_minus_end
            self._valid = output
            return output
        ## return True if all valid checks passed
        self._valid_cache = True
        return True
    ## returns (None, None) if not a valid primer pair, else returns (<str of plus primer seq>, <str of minus primer seq>)
    def as_seq_tuple(self):
        return (self.primer_plus().seq, self.primer_minus().seq)
    ## returns total penalty of primers
    def penalty(self):
        return (self.primer_A.get_stat("PENALTY") + self.primer_B.get_stat("PENALTY"))

class CollapsedPrimer3(CollapsedPrimers):
    def __init__(self, design_primers_output):
        super().__init__()
        self.params = None
        self.target = None
        ## store raw output
        self.raw = design_primers_output
        # print("primer left:", self.raw["PRIMER_LEFT"])
        # print("primer right:", self.raw["PRIMER_RIGHT"])
        self.off_target_checker = None
        self._parse()
    def _parse(self):
        ## custom stats START and END are 0-indexed, end-exclusive
        for i, primer in enumerate(self.raw["PRIMER_LEFT"]):
            start, length = primer["COORDS"]
            new_primer_obj = Primer(f"PRIMER_LEFT_{i}",
                                    {**primer, **{"START": start, "END": start+length, "STRAND": '+'}})
            self.add_primer(new_primer_obj, strand = 1)
        for i, primer in enumerate(self.raw["PRIMER_RIGHT"]):
            end, length = primer["COORDS"]
            new_primer_obj = Primer(f"PRIMER_RIGHT_{i}",
                                    {**primer, **{"START": end-length-1, "END": end, "STRAND": '-'}})
            self.add_primer(new_primer_obj, strand = -1)
        return
    ## yields PrimerPairPrimer3 objects
    def _primer_pairs(self, check_valid = True,
                      min_amplicon_size = 0, max_amplicon_size = float("Inf"), include_pos = None,
                      primers_off_target_checker = None,
                      min_acceptable_off_target_amplicon_size = 0, primer3_args = None,
                      full = False, specific_in_template = 1):
        if include_pos is None: include_pos = self.target
        logging.debug(' '.join(map(str, [min_amplicon_size, max_amplicon_size, include_pos])))
        ## get invalid primer pairs (sequences) per primers_off_target_checker
        if primers_off_target_checker is None:
            primers_off_target_checker = self.off_target_checker
        if check_valid and primers_off_target_checker is not None:
            invalid_pairs_off_target = set(primers_off_target_checker.invalid_primer_pairs(
                min_acceptable_off_target_amplicon_size = min_acceptable_off_target_amplicon_size,
                seq = True, specific_in_template = specific_in_template
            ))
            logging.info(''.join(map(str, ["  # invalid pairs by off-target: ", len(invalid_pairs_off_target)])))
            logging.debug(''.join(map(str, ["    # primer left in invalid: ",
                                            len(set(x[0] for x in invalid_pairs_off_target))])))
            logging.debug(''.join(map(str, ["    # primer right in invalid: ",
                                            len(set(x[1] for x in invalid_pairs_off_target))])))
            logging.debug(''.join(map(str, ["    # primer total in invalid: ",
                                            len(set(itertools.chain(*invalid_pairs_off_target)))])))
        else:
            invalid_pairs_off_target = set()
        logging.debug(''.join(map(str, ["  # invalid primer pairs: ", len(invalid_pairs_off_target)])))
        is_invalid_off_target_pair = lambda primer_pair: primer_pair.as_seq_tuple() in invalid_pairs_off_target
        ## function to check if primer pair is valid per primer3 QC
        if primer3_args is None:
            is_valid_primer3_pair = lambda primer_pair: True
        else:
            is_valid_primer3_pair = lambda primer_pair: (
                primer3.bindings.design_primers(
                    primer3_args[0],
                    {**primer3_args[1], **{"SEQUENCE_PRIMER": primer_pair.primer_plus().seq,
                                           "SEQUENCE_PRIMER_REVCOMP": primer_pair.primer_minus().seq}}
                )["PRIMER_PAIR_NUM_RETURNED"] > 0)
        ## get valid primer pairs
        # i = 0
        logging.debug(''.join(map(str, ["  # primer left: ", len(self.primers_left),
                                        "; # primer right: ", len(self.primers_right)])))
        logging.debug(''.join(map(str, ["  # unique primer seqs: ", len(self.collapsed_seqs)])))
        ## if not full, keep only first instance of primer sequence encountered
        primers_left = self.primers_left if full else (self.primers_left_unique())
        primers_right = self.primers_right if full else (self.primers_right_unique())
        logging.info(''.join(map(str, ["  # primer left for screening: ", len(primers_left),
                                       "; # primer right for screening: ", len(primers_right)])))
        ## iterate through combinations
        track_primer_pair_failure = {"combos": 0, "invalid_off-target": 0,
                                     "invalid_on-target": 0, "invalid_primer3": 0}
        for id_left in primers_left:
            primer_left = self.primers[id_left]
            for id_right in primers_right:
                # i += 1
                # if i % 500 == 0: print(i)
                primer_right = self.primers[id_right]
                primer_pair = PrimerPairPrimer3(primer_left, primer_right)
                primer_pair_seq = primer_pair.as_seq_tuple()
                track_primer_pair_failure["combos"] += 1
                ## if invalid pair, continue to next right primer
                if check_valid:
                    if is_invalid_off_target_pair(primer_pair):
                        track_primer_pair_failure["invalid_off-target"] += 1
                        continue
                    if (include_pos is not None
                        and not primer_pair.valid(min_amplicon_size = min_amplicon_size,
                                                  max_amplicon_size = max_amplicon_size,
                                                  include_pos = include_pos)):
                        track_primer_pair_failure["invalid_on-target"] += 1
                        continue
                    if not is_valid_primer3_pair(primer_pair):
                        track_primer_pair_failure["invalid_primer3"] += 1
                        continue
                yield (primer_pair)
        # print(primer_pair.as_seq_tuple(), list(invalid_pairs_off_target)[:5])
        logging.info(''.join(map(str, ["  # combos: ", track_primer_pair_failure["combos"]])))
        logging.info(''.join(map(str, ["    # invalid after off-target: ", track_primer_pair_failure["invalid_off-target"]])))
        logging.info(''.join(map(str, ["    # invalid on-target: ", track_primer_pair_failure["invalid_on-target"]])))
        logging.info(''.join(map(str, ["    # invalid primer3: ", track_primer_pair_failure["invalid_primer3"]])))
        return
    ## returns list of PrimerPairPrimer3 objects
    ## args and kwargs for self._primer_pairs
    def primer_pairs(self, *args, **kwargs):
        output = {}
        for primer_pair in self._primer_pairs(*args, **kwargs):
            primer_pair_seq = primer_pair.as_seq_tuple()
            output[primer_pair_seq] = output.get(primer_pair_seq, []) + [primer_pair]
        ## flatten output and return
        return list(itertools.chain(*output.values()))
    ## args and kwargs for self._primer_pairs
    def has_primer_pairs(self, *args, **kwargs):
        for primer_pair in self._primer_pairs(*args, **kwargs):
            return True ## return True if self._primer_pairs yields ANYTHING
        return False

class PrimersOffTargetChecker():
    def __init__(self, *primers_fastas):
        self.primers = {} ## nested dict indexed by subject_id, strand, pos (tuple of (min, max)), and finally a set of query_id
        self.primers_map = {record.id: str(record.seq)
                            for primers_fasta in primers_fastas
                            for record in SeqIO.parse(primers_fasta, "fasta")} ## dict of {<primer id>: <primer seq>}
        self.forgive_coords = {}
        self.template_coords = {}
    def seqids(self):
        return set(self.primers_map.keys())
    def seqs(self):
        return set(self.primers_map.values())
    def _merge_primers(self, other, prefix = ''):
        for subject_id, subject_dat in other.primers.items():
            updated_subject_dat = self.primers.get(subject_id, {})
            for strand, strand_dat in subject_dat.items():
                updated_strand_dat = updated_subject_dat.get(strand, {})
                for pos, query_ids in strand_dat.items():
                    updated_strand_dat[pos] = updated_strand_dat.get(pos, set()).union((prefix+query_id) for query_id in query_ids)
                updated_subject_dat[strand] = updated_strand_dat
            self.primers[subject_id] = updated_subject_dat
        return
    def _merge_primers_map(self, other, prefix = ''):
        self.primers_map = {**self.primers_map,
                            **{prefix+primer_id: primer_seq for primer_id, primer_seq in other.primers_map.items()}}
        return
    ## this isn't quite merging. overlapping ranges will NOT be merged. (however, duplicated ranges will be discarded)
    def _merge_forgive_coords(self, other):
        for subject_id, coords_iter in other.forgive_coords.items():
            self.forgive_coords[subject_id] = self.forgive_coords.get(subject_id, set()).union(coords_iter)
        return
    ## this isn't quite merging. overlapping ranges will NOT be merged. (however, duplicated ranges will be discarded)
    def _merge_template_coords(self, other):
        for subject_id, coords_iter in other.template_coords.items():
            self.template_coords[subject_id] = self.template_coords.get(subject_id, set()).union(coords_iter)
        return
    def merge(self, other, prefix = ''):
        ## merge PrimersOffTargetChecker.primers
        self._merge_primers(other, prefix = prefix)
        ## merge primers_map
        self._merge_primers_map(other, prefix = prefix)
        ## merge forgive_coords
        self._merge_forgive_coords(other)
        ## merge template_coords
        self._merge_template_coords(other)
        return
    def _subset_primers(self, seqids):
        seqids = set(seqids)
        output = {}
        for subject_id, subject_dat in self.primers.items():
            new_subject_dat = {}
            for strand, strand_dat in subject_dat.items():
                new_strand_dat = {}
                for pos, query_ids in strand_dat.items():
                    intersect = seqids & query_ids
                    if intersect:
                        new_strand_dat[pos] = intersect
                if new_strand_dat:
                    new_subject_dat[strand] = new_strand_dat
            if new_subject_dat:
                output[subject_id] = new_subject_dat
        return output
    def _subset_primers_map(self, seqids):
        seqids = set(seqids)
        return {seqid: seq for seqid, seq in self.primers_map.items() if seqid in seqids}
    ## this returns a new PrimersOffTargetCheckef object
    def subset(self, seqs = [], seqids = []):
        ## get all seqids to keep
        seqids = set(seqid for seqid in seqids if seqid in self.primers_map)
        for seq in set(seqs):
            seqids.update(known_seqid for known_seqid, known_seq in self.primers_map.items() if known_seq == seq)
        ## create new object
        new = PrimersOffTargetChecker()
        ## subset primers
        new.primers = self._subset_primers(seqids)
        ## subset primers_map
        new.primers_map = self._subset_primers_map(seqids)
        ## deep copy forgive_coords & template_coords
        new.forgive_coords = copy.deepcopy(self.forgive_coords)
        new.template_coords = copy.deepcopy(self.template_coords)
        return new
    def id_to_seq(self, seqid):
        return self.primers_map.get(seqid, None)
    def add_primer(self, query_id, subject_id, strand, start, end):
        pos = (min(start, end)-1, max(start, end)) ## make it 0-indexed AND end-exclusive
        primers_subject = self.primers.get(subject_id, {})
        primers_strand = primers_subject.get(strand, {})
        primers_pos = primers_strand.get(pos, set())
        primers_pos.update({query_id})
        primers_strand[pos] = primers_pos
        primers_subject[strand] = primers_strand
        self.primers[subject_id] = primers_subject
        return
    ## returns:
    ## [1] if subject_id!=None: dict of {(<start>,<end>): <set of primer ids>}
    ## [2] if subject_id==None: dict of {<subjec id>:{(<start>,<end>): <set of primer ids>}}
    def primers_pos_plus(self, subject_id = None):
        if subject_id:
            primer_strands = self.primers.get(subject_id, {})
            return primer_strands.get(1, primer_strands.get('+', primer_strands.get("Plus", {})))
        else:
            output = {}
            for subject_id, primer_strands in self.primers.values():
                output = {**output,
                          **{subject_id:
                             primer_strands.get(1, primer_strands.get('+', primer_strands.get("Plus", {})))}}
            return output
    def primers_pos_minus(self, subject_id = None):
        if subject_id:
            primer_strands = self.primers.get(subject_id, {})
            return primer_strands.get(-1, primer_strands.get('-', primer_strands.get("Minus", {})))
        else:
            output = {}
            for primer_strands in self.primers.values():
                output = {**output,
                          **primer_strands.get(-1, primer_strands.get('-', primer_strands.get("Minus", {})))}
            return output
    ## forgive_coords is a dict of {<subject_id>: [(<start>, <end>), (<start>, <end>)]}
    ## returns a function that returns true if all provided numeric values are between ONE range for a given subject_id
    ## if pos are spread over multiple ranges or subject_ids, returns False
    def _make_within_forgive_coords(self, forgive_coords = None):
        if forgive_coords is None:
            forgive_coords = self.forgive_coords
        def helper(subject_id, *pos):
            coords = forgive_coords.get(subject_id)
            if coords:
                for start, end in coords:
                    if all(map(lambda p: p>=start and p<=end, pos)):
                        return True
            return False
        return helper
    ## forgive_coords is a dict of {<subject_id>: [(<start>, <end>), (<start>, <end>)]}
    ## forgive_coords tells to function to ignore any off-target amplicons within the region (both primers of a pair must be within)
    ## returns:
    ## [A] IF seq=False: set of tuple of invalid +- primer combinations {(<plus_primer_id>, <minus_primer_id>)}
    ## [B] IF seq=True: set of tuple of invalid +- primer combinations {(<plus_primer_seq>, <minus_primer_seq>)}
    ## TODO: fix!! this function is allowing primers that aren't specific in template to pass!!
    ## if specific_in_template=0: ignore specificity in template
    ## elif specific_in_template=1: at least one primer of a pair must be specific in template
    ## elif specific_in_template>1: >1 primers of a pair (in practice this is just 2 primers) must be specific in template
    def invalid_primer_pairs(self, min_acceptable_off_target_amplicon_size = 0,
                             forgive_coords = None, seq = False, specific_in_template = 1):
        ## generate filtering functions
        within_forgive_coords = self._make_within_forgive_coords(forgive_coords = forgive_coords)
        if specific_in_template > 0:
            primer_id_specific_in_template = set(self.uniquely_binding_primers_in_template(seq = False))
            is_specific_in_template = lambda primer_id: primer_id in primer_id_specific_in_template
        else:
            is_specific_in_template = lambda primer_id: True
        ## generate output formatting function
        if seq: format_primer_pair = lambda plus_id, minus_id: (self.id_to_seq(plus_id), self.id_to_seq(minus_id))
        else: format_primer_pair = lambda plus_id, minus_id: (plus_id, minus_id)
        ## tracker & updater for invalid combos
        primers_invalid = set()
        def update_invalid(plus_ids, minus_ids):
            for plus_id in plus_ids:
                for minus_id in minus_ids:
                    primers_invalid.update({format_primer_pair(plus_id, minus_id)})
            return
        ## iterate through all subjects
        for subject_id, primer_strands in self.primers.items():
            if len(primer_strands) < 2: continue
            primers_plus = self.primers_pos_plus(subject_id = subject_id)
            primers_minus = self.primers_pos_minus(subject_id = subject_id)
            ## iterate through all combinations of + and - primers for the subject_id in order of position
            for pos_plus in sorted(primers_plus):
                # ## DISALLOW: if not within forgive_coords
                # if not within_forgive_coords(subject_id, *pos_plus):
                #     ## add all possible combinations of current pos_plus primer and primers
                #     ## for all pos_minus to primers_invalid
                #     for pos_minus in primers_minus:
                #         update_invalid(primers_plus[pos_plus], primers_minus[pos_minus])
                #     continue
                ## DISALLOW: if specific_in_template=True and pos_plus primer is not specific in template
                plus_ids_not_specific_in_template = set(
                    plus_id for plus_id in primers_plus[pos_plus]
                    if not is_specific_in_template(plus_id)
                )
                ## screen plus-minus combos
                for pos_minus in sorted(primers_minus):
                    ## ALLOW: if + and - primers occupy same position
                    ## (or if + pos is greater than - pos), move on to next - primer
                    if pos_plus[0] >= pos_minus[0]:
                        continue
                    ## ALLOW: if + and - primers are >= min_acceptable_off_target_amplicon_size apart
                    elif pos_minus[0] - pos_plus[-1] >= min_acceptable_off_target_amplicon_size:
                        continue
                    ## CONDITIONAL ALLOW: if within forgive_coords
                    elif within_forgive_coords(subject_id, *pos_plus, *pos_minus):
                        minus_ids_not_specific_in_template = set(
                            minus_id for minus_id in primers_minus[pos_minus]
                            if not is_specific_in_template(minus_id)
                        )
                        ## DISALLOW: if specific_in_template> 1
                        ##           and >=1 pos_plus primer and/or >= 1 pos_minus primer is not specific in template
                        if specific_in_template > 1:
                            if plus_ids_not_specific_in_template: ## disallow combos with non-specific + primer
                                update_invalid(plus_ids_not_specific_in_template, primers_minus[pos_minus])
                            if minus_ids_not_specific_in_template: ## disallow combos with non-specific - primer
                                update_invalid(primers_plus[pos_plus], minus_ids_not_specific_in_template)
                        ## DISALLOW: if specific_in_template=1
                        ##           and >=1 pos_plus primer AND >= 1 pos_minus primer are not specific in template
                        elif (specific_in_template == 1
                              and plus_ids_not_specific_in_template
                              and minus_ids_not_specific_in_template):
                            ## disallow combos where BOTH primers are non-specific in template
                            update_invalid(plus_ids_not_specific_in_template, minus_ids_not_specific_in_template)
                        ## ALLOW otherwise
                        continue
                    ## DISALLOW otherwise
                    update_invalid(pos_plus, pos_minus)
        return primers_invalid
    ## get number of unique binding positions for all primers
    def primer_binding_counts(self):
        output = get_count_dict(
            [primer_id
             for subject_dat in self.primers.values()
             for strand_dat in subject_dat.values()
             for pos_dat in strand_dat.values()
             for primer_id in pos_dat]
        )
        return output
    ## returns dict of {(<template chrom>, (<template start>, <template end>)): {<primer_id>: <binding count in template>}
    ## where template chrom, start, end are the positions of the template in the genomic background
    def primer_binding_counts_in_template(self):
        output = {}
        for subject_id, coords_multi in self.template_coords.items():
            for coords in coords_multi:
                template_raw_primer_instances = []
                for strand, strand_dat in self.primers[subject_id].items():
                    for pos, primer_ids in strand_dat.items():
                        if pos[0] >= coords[0] and pos[1] <= coords[1]:
                            template_raw_primer_instances.extend(list(primer_ids))
                output[(subject_id, coords)] = get_count_dict(template_raw_primer_instances)
        return output
    ## returns dict of {<primer_id>: <maximum binding count in any template>}
    def max_primer_binding_counts_per_template(self):
        primer_binding_counts_in_template = self.primer_binding_counts_in_template()
        primer_ids = set(itertools.chain(*[x.keys() for x in primer_binding_counts_in_template.values()]))
        output = {
            primer_id: max(template_dat.get(primer_id, 0)
                           for template_dat in primer_binding_counts_in_template.values())
            for primer_id in primer_ids
        }
        # print("max primer binding counts per template:", output)
        return output
    ## returns list of primer_ids of primers that only have one binding site
    def uniquely_binding_primers(self, seq = False):
        primer_ids = [primer_id
                      for primer_id, binding_count
                      in self.primer_binding_counts().items()
                      if binding_count == 1]
        if seq: return [self.id_to_seq(primer_id) for primer_id in primer_ids]
        else: return primer_ids
    ## returns list of primer_ids of primers that have no more than one binding site within any template region
    def uniquely_binding_primers_in_template(self, seq = False):
        primer_ids = [primer_id
                      for primer_id, binding_count
                      in self.max_primer_binding_counts_per_template().items()
                      if binding_count == 1]
        # print(primer_ids)
        if seq: return [self.id_to_seq(primer_id) for primer_id in primer_ids]
        else: return primer_ids


########################
##  BLAST FORMATTING  ##
########################

CHAR_UNALIGNED = ' '
CHAR_DEL = 'd'
CHAR_INS = 'i'
CHAR_MISMATCH = 'm'
CHAR_MATCH = '.'
CHAR_GAP = 'g' ## used in MINORg off-target regex
BTOP_GAP = '-' ## used in raw btop pattern (generated by BLAST)

class BlastHSP:
    """
    Class that binds HSP object with QueryResult and Hit as query and subject respectively.
    """
    def __init__(self, hsp, query, subject):
        self.hsp = hsp
        self.query = query
        self.subject = subject
        self.qlen = None if "seq_len" not in dir(self.query) else self.query.seq_len
        self.sbtop = None
        self.tbtop = None
        self.tbtop_rvs = None
        self.parse_btop()
    
    @property
    def btop(self):
        return self.hsp.btop
    
    def parse_btop(self):
        if "btop" in dir(self.hsp):
            self.sbtop = self.expand_btop_str(include_unaligned = True)
            self.tbtop = self.expand_btop_tuple(include_unaligned = True)
            self.tbtop_rvs = self.expand_btop_tuple(include_unaligned = True, rvs_index = True)
    
    def expand_btop_str(self, include_unaligned = False):
        """
        Expands btop to a string of characters of length equal to alignment, where:
            '.' is a match,
            'm' is a mismatch,
            'i' is an insertion in the query,
            'd' is a deletion in the query.
        
        Arguments:
            include_unaligned (bool): add space character for each unaligned position at 3' and 5' ends
        
        Returns
        -------
        str
            Expanded btop pattern
        """
        output = ''
        btop = self.btop
        ## processed characters are removed from btop until none are left
        while btop:
            num = re.search('^\d+', btop)
            if num: ## if match
                num = num.group(0)
                output += CHAR_MATCH*int(num)
                btop = btop[len(num):]
            else: ## if gap or mismatch
                unit = btop[:2]
                if unit[0] == BTOP_GAP: output += CHAR_DEL ## deletion in query
                elif unit[1] == BTOP_GAP: output += CHAR_INS ## insertion in query
                else: output += CHAR_MISMATCH ## mismatch
                btop = btop[2:]
        if include_unaligned:
            output = (CHAR_UNALIGNED*self.hsp.query_start) + \
                     output + \
                     (CHAR_UNALIGNED*(self.qlen-self.hsp.query_end))
        return output
    
    def expand_btop_tuple(self, rvs_index = False, include_unaligned = False):
        """
        Calls str_expand_btop and groups expanded btop to a tuple of characters
        of length equal to aligned query, where:
            '.' is a match,
            'm' is a mismatch,
            'i' is an insertion in the query,
            'd' is a deletion in the query (deletion between N and N+1 will be grouped with position N+1 if rvs_index=False else N).

        Arguments:
            rvs_index (bool): default=False. A deletion between positions N and N+1 will be grouped with:
                position N+1 if rvs_index=False;
                position N, if rvs_index=True.

        Returns
        -------
        tuple of str
            Expanded btop pattern
        """
        str_expanded = self.expand_btop_str()
        output = []
        non_del = {CHAR_MATCH, CHAR_MISMATCH, CHAR_INS}
        ## processed characters are removed from str_expanded until none are left
        while str_expanded:
            c = str_expanded[0]
            if c in non_del: ## if not deletion
                output.append(c)
                str_expanded = str_expanded[1:]
            elif rvs_index: ## if deletion in reversed index
                output[-1] = output[-1] + c
                str_expanded = str_expanded[1:]
            else: ## if deletion in index counting up
                output.append(str_expanded[:2])
                str_expanded = str_expanded[2:]
        if include_unaligned:
            output = ([CHAR_UNALIGNED]*self.hsp.query_start) + output + \
                     ([CHAR_UNALIGNED]*(self.qlen-self.hsp.query_end))
        return tuple(output)
    
    def splice_expanded_btop(self, start = None, end = None, fmt = tuple,
                             rvs_index = False, include_unaligned = True):
        """
        Returns expanded btop spliced to specified range and in specified format
        
        Arguments:
            start (int): optional, start position (default=0)
            end (int): optional, end position (non-inclusive) (default=<query length (self.qlen))
            fmt (type): str or tuple; return expanded btop type
            rvs_index (bool): used when fmt=tuple;
                deletion between N and N+1 will be grouped with position N+1 if rvs_index=False else N
            include unaligned (bool): add space character for each unaligned position at 3' and 5' ends
        
        Returns
        -------
        str or tuple
            Expanded btop spliced to specified range
        """
        if start is None: start = 0
        if end is None: end = self.qlen
        ## originally: incr (int): default=1; splice increment (use -1 when start < end; use 1 when start > end)
        incr = 1 if (end >= start) else -1
        btop = self.tbtop_rvs if rvs_index else self.tbtop
        btop = btop[start:end:incr]
        if not include_unaligned:
            btop = type(btop)(x for x in btop if x != CHAR_UNALIGNED)
        if fmt == str: return ''.join(btop)
        else: return btop
        
    def _num_char_in_range(self, *char, start = None, end = None, rvs_index = False):
        """
        Returns number of times a given set of characters appears within a given range
        
        Arguments:
            *char (str): characters to search for
            start (int): optional, start position (default=0)
            end (int): optional, end position (non-inclusive) (default=<query length (self.qlen))
            rvs_index (bool): used when fmt=tuple;
                deletion between N and N+1 will be grouped with position N+1 if rvs_index=False else N
        
        Returns
        -------
        int
            Number of times characters appear within range
        """
        sbtop = self.splice_expanded_btop(start, end, fmt = str, rvs_index = rvs_index,
                                          include_unaligned = True)
        output = 0
        for c in char:
            output += sbtop.count(c)
        return output
    
    def num_unaligned(self, start = None, end = None):
        """
        Returns number of unaligned positions in query
        
        Arguments:
            start (int): optional, start position (default=0)
            end (int): optional, end position (non-inclusive) (default=<query length (self.qlen))
        
        Returns
        -------
        int
            Number of unaligned positions within range
        """
        return self._num_char_in_range(CHAR_UNALIGNED, start = start, end = end)
    def num_insertion(self, start = None, end = None):
        """
        Returns number of insertions relative to subject
        (that is, one or more bases that is/are present in the subject but not in the query).
        
        Arguments:
            start (int): optional, start position (default=0)
            end (int): optional, end position (non-inclusive) (default=<query length (self.qlen))
        
        Returns
        -------
        int
            Number of insertions within range
        """
        return self._num_char_in_range(CHAR_INS, start = start, end = end)
    
    def num_deletion(self, start = None, end = None, rvs_index = False):
        """
        Returns number of deletions relative to subject
        (that is, one or more bases that is/are present in the query but not in the subject).
        
        Arguments:
            start (int): optional, start position (default=0)
            end (int): optional, end position (non-inclusive) (default=<query length (self.qlen))
            rvs_index (bool): used when fmt=tuple;
                deletion between N and N+1 will be grouped with position N+1 if rvs_index=False else N
        
        Returns
        -------
        int
            Number of deletions within range
        """
        return self._num_char_in_range(CHAR_DEL, start = start, end = end, rvs_index = rvs_index)
    
    def num_mismatch(self, start = None, end = None):
        """
        Returns number of mismatches.
        
        Arguments:
            start (int): optional, start position (default=0)
            end (int): optional, end position (non-inclusive) (default=<query length (self.qlen))
        
        Returns
        -------
        int
            Number of mismatches within range
        """
        return self._num_char_in_range(CHAR_MISMATCH, start = start, end = end)
    
    def num_gap(self, start = None, end = None, rvs_index = False):
        """
        Returns number of gaps.
        
        Arguments:
            start (int): optional, start position (default=0)
            end (int): optional, end position (non-inclusive) (default=<query length (self.qlen))
            rvs_index (bool): used when fmt=tuple;
                deletion between N and N+1 will be grouped with position N+1 if rvs_index=False else N
        
        Returns
        -------
        int
            Number of gaps within range
        """
        return self._num_char_in_range(CHAR_DEL, CHAR_INS, start = start, end = end, rvs_index = rvs_index)

class BlastResult:
    """
    Generator that reads blast-tab format and yields BlastHSP,
    which stores a HSP object with its associated QueryResult and Hit objects 
    as attributes query and subject respectively.
    """
    def __init__(self, filename, fmt, **kwargs):
        """
        Create a BlastResult object.
        
        Arguments:
            filename (str): path to file
            fmt (str): file format (e.g. 'blast-tab', 'blast-xml')
            **kwargs: additional arguments for SearchIO (notably, fields)
        """
        self.filename = filename
        self.fmt = fmt
        self.kwargs = kwargs
    def __iter__(self):
        for query_result in SearchIO.parse(self.filename, self.fmt, **self.kwargs):
            for hit in query_result:
                for hsp in hit:
                    yield BlastHSP(hsp, query_result, hit)
