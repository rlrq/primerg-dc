import os
import shutil
import tempfile
from Bio import Seq
from pyfaidx import Fasta, Faidx, MutableFastaRecord, FastaRecord, OrderedDict

class rFaidx(Faidx):
    """
    Faidx object that writes index into a temporary directory if file location is not writable
    """
    @property
    def indexname(self):
        return self._indexname
    @indexname.setter
    def indexname(self, value):
        if os.access(os.path.dirname(value), os.W_OK):
            self._indexname = value
        else:
            tmpdir = tempfile.mkdtemp()
            self._tmpdirs.append(tmpdir)
            self._indexname = os.path.join(tmpdir, os.path.basename(value))
        return
    ## remove any tmp files
    def __del__(self):
        for tmpdir in self._tmpdirs:
            shutil.rmtree(tmpdir)
        return
    def __init__(self, *args, **kwargs):
        self._indexname = None
        self._tmpdirs = []
        super().__init__(*args, **kwargs)
        return

## basically pyfaidx.Fasta class, but returns Bio.Seq.Seq item when sliced and filename for __repr__
class IndexedFasta(Fasta):
    """
    Basically a pyfaidx.Fasta class that:
        1. returns Bio.Seq.Seq item when sliced
        2. returns IndexedFasta(<filename>) for __repr__
        3. is capable of writing index into a temporary directory if file location is not writable
    """
    def __init__(self,
                 fasta,
                 default_seq=None,
                 key_function=lambda x: x,
                 as_raw=False,
                 strict_bounds=False,
                 read_ahead=None,
                 mutable=False,
                 split_char=None,
                 filt_function=lambda x: True,
                 one_based_attributes=True,
                 read_long_names=False,
                 duplicate_action="stop",
                 sequence_always_upper=False,
                 rebuild=True,
                 build_index=True):
        if isinstance(fasta, self.__class__):
            self.__dict__ = copy.copy(fasta.__dict__) ## no need for deep copy
        else:
            if isinstance(fasta, Fasta):
                filename = fasta.filename
            else:
                filename = fasta
            ## everything below is a reproduction of pyfaidx.Fasta's __init__ method aside from
            ## the substitution of 'Faidx' class with rFaidx' class
            self.filename = filename
            self.mutable = mutable
            self.faidx = rFaidx( ## rFaidx enables indexing of files in non-writable locations
                filename,
                key_function=key_function,
                as_raw=as_raw,
                default_seq=default_seq,
                strict_bounds=strict_bounds,
                read_ahead=read_ahead,
                mutable=mutable,
                split_char=split_char,
                filt_function=filt_function,
                one_based_attributes=one_based_attributes,
                read_long_names=read_long_names,
                duplicate_action=duplicate_action,
                sequence_always_upper=sequence_always_upper,
                rebuild=rebuild,
                build_index=build_index)        
            _record_constructor = MutableFastaRecord if self.mutable else FastaRecord
            self.records = OrderedDict([(rname, _record_constructor(rname, self))
                                        for rname in self.faidx.index.keys()])
        return
    def __repr__(self):
        return f"IndexedFasta({self.filename})"
    def get_seq(self, *args, **kwargs):
        pyfaidx_seq = super().get_seq(*args, **kwargs)
        return Seq.Seq(pyfaidx_seq.seq)
    def get_spliced_seq(self, *args, **kwargs):
        pyfaidx_seq = super().get_spliced_seq(*args, **kwargs)
        return Seq.Seq(pyfaidx_seq.seq)


# class IndexedFasta:
#     '''
#     chunksize must be larger than the maximum length of a single line
#     '''
#     def __init__(self, filename, indices = None, chunksize = 10000, **kwargs_for_read):
#         self.filename = filename
#         self.chunksize = chunksize
#         self.indices = indices
#         self.kwargs_for_read = kwargs_for_read
#         if indices is None:
#             self.index()
#     def _chunk_bin(self, pos):
#         return pos // self.chunksize
#     def _chunk_offset(self, pos):
#         return pos % self.chunksize
#     def index(self, chunksize = None):
#         if chunksize is not None:
#             self.chunksize = chunksize
#         indices = {}
#         seqid = seqpos = None
#         with open(self.filename, 'r', **self.kwargs_for_read) as f:
#             while True:
#                 prev_pos = f.tell()
#                 line = f.readline()
#                 if not line: break ## EOF
#                 line = line.rstrip()
#                 if not line: continue ## empty line
#                 if line[0] == '>':
#                     print("new chrom!")
#                     seqid = line[1:]
#                     indices[seqid] = [f.tell()]
#                     seqpos = 0
#                     continue
#                 new_seqpos = seqpos + len(line)
#                 ## if the position to index is within this line
#                 if self._chunk_bin(seqpos) < self._chunk_bin(new_seqpos):
#                     ## get sequence in line up until the position to index
#                     target_offset = self.chunksize - self._chunk_offset(seqpos)
#                     target_seq = line[:target_offset]
#                     ## return to start of line and increment position in line by
#                     ##   tmp_offset number of bytes until position to index is found
#                     ##   - presumably, char size in other encodings will be multiple of char size in ascii
#                     f.seek(prev_pos)
#                     while True:
#                         incr_seq = f.read(target_offset)
#                         ## if sequence until this position == target_seq, store position and break
#                         if incr_seq == target_seq:
#                             indices[seqid].append(f.tell())
#                             break
#                         else:
#                             incr_seq += f.read(target_offset)
#                 seqpos = new_seqpos
#         self.indices = indices
#     def get(self, seqid, start, end):
#         '''
#         0-indexed, start-inclusive end-exclusive
#         '''
#         ## check validity of inputs
#         if seqid not in self.indices:
#             raise Exception(f"'{seqid}' does not exist.")
#         is_valid_pos = lambda i: i >= 0 and int(i) == i
#         invalid = [i for i in (start, end) if not is_valid_pos(i)]
#         if invalid: raise Exception(f"Invalid indices: {','.join(map(str, invalid))}")
#         if start == end: raise Exception("Start and end coordinates cannot be identical.")
#         if start > end: raise Exception("Start coordinate must be smaller than end coordinate.")
#         ## get bin for each coordinate
#         start_bin = self._chunk_bin(start)
#         start_offset = self._chunk_offset(start)
#         end_offset = end - start + start_offset
#         print("start_bin:", start_bin, "start_offset:", start_offset, "end_offset:", end_offset)
#         ## jump to bin and start extracting
#         output = ''
#         offset = 0
#         with open(self.filename, 'r', **self.kwargs_for_read) as f:
#             f.seek(self.indices[seqid][start_bin])
#             while True:
#                 tmp_seq = f.readline().rstrip()
#                 new_offset = offset + len(tmp_seq)
#                 ## exit loop if new molecule encountered
#                 if tmp_seq[0] == '>' or not tmp_seq: break
#                 ## if start position has been encountered, set start position as 0 (start of line)
#                 if output:
#                     tmp_start
#                 else:
#                     ## if start position within line, get start position
#                     if offset <= start_offset < new_offset:
#                         tmp_start = start_offset - offset
#                     ## else move to next line
#                     else:
#                         offset = new_offset
#                         continue
#                 ## extend output sequence accordingly
#                 ## if end position in line, get end position, slice, and break
#                 if end_offset < new_offset:
#                     print(end_offset, new_offset, offset)
#                     output += tmp_seq[tmp_start:end_offset - offset]
#                     break
#                 ## if end position not in line, set end position as end of line
#                 else:
#                     output += tmp_seq[tmp_start:]
#                 offset = new_offset
#         return Seq.Seq(output)
