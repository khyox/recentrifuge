"""Bio.SeqIO support for the "lmat output" (.out extension) file format

You are expected to use this module via the Bio.SeqIO functions.
This module is for reading and writing LMAT output files as SeqRecord
objects.  The code is partly inspired by Peter Cock FASTA Biopython module
"""

from typing import Tuple, TextIO, Iterator, Generator, Optional, Callable, cast, Any, Dict

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.Interfaces import SequenceWriter

__docformat__ = "restructuredtext en"


def simple_lmat_out_parser(handle: TextIO) -> Generator[Tuple[str, str, str, str, str], None, None]:
    """Generator function to iterate LMAT output records (as string tuples)

    For each record a tuple of five strings is returned:
    - the FASTA title line
    - the sequence
    - scoring statistics (average, standard deviation, # of k-mers)
    - list of taxid,score pairs
    - final LMAT taxid call with score and match type. Match types are
      MultiMatch (lowest common ancestor), DirectMatch (best match),
      NoMatch (no database match).
    """
    # Skip any text before the first record (e.g. blank lines, comments)
    while True:
        rawline = handle.readline()
        if rawline == "":
            return  # Premature end of file, or just empty?
        else:
            break
    while True:
        line = (rawline.rstrip('\n')).split('\t')
        try:
            (title, sequence, stats, lists, finalcall) = line
        except ValueError:
            return  # Stop Iteration
        yield title, sequence, stats, lists, finalcall
        rawline = handle.readline()


def lmat_out_iterator(handle: TextIO) -> Iterator[SeqRecord]:
    """Generator function to iterate LMAT output records (as SeqRecord objects)

    Arguments:
     - handle - input file
     - alphabet - optional alphabet
    """
    for (title, sequence,
         stats, candidates, finalcall) in simple_lmat_out_parser(handle):
        try:
            first_word = title.split(None, 1)[0]
        except IndexError:
            if title:  # Non-empty title but no words found is unexpected
                raise ValueError(f'Unexpected title format: {repr(title)}')
            # Should we use SeqRecord default for no ID?
            first_word = ""
        avg, std, kms = stats.split()
        statsdict: Dict[str, float | int] = {'averg': float(avg),
                                              'stdev': float(std),
                                              'kmers': int(kms)}
        final_taxid, final_score, final_match = finalcall.split()
        candids = candidates.split()
        candidict: Dict[str, float] = {}
        for i in range(0, len(candids), 2):
            candidict[candids[i]] = float(candids[i+1])
        # Create SeqRecord and update annotations; cast to Dict[str, Any]
        # to allow flexible value types (Biopython's _AnnotationsDict is strict)
        record = SeqRecord(Seq(sequence),
                           id=first_word, name=first_word, description=title)
        annot = cast(Dict[str, Any], record.annotations)
        annot['stats'] = stats
        annot['candidates'] = candidates
        annot['final_taxid'] = final_taxid
        annot['final_score'] = float(final_score)
        annot['final_match'] = final_match
        annot['statsdict'] = statsdict
        annot['candidict'] = candidict
        annot['tags'] = set()
        yield record


class LmatOutWriter(SequenceWriter):
    """Class to write LMAT output files"""
    
    record2title: Optional[Callable[[SeqRecord], str]]
    
    def __init__(self, handle: Any,
                 record2title: Optional[Callable[[SeqRecord], str]] = None) -> None:
        """Create a LMAT output writer

        Arguments:

         - handle - Handle to an output file, e.g. as returned
           by open(filename, 'w')
         - record2title - Optional function to return the text to be
           used for the title line of each record.  By default
           a combination of the record.id and record.description
           is used.  If the record.description starts with the
           record.id, then just the record.description is used.

        You can either use::

            handle = open(filename, 'w')
            writer = LmatOutWriter(handle)
            writer.write_file(myRecords)
            handle.close()

        Or, follow the sequential file writer system, for example::

            handle = open(filename, 'w')
            writer = LmatOutWriter(handle)
            writer.write_header() # does nothing for LMAT output files
            ...
            Multiple writer.write_record() and/or writer.write_records() calls
            ...
            writer.write_footer() # does nothing for LMAT out files
            handle.close()

        """
        SequenceWriter.__init__(self, handle)
        self.record2title = record2title

    def write_record(self, record: SeqRecord) -> None:
        """Write a single LMAT output record to the file"""
        # Validate writer state (replaces asserts for HPC code)
        if not getattr(self, '_header_written', True):
            raise RuntimeError('Header must be written before records')
        if getattr(self, '_footer_written', False):
            raise RuntimeError('Cannot write records after footer')
        self._record_written = True  # type: ignore[attr-defined]

        if self.record2title:
            title = self.clean(self.record2title(record))
        else:
            rec_id = self.clean(record.id or '')
            description = self.clean(record.description or '')
            if description and description.split(None, 1)[0] == rec_id:
                # The description includes the rec_id at the start
                title = description
            elif description:
                title = '%s %s' % (rec_id, description)
            else:
                title = rec_id

        # Validate title format (replaces asserts)
        if '\n' in title or '\r' in title:
            raise ValueError(f'Title contains newline characters: {repr(title)}')
        
        handle = cast(TextIO, self.handle)
        handle.write('%s\t' % title)

        data = self._get_seq_string(record)  # type: ignore[attr-defined]  # Catches sequence being None

        # Validate data format (replaces asserts)
        if '\n' in data or '\r' in data:
            raise ValueError(f'Sequence data contains newline characters')

        handle.write(data + '\t')
        handle.write(str(record.annotations['stats']) + '\t')
        handle.write(str(record.annotations['candidates']) + '\t')
        handle.write(str(record.annotations['finalcall']) + '\n')
