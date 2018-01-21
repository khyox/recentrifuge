"""Bio.SeqIO support for the "lmat output" (.out extension) file format

You are expected to use this module via the Bio.SeqIO functions.
This module is for reading and writing LMAT output files as SeqRecord
objects.  The code is partly inspired by Peter Cock FASTA Biopython module

"""

from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.Interfaces import SequentialSequenceWriter

__docformat__ = "restructuredtext en"


def simple_lmat_out_parser(handle):
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


def lmat_out_iterator(handle, alphabet=single_letter_alphabet):
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
            assert not title, repr(title)
            # Should we use SeqRecord default for no ID?
            first_word = ""
        avg, std, kms = stats.split()
        statsdict = {'averg':float(avg), 'stdev':float(std), 'kmers':int(kms)}
        final_taxid, final_score, final_match = finalcall.split()
        candids = candidates.split()
        candidict = {}
        for i in range(0, len(candids), 2):
            candidict[candids[i]] = float(candids[i+1])
        yield SeqRecord(Seq(sequence, alphabet),
                        id=first_word, name=first_word, description=title,
                        annotations={'stats': stats,
                                     'candidates': candidates,
                                     'final_taxid': final_taxid,
                                     'final_score': float(final_score),
                                     'final_match': final_match,
                                     'statsdict': statsdict,
                                     'candidict': candidict,
                                     'tags': set()
                                     })


class LmatOutWriter(SequentialSequenceWriter):
    """Class to write LMAT output files"""
    def __init__(self, handle, record2title=None):
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
        SequentialSequenceWriter.__init__(self, handle)
        self.record2title = record2title

    def write_record(self, record):
        """Write a single LMAT output record to the file"""
        assert self._header_written
        assert not self._footer_written
        self._record_written = True

        if self.record2title:
            title = self.clean(self.record2title(record))
        else:
            rec_id = self.clean(record.id)
            description = self.clean(record.description)
            if description and description.split(None, 1)[0] == rec_id:
                # The description includes the rec_id at the start
                title = description
            elif description:
                title = '%s %s' % (rec_id, description)
            else:
                title = rec_id

        assert '\n' not in title
        assert '\r' not in title
        self.handle.write('%s\t' % title)

        data = self._get_seq_string(record)  # Catches sequence being None

        assert '\n' not in data
        assert '\r' not in data

        self.handle.write(data + '\t')
        self.handle.write(record.annotations['stats'] + '\t')
        self.handle.write(record.annotations['candidates'] + '\t')
        self.handle.write(record.annotations['finalcall'] + '\n')
