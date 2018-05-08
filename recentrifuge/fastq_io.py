"""Bio.SeqIO quick support for FASTQ files

You are expected to use this module via the Bio.SeqIO functions.
This module is for reading and writing FASTQ output files as SeqRecord
objects, but omitting some checks included in the Biopython method by Peter
Cock. These checks were very useful in the "olden times" but, currently,
with huge FASTQ files using standardized PHRED quality scores, they can be
omitted, which improves the code performance.

"""

from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.Interfaces import SequentialSequenceWriter
from Bio.SeqIO.QualityIO import FastqGeneralIterator

__docformat__ = "restructuredtext en"


def quick_fastq_iterator(handle, alphabet=single_letter_alphabet):
    """Parse Illumina 1.3 to 1.7 FASTQ files without decoding the qualities
    to improve the performance.
    """
    for title, sequence, quality in FastqGeneralIterator(handle):
        first_word = title.split()[0]
        yield SeqRecord(Seq(sequence, alphabet),
                        id=first_word, name=first_word, description=title,
                        annotations={'quality': quality})


class QuickFastqWriter(SequentialSequenceWriter):
    """Class to write standard FASTQ format files with sequences
    previously read by QuickFastqIterator function.

    Though you can use this class directly, you are strongly encouraged
    to use the Bio.SeqIO.write() function instead, via the format name
    "quickfastq". For example, this code reads a standard FASTQ file
    and copy it:

    >>> from Bio import SeqIO
    >>> record_iterator = SeqIO.parse("example.fastq", "quickfastq")
    >>> with open("copy.fastq", "w") as out_handle:
    ...     SeqIO.write(record_iterator, out_handle, "quickfastq")
    """

    def write_record(self, record):
        """Quickly write a single FASTQ record to the file."""

        self.handle.write(f'@{record.description}\n{str(record.seq)}\n+'
                          f'\n{record.annotations["quality"]}\n')
