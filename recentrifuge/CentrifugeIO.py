"""Bio.SeqIO support for Centrifuge output

You are expected to use this module via the Bio.SeqIO functions."""

from Bio.Alphabet import single_letter_alphabet
from Bio.Seq import UnknownSeq
from Bio.SeqRecord import SeqRecord

__docformat__ = "restructuredtext en"


def simple_out_parser(handle):
    """Generator to iterate Centrifuge output (as string tuples)

    For each record a tuple of five strings is returned:
    - readID: read ID from a raw sequencing read.
    - seqID: sequence ID of the genomic sequence, where the read is
        classified.
    - taxID: taxonomic ID of the genomic sequence in the second column.
    - score: score for the classification (weighted sum of hits).
    - 2ndBestScore: score for the next best classification.
    - hitLength: a pair of two numbers:
        (1) an approximate number of base pairs of the read that match
            the genomic sequence,
        (2) the length of a read or the combined length of mate pairs.
    - queryLength: same pair of numbers as in the previous column.
    - numMatches: number of classifications for this read, indicating
        how many assignments were made
    """
    # Skip any text before the first record (e.g. blank lines, comments)
    while True:
        raw_line = handle.readline()
        if raw_line == '':
            return  # Premature end of file, or just empty?
        elif 'readID' in raw_line:
            pass  # title line
        else:
            break
    while True:
        line = (raw_line.rstrip('\n')).split('\t')
        try:
            (read_id, seq_id, tax_id, score, second_score, hit_length,
             query_length, num_matches) = line
        except ValueError:
            return  # Stop Iteration
        yield (read_id, seq_id, tax_id, score, second_score,
               hit_length, query_length, num_matches)
        raw_line = handle.readline()


def cfg_out_iterator(handle, alphabet=single_letter_alphabet):
    """Generator to iterate Centrifuge output (as SeqRecord objects)

    Arguments:
     - handle - input file
     - alphabet - optional alphabet
    """
    for (read_id, seq_id, tax_id, score, second_score,
         hit_length, query_length, num_matches) in simple_out_parser(handle):
        try:
            first_word = read_id.split(None, 1)[0]
        except IndexError:
            assert not read_id, repr(read_id)
            # Should we use SeqRecord default for no ID?
            first_word = ""
        # From Centrifuge score get the "single hit equivalent length"
        try:
            adapted_score = float(score) ** 0.5 + 15
        except ValueError:
            print(f'Error parsing score ({score}) for taxid {tax_id}'
                  f' in {handle}...')
            raise
        try:
            adapted_2nd_score = float(second_score) ** 0.5 + 15
        except ValueError:
            print(f'Error parsing score ({second_score}) for taxid {tax_id}'
                  f' in {handle}...')
            raise
        yield SeqRecord(UnknownSeq(0, alphabet),
                        id=first_word,
                        name=first_word,
                        description=read_id,
                        dbxrefs=[seq_id],
                        annotations={'taxID': tax_id,
                                     'score': adapted_score,
                                     '2ndBestScore': adapted_2nd_score,
                                     'hitLength': hit_length,
                                     'queryLength': query_length,
                                     'numMatches': int(num_matches),
                                     })
