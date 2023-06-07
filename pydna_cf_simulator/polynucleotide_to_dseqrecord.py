from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
from pydna.utils import rc


def polynucleotide_to_dseqrecord(poly):
    # Grab the 'sequences' portion of the Polynucleotide
    coding_strand = poly.sequence
    # Calculate the reverse complement
    complementary_strand = rc(coding_strand)

    # Handle 5' extension
    ext_5 = poly.ext5
    if not ext_5.startswith('-'):
        coding_strand = ext_5 + coding_strand
        overhang = len(ext_5)
    else:
        ext_5 = ext_5[1:]
        complementary_strand = complementary_strand + rc(ext_5)
        overhang = -len(ext_5)

    # Handle 3' extension
    ext_3 = poly.ext3
    if not ext_3.startswith('-'):
        ext_3 = rc(ext_3)
        complementary_strand = ext_3 + complementary_strand
    else:
        ext_3 = ext_3[1:]
        coding_strand += ext_3

    # Create Dseqrecord
    sequence = Dseq(coding_strand, complementary_strand, overhang)
    return Dseqrecord(sequence)
