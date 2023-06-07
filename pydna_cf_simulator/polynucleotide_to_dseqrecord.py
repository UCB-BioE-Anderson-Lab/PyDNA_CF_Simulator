import re
from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
from pydna.utils import rc


def polynucleotide_to_dseqrecord(poly):
    # Convert the sequence and extensions to uppercase
    sequence = poly.sequence.upper()
    ext5 = poly.ext5.upper() if poly.ext5 else ''
    ext3 = poly.ext3.upper() if poly.ext3 else ''
    
    # Check the sequence data
    if not re.match(r'^[ATCGNRKYSWBVHDM]+$', poly.sequence):
        raise ValueError("Invalid characters in sequence. Sequence must only contain the characters ATCGNRKYSWBVHDM.")

    if poly.ext5 and not re.match(r'^(-)?[ATCG]+$', poly.ext5):
        raise ValueError("Invalid characters in ext5. ext5 must only contain the characters ATCG and optionally start with a '-'.")

    if poly.ext3 and not re.match(r'^(-)?[ATCG]+$', poly.ext3):
        raise ValueError("Invalid characters in ext3. ext3 must only contain the characters ATCG and optionally start with a '-'.")

    # Check that circular sequences have no sticky ends
    if poly.is_circular and (poly.ext5 or poly.ext3):
        raise Exception("Circular Polynucleotide cannot have overhangs")

    # Grab the 'sequences' portion of the Polynucleotide
    coding_strand = poly.sequence
    
    # Calculate the reverse complement
    complementary_strand = rc(coding_strand)

    # Handle 5' extension
    ext_5 = poly.ext5
    if not ext_5.startswith('-'):
        coding_strand = ext_5 + coding_strand
        overhang = -len(ext_5)
    else:
        ext_5 = ext_5[1:]
        complementary_strand = complementary_strand + rc(ext_5)
        overhang = len(ext_5)

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

    # Loop it if it's circular
    if poly.is_circular:
        sequence = sequence.looped()

    return Dseqrecord(sequence)
