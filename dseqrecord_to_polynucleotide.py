from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
from pydna.utils import rc
from .polynucleotide import Polynucleotide


def dseqrecord_to_polynucleotide(dseqrecord, mod_ext5, mod_ext3):
    # Extract the sequence and its complement
    sequence = dseqrecord.seq.watson
    complement = rc(sequence)

    # Determine the ext_5
    overhang = dseqrecord.seq.ovhg
    if overhang > 0:
        ext_5 = sequence[:overhang]
        sequence = sequence[overhang:]
    elif overhang < 0:
        ext_5 = '-' + sequence[overhang:]
        sequence = sequence[:overhang]
    else:
        ext_5 = ''

    # Determine the ext_3
    length_difference = len(sequence) - len(complement)
    if length_difference > 0 and overhang >= 0:
        ext_3 = '-' + sequence[-length_difference:]
        sequence = sequence[:-length_difference]
    elif length_difference < 0 and overhang < 0:
        ext_3 = complement[:length_difference]
        complement = complement[length_difference:]
    else:
        ext_3 = ''

    # Create Polynucleotide
    polynucleotide = Polynucleotide(sequence, ext_5, ext_3, True, dseqrecord.seq.circular(), mod_ext5, mod_ext3)
    return polynucleotide