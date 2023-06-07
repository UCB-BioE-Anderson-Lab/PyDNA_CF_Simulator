import pytest
from pydna_cf_simulator.polynucleotide import Polynucleotide
from pydna_cf_simulator.polynucleotide_to_dseqrecord import polynucleotide_to_dseqrecord
from pydna_cf_simulator.dseqrecord_to_polynucleotide import dseqrecord_to_polynucleotide

test_examples = [
    ('Blunt both ends', 'ATGCGCTGAC', '', '', False),
    ('Long top strand', 'ATGCGCTGAC', 'CGAA', '-TC', False),
    ('Long bottom strand', 'ATGCGCTGAC', '-CGAA', 'TC', False),
    ('Lefty blunt, righty 5\' extension', 'ATGCGCTGAC', '', '-TC', False),
    ('Lefty blunt, righty 3\' extension', 'ATGCGCTGAC', '', 'TC', False),
    ('Righty blunt, lefty 5\' extension', 'ATGCGCTGAC', 'CGAA', '', False),
    ('Righty blunt, lefty 3\' extension', 'ATGCGCTGAC', '-CGAA', '', False),
    ('Circular sequence with no overhangs', 'ATGCGCTGAC', '', '', True),
    ('Degenerate nucleotides', 'ATGCNRYKSWBVHDM', '', '', False),
]

exception_test_examples = [
    ('Circular sequence with overhangs (ext5)', 'ATGCGCTGAC', 'CGAA', '', True),
    ('Circular sequence with overhangs (ext3)', 'ATGCGCTGAC', '', 'TC', True),
    ('Invalid characters in sequence', 'ATGCGCTGACX', '', '', False),
    ('Whitespace characters', 'ATGC ATCG', '', '', False),
    ('Special characters', 'ATGC!', '', '', False),
    ('Numeric characters', 'ATGC4', '', '', False),
    ('Degenerate nucleotides in left extension', 'ATGC', 'NAAA', '', False),
    ('Degenerate nucleotides in righty extension', 'ATGC', '', 'NAAA', False),
    ('Whitespace characters in lefty extension', 'A TGC', 'ATGC', '', False),
    ('Whitespace characters in righty extension', 'ATGC', 'A TGC', '', False),
    ('Special characters in lefty extension', '!ATGC', 'ATCG', '', False),
    ('Special characters in rigty extension', 'ATGC', '!ATCG', '', False),
    ('Numeric characters in lefty extension', '1ATCG', 'ATCG', '', False),
    ('Numeric characters in righty extension', 'ATGC', '1ATCG', '', False),
    ('Empty sequence with non-empty extensions', '', 'ATGC', 'CGTA', False),
]

@pytest.mark.parametrize('description,sequence,ext_5,ext_3,is_circular', test_examples)
def test_round_trip_conversion(description, sequence, ext_5, ext_3, is_circular):
    # Create original Polynucleotide
    original_poly = Polynucleotide(sequence, ext_5, ext_3, True, is_circular, 'phosphate', 'phosphate')

    # Convert to Dseqrecord
    dseqrecord = polynucleotide_to_dseqrecord(original_poly)
    print(f'Description: {description}')
    print(f'Original Polynucleotide: {original_poly}')
    print(f'Converted Dseq overhang: {dseqrecord.seq.ovhg}')
    print(f'Forward strand: 5\'-{dseqrecord.seq.watson}-3\'')
    print(f'Reverse strand: 5\'-{dseqrecord.seq.crick}-3\'')

    # Convert back to Polynucleotide
    converted_poly = dseqrecord_to_polynucleotide(dseqrecord, 'phosphate', 'phosphate')
    print(f'Converted Polynucleotide: {converted_poly}')

    # Assert that the converted Polynucleotide is equal to the original
    assert converted_poly == original_poly
@pytest.mark.parametrize('description,sequence,ext_5,ext_3,is_circular', exception_test_examples)
def test_exception_cases(description, sequence, ext_5, ext_3, is_circular):
    with pytest.raises(Exception):
        original_poly = Polynucleotide(sequence, ext_5, ext_3, True, is_circular, 'phosphate', 'phosphate')
        dseqrecord = polynucleotide_to_dseqrecord(original_poly)
