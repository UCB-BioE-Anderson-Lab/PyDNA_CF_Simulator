import pytest
from pydna_cf_simulator.parse_CF_shorthand import parse_CF_shorthand
from pydna_cf_simulator.construction_file import ConstructionFile, PCR, GoldenGate, Transform, Digest, Ligate, Gibson
from pydna_cf_simulator.polynucleotide import Polynucleotide


def test_parse_pcr():
    cf = parse_CF_shorthand('PCR ForwardPrimer ReversePrimer Template ProductName')
    assert isinstance(cf, ConstructionFile)
    assert len(cf.steps) == 1
    assert isinstance(cf.steps[0], PCR)
    assert cf.steps[0].forward_oligo == 'ForwardPrimer'
    assert cf.steps[0].reverse_oligo == 'ReversePrimer'
    assert cf.steps[0].template == 'Template'
    assert cf.steps[0].output == 'ProductName'

def test_parse_golden_gate():
    cf = parse_CF_shorthand('GoldenGate Fragment1 Fragment2 BsaI ProductName')
    assert isinstance(cf, ConstructionFile)
    assert len(cf.steps) == 1
    assert isinstance(cf.steps[0], GoldenGate)
    assert 'Fragment1' in cf.steps[0].dnas
    assert 'Fragment2' in cf.steps[0].dnas
    assert cf.steps[0].enzyme == 'BsaI'
    assert cf.steps[0].output == 'ProductName'

def test_parse_digest():
    cf = parse_CF_shorthand('Digest DNA EcoRI,BamHI 1 ProductName')
    assert isinstance(cf, ConstructionFile)
    assert len(cf.steps) == 1
    assert isinstance(cf.steps[0], Digest)
    assert cf.steps[0].dna == 'DNA'
    assert 'EcoRI' in cf.steps[0].enzymes
    assert 'BamHI' in cf.steps[0].enzymes
    assert cf.steps[0].fragSelect == 1
    assert cf.steps[0].output == 'ProductName'

def test_parse_ligate():
    cf = parse_CF_shorthand('Ligate Fragment1 Fragment2 ProductName')
    assert isinstance(cf, ConstructionFile)
    assert len(cf.steps) == 1
    assert isinstance(cf.steps[0], Ligate)
    assert 'Fragment1' in cf.steps[0].dnas
    assert 'Fragment2' in cf.steps[0].dnas
    assert cf.steps[0].output == 'ProductName'

def test_parse_gibson():
    cf = parse_CF_shorthand('Gibson Fragment1 Fragment2 ProductName')
    assert isinstance(cf, ConstructionFile)
    assert len(cf.steps) == 1
    assert isinstance(cf.steps[0], Gibson)
    assert 'Fragment1' in cf.steps[0].dnas
    assert 'Fragment2' in cf.steps[0].dnas
    assert cf.steps[0].output == 'ProductName'

def test_parse_transform():
    cf = parse_CF_shorthand('Transform DNA Strain Amp ProductName')
    assert isinstance(cf, ConstructionFile)
    assert len(cf.steps) == 1
    assert isinstance(cf.steps[0], Transform)
    assert cf.steps[0].dna == 'DNA'
    assert cf.steps[0].strain == 'Strain'
    assert cf.steps[0].antibiotics == ['Amp']
    assert cf.steps[0].temperature == None
    assert cf.steps[0].output == 'ProductName'

def test_parse_transform_with_temperature():
    cf = parse_CF_shorthand('Transform DNA Strain Amp 37 ProductName')
    assert isinstance(cf, ConstructionFile)
    assert len(cf.steps) == 1
    assert isinstance(cf.steps[0], Transform)
    assert cf.steps[0].dna == 'DNA'
    assert cf.steps[0].strain == 'Strain'
    assert cf.steps[0].antibiotics == ['Amp']
    assert cf.steps[0].temperature == 37
    assert cf.steps[0].output == 'ProductName'

def test_parse_multiple_steps():
    cf = parse_CF_shorthand('PCR ForwardPrimer ReversePrimer Template ProductName\nGoldenGate Fragment1 Fragment2 BsaI ProductName')
    assert isinstance(cf, ConstructionFile)
    assert len(cf.steps) == 2
    assert isinstance(cf.steps[0], PCR)
    assert isinstance(cf.steps[1], GoldenGate)

def test_parse_oligo_sequence():
    cf = parse_CF_shorthand('P6libF CCAAAGGTCTCATTATANNNNNNNNNNNNNNNNTGTCAANNNNGAACCCAGGACTCCTCGAAGTCGTTCTTAAGACAAC')
    assert isinstance(cf, ConstructionFile)
    assert len(cf.sequences) == 1
    assert isinstance(cf.sequences['P6libF'], Polynucleotide)
    assert cf.sequences['P6libF'].sequence == 'CCAAAGGTCTCATTATANNNNNNNNNNNNNNNNTGTCAANNNNGAACCCAGGACTCCTCGAAGTCGTTCTTAAGACAAC'
    assert cf.sequences['P6libF'].is_double_stranded == False
    assert cf.sequences['P6libF'].is_circular == False
    
def test_parse_plasmid_sequence():
    long_sequence = 'A' * 1000
    cf = parse_CF_shorthand(f'pPlas1 {long_sequence}')
    assert isinstance(cf, ConstructionFile)
    assert len(cf.sequences) == 1
    assert isinstance(cf.sequences['pPlas1'], Polynucleotide)
    assert cf.sequences['pPlas1'].sequence == long_sequence
    assert cf.sequences['pPlas1'].is_double_stranded == True
    assert cf.sequences['pPlas1'].is_circular == True

def test_parse_multiple_sequences():
    cf = parse_CF_shorthand('oligo P6libF CCAAAGGTCTCATTATANNNNNNNNNNNNNNNNTGTCAANNNNGAACCCAGGACTCCTCGAAGTCGTTCTTAAGACAAC\noligo P6libR CAGTTGGTCTCAATAATNNNNNNANNNNGTTAGTATTTCTCCTCGTCTACGGTTAACTGATACTC')
    assert isinstance(cf, ConstructionFile)
    assert len(cf.sequences) == 2
    assert isinstance(cf.sequences['P6libF'], Polynucleotide)
    assert cf.sequences['P6libF'].sequence == 'CCAAAGGTCTCATTATANNNNNNNNNNNNNNNNTGTCAANNNNGAACCCAGGACTCCTCGAAGTCGTTCTTAAGACAAC'
    assert cf.sequences['P6libF'].is_double_stranded == False
    assert cf.sequences['P6libF'].is_circular == False
    assert isinstance(cf.sequences['P6libR'], Polynucleotide)
    assert cf.sequences['P6libR'].sequence == 'CAGTTGGTCTCAATAATNNNNNNANNNNGTTAGTATTTCTCCTCGTCTACGGTTAACTGATACTC'
    assert cf.sequences['P6libR'].is_double_stranded == False
    assert cf.sequences['P6libR'].is_circular == False

def test_parse_oligo():
    cf = parse_CF_shorthand('oligo oligoname ATCAGATTT')
    assert isinstance(cf, ConstructionFile)
    assert len(cf.sequences) == 1
    assert isinstance(cf.sequences['oligoname'], Polynucleotide)
    assert cf.sequences['oligoname'].sequence == 'ATCAGATTT'
    assert cf.sequences['oligoname'].is_double_stranded == False
    assert cf.sequences['oligoname'].is_circular == False

def test_parse_plasmid():
    cf = parse_CF_shorthand('plasmid plasmidname ATCAGATTT')
    assert isinstance(cf, ConstructionFile)
    assert len(cf.sequences) == 1
    assert isinstance(cf.sequences['plasmidname'], Polynucleotide)
    assert cf.sequences['plasmidname'].sequence == 'ATCAGATTT'
    assert cf.sequences['plasmidname'].is_double_stranded == True
    assert cf.sequences['plasmidname'].is_circular == True

def test_parse_dsdna():
    cf = parse_CF_shorthand('dsdna dsdnaname ATCAGATTT')
    assert isinstance(cf, ConstructionFile)
    assert len(cf.sequences) == 1
    assert isinstance(cf.sequences['dsdnaname'], Polynucleotide)
    assert cf.sequences['dsdnaname'].sequence == 'ATCAGATTT'
    assert cf.sequences['dsdnaname'].is_double_stranded == True
    assert cf.sequences['dsdnaname'].is_circular == False

def test_parse_invalid_sequence_format():
    with pytest.raises(ValueError) as e:
        parse_CF_shorthand('P6libF CCAAAGGTCTCATTATANNNNNNNNNNNNNNNNTGTCAXXXGAACCCAGGACTCCTCGAAGTCGTTCTTAAGACAAC')
    assert isinstance(e.value, ValueError)

def test_parse_invalid_number_of_arguments_pcr():
    with pytest.raises(ValueError) as e:
        parse_CF_shorthand('PCR ForwardPrimer ReversePrimer')
    assert isinstance(e.value, ValueError)

def test_parse_invalid_number_of_arguments_digest():
    with pytest.raises(ValueError) as e:
        parse_CF_shorthand('Digest DNA EcoRI,BamHI 1 ProductName ExtraArg')
    assert isinstance(e.value, ValueError)

def test_parse_invalid_number_of_arguments_ligate():
    with pytest.raises(ValueError) as e:
        parse_CF_shorthand('Ligate Fragment1 ')
    assert isinstance(e.value, ValueError)

def test_parse_invalid_number_of_arguments_golden_gate():
    with pytest.raises(ValueError) as e:
        parse_CF_shorthand('GoldenGate Fragment1')
    assert isinstance(e.value, ValueError)

def test_parse_invalid_number_of_arguments_gibson():
    with pytest.raises(ValueError) as e:
        parse_CF_shorthand('Gibson Fragment1')
    assert isinstance(e.value, ValueError)

def test_parse_invalid_number_of_arguments_transform():
    with pytest.raises(ValueError) as e:
        parse_CF_shorthand('Transform DNA Strain Amp')
    assert isinstance(e.value, ValueError)

def test_parse_invalid_temperature_transform():
    with pytest.raises(ValueError) as e:
        parse_CF_shorthand('Transform DNA Strain Amp InvalidTemperature ProductName')
    assert isinstance(e.value, ValueError)

def test_parse_unrecognized_enzyme_digest():
    with pytest.raises(ValueError) as e:
        parse_CF_shorthand('Digest DNA EcoRI,NotAnEnzyme,XbaI 1 ProductName')
    assert isinstance(e.value, ValueError)

def test_parse_invalid_sequence_number():
    with pytest.raises(ValueError) as e:
        parse_CF_shorthand('Digest DNA EcoRI,BamHI -1 ProductName')
    assert isinstance(e.value, ValueError)
