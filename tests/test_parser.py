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
    assert 'Fragment1' in cf.steps[0].inputs
    assert 'Fragment2' in cf.steps[0].inputs
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

def test_parse_sequences():
    cf = parse_CF_shorthand('P6libF CCAAAGGTCTCATTATANNNNNNNNNNNNNNNNTGTCAANNNNGAACCCAGGACTCCTCGAAGTCGTTCTTAAGACAAC')
    assert isinstance(cf, ConstructionFile)
    assert len(cf.sequences) == 1
    assert isinstance(cf.sequences['P6libF'], Polynucleotide)
    assert cf.sequences['P6libF'].sequence == 'CCAAAGGTCTCATTATANNNNNNNNNNNNNNNNTGTCAANNNNGAACCCAGGACTCCTCGAAGTCGTTCTTAAGACAAC'
