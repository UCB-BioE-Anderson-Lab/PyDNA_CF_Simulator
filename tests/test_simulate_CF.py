import pytest
from pydna_cf_simulator.simulate_CF import simulate_CF
from pydna_cf_simulator.construction_file import ConstructionFile, PCR, Digest, Ligate, GoldenGate, Gibson, Transform
from pydna_cf_simulator.polynucleotide import Polynucleotide


def test_simulate_CF_PCR():
    # Define sequences
    forward = Polynucleotide('CCGCAACACACTTAACCTTG', '', '', False, False, 'hydroxyl', 'hydroxyl')
    reverse = Polynucleotide('GTGGTTGTGGCCGGTCAAATC', '', '', False, False, 'hydroxyl', 'hydroxyl')
    template = Polynucleotide('CCGCAACACACTTAACCTTGGCGTCGGGATACGTACATTGGAGAACGGTTGGCTGTACGGACTTAATACTTTTTATGATAATGATTTGACCGGCCACAACCACCG', '', '', True, False, 'phosphate', 'phosphate')

    # Define steps
    step = PCR('forward', 'reverse', 'template', 'product')

    # Define ConstructionFile
    cf = ConstructionFile([step], {'forward': forward, 'reverse': reverse, 'template': template})

    # Simulate ConstructionFile
    result = simulate_CF(cf)

    # Assert that the product is correct
    expected = Polynucleotide('CCGCAACACACTTAACCTTGGCGTCGGGATACGTACATTGGAGAACGGTTGGCTGTACGGACTTAATACTTTTTATGATAATGATTTGACCGGCCACAACCAC', '', '', True, False, 'hydroxyl', 'hydroxyl')
