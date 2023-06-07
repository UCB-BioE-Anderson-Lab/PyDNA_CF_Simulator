import pytest
from pydna_cf_simulator.simulate_CF import simulate_CF
from pydna_cf_simulator.construction_file import ConstructionFile, PCR, Digest, Ligate, GoldenGate, Gibson, Transform
from pydna_cf_simulator.polynucleotide import Polynucleotide

def test_simulate_PCR():
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
    assert result['product'] == expected
    
def test_simulate_PCR_with_modifications():
    # Define sequences
    forward = Polynucleotide('CCGCAACACACTTAACCTTG', '', '', False, False, 'phosphate', '')
    reverse = Polynucleotide('GTGGTTGTGGCCGGTCAAATC', '', '', False, False, 'phosphate', '')
    template = Polynucleotide('CCGCAACACACTTAACCTTGGCGTCGGGATACGTACATTGGAGAACGGTTGGCTGTACGGACTTAATACTTTTTATGATAATGATTTGACCGGCCACAACCACCG', '', '', True, False, 'phosphate', 'phosphate')

    # Define steps
    step = PCR('forward', 'reverse', 'template', 'product')

    # Define ConstructionFile
    cf = ConstructionFile([step], {'forward': forward, 'reverse': reverse, 'template': template})

    # Simulate ConstructionFile
    result = simulate_CF(cf)

    # Assert that the product is correct
    expected = Polynucleotide('CCGCAACACACTTAACCTTGGCGTCGGGATACGTACATTGGAGAACGGTTGGCTGTACGGACTTAATACTTTTTATGATAATGATTTGACCGGCCACAACCAC', '', '', True, False, 'phosphate', 'phosphate')
    assert result['product'] == expected
    
def test_simulate_PCR_with_extensions():
    # Define sequences
    forward = Polynucleotide('CCATAGGATCCCCGCAACACACTTAACCTTG', '', '', False, False, 'hydroxyl', 'hydroxyl')
    reverse = Polynucleotide('GAGTCGAATTCGTGGTTGTGGCCGGTCAAATC', '', '', False, False, 'hydroxyl', 'hydroxyl')
    template = Polynucleotide('CCGCAACACACTTAACCTTGGCGTCGGGATACGTACATTGGAGAACGGTTGGCTGTACGGACTTAATACTTTTTATGATAATGATTTGACCGGCCACAACCACCG', '', '', True, False, 'phosphate', 'phosphate')

    # Define steps
    step = PCR('forward', 'reverse', 'template', 'product')

    # Define ConstructionFile
    cf = ConstructionFile([step], {'forward': forward, 'reverse': reverse, 'template': template})

    # Simulate ConstructionFile
    result = simulate_CF(cf)

    # Assert that the product is correct
    expected = Polynucleotide('CCATAGGATCCCCGCAACACACTTAACCTTGGCGTCGGGATACGTACATTGGAGAACGGTTGGCTGTACGGACTTAATACTTTTTATGATAATGATTTGACCGGCCACAACCACGAATTCGACTC', '', '', True, False, 'hydroxyl', 'hydroxyl')
    assert result['product'] == expected
    
def test_simulate_circular_template_PCR():
    # Define sequences
    forward = Polynucleotide('CCGCAACACACTTAACCTTG', '', '', False, False, 'hydroxyl', 'hydroxyl')
    reverse = Polynucleotide('GTGGTTGTGGCCGGTCAAATC', '', '', False, False, 'hydroxyl', 'hydroxyl')
    template = Polynucleotide('CCGCAACACACTTAACCTTGGCGTCGGGATACGTACATTGGAGAACGGTTGGCTGTACGGACTTAATACTTTTTATGATAATGATTTGACCGGCCACAACCACCG', '', '', True, True, '', '')

    # Define steps
    step = PCR('forward', 'reverse', 'template', 'product')

    # Define ConstructionFile
    cf = ConstructionFile([step], {'forward': forward, 'reverse': reverse, 'template': template})

    # Simulate ConstructionFile
    result = simulate_CF(cf)

    # Assert that the product is correct
    expected = Polynucleotide('CCGCAACACACTTAACCTTGGCGTCGGGATACGTACATTGGAGAACGGTTGGCTGTACGGACTTAATACTTTTTATGATAATGATTTGACCGGCCACAACCAC', '', '', True, False, 'hydroxyl', 'hydroxyl')
    assert result['product'] == expected
    
def test_simulate_permuted_circular_template_PCR():
    # Define sequences
    forward = Polynucleotide('CCGCAACACACTTAACCTTG', '', '', False, False, 'hydroxyl', 'hydroxyl')
    reverse = Polynucleotide('GTGGTTGTGGCCGGTCAAATC', '', '', False, False, 'hydroxyl', 'hydroxyl')
    template = Polynucleotide('GTACGGACTTAATACTTTTTATGATAATGATTTGACCGGCCACAACCACCGCCGCAACACACTTAACCTTGGCGTCGGGATACGTACATTGGAGAACGGTTGGCT', '', '', True, True, '', '')

    # Define steps
    step = PCR('forward', 'reverse', 'template', 'product')

    # Define ConstructionFile
    cf = ConstructionFile([step], {'forward': forward, 'reverse': reverse, 'template': template})

    # Simulate ConstructionFile
    result = simulate_CF(cf)

    # Assert that the product is correct
    expected = Polynucleotide('CCGCAACACACTTAACCTTGGCGTCGGGATACGTACATTGGAGAACGGTTGGCTGTACGGACTTAATACTTTTTATGATAATGATTTGACCGGCCACAACCAC', '', '', True, False, 'hydroxyl', 'hydroxyl')
    assert result['product'] == expected
    
def test_simulate_digest():
    # Define sequences
    sequence = Polynucleotide('GAGTCGAATTCATACGAGGGATCCAATCG', '', '', True, False, 'hydroxyl', 'hydroxyl')

    # Define steps
    step = Digest('sequence', 'EcoRI,BamHI', 1, 'product')

    # Define ConstructionFile
    cf = ConstructionFile([step], {'sequence': sequence})

    # Simulate ConstructionFile
    result = simulate_CF(cf)

    # Assert that the product is correct
    expected = Polynucleotide('CATACGAGG', 'AATT', 'GATC', True, False, 'phosphate', 'phosphate')
    assert result['product'] == expected
    
def test_simulate_digest_circular1(): # Fragment in middle of sequence (index 1)
    # Define sequences
    sequence = Polynucleotide('GAGTCGAATTCATACGAGGGATCCAATCG', '', '', True, True, '', '')

    # Define steps
    step = Digest('sequence', 'BamHI,EcoRI', 1, 'product')

    # Define ConstructionFile
    cf = ConstructionFile([step], {'sequence': sequence})

    # Simulate ConstructionFile
    result = simulate_CF(cf)

    # Assert that the product is correct
    expected = Polynucleotide('CATACGAGG', 'AATT', 'GATC', True, False, 'phosphate', 'phosphate')
    assert result['product'] == expected

def test_simulate_digest_circular2():  # Input DNA starts with restriction site (index 1)
    # Define sequences
    sequence = Polynucleotide('GAATTCATACGAGGGATCCAATCGGAGTC', '', '', True, True, '', '')

    # Define steps
    step = Digest('sequence', 'BamHI,EcoRI', 1, 'product')

    # Define ConstructionFile
    cf = ConstructionFile([step], {'sequence': sequence})

    # Simulate ConstructionFile
    result = simulate_CF(cf)

    # Assert that the product is correct
    expected = Polynucleotide('CATACGAGG', 'AATT', 'GATC', True, False, 'phosphate', 'phosphate')
    assert result['product'] == expected

def test_simulate_digest_circular3():  # Input DNA starts with sticky end  (index 0)
    # Define sequences
    sequence = Polynucleotide('AATTCATACGAGGGATCCAATCGGAGTCG', '', '', True, True, '', '')

    # Define steps
    step = Digest('sequence', 'BamHI,EcoRI', 0, 'product')

    # Define ConstructionFile
    cf = ConstructionFile([step], {'sequence': sequence})

    # Simulate ConstructionFile
    result = simulate_CF(cf)

    # Assert that the product is correct
    expected = Polynucleotide('CATACGAGG', 'AATT', 'GATC', True, False, 'phosphate', 'phosphate')
    assert result['product'] == expected

def test_simulate_ligate_linear():
    # Define sequences
    polyL = Polynucleotide('GTATACCCA', '', 'GATC', True, False, 'hydroxyl', 'phosphate')
    polyR = Polynucleotide('CATATGCAG', 'GATC', '', True, False, 'phosphate', 'hydroxyl')

    # Define steps
    step = Ligate(['polyL', 'polyR'], 'product')

    # Define ConstructionFile
    cf = ConstructionFile([step], {'polyL': polyL, 'polyR': polyR})

    # Simulate ConstructionFile
    result = simulate_CF(cf)

    # Assert that the product is correct
    expected = Polynucleotide('GTATACCCAGATCCATATGCAG', '', '', True, False, 'hydroxyl', 'hydroxyl')
    assert result['product'] == expected

def test_simulate_ligate_circular():
    # Define sequences
    polyL = Polynucleotide('GTATACCCA', 'AATT', 'GATC', True, False, 'phosphate', 'phosphate')
    polyR = Polynucleotide('CATATGCAG', 'GATC', 'AATT', True, False, 'phosphate', 'phosphate')

    # Define steps
    step = Ligate(['polyL', 'polyR'], 'product')

    # Define ConstructionFile
    cf = ConstructionFile([step], {'polyL': polyL, 'polyR': polyR})

    # Simulate ConstructionFile
    result = simulate_CF(cf)

    # Assert that the product is correct
    expected = Polynucleotide('AATTGTATACCCAGATCCATATGCAG', '', '', True, True, '', '')
    assert result['product'] == expected

def test_simulate_ligate_single_circular():
    # Define sequences
    poly = Polynucleotide('GTATACCCA', 'GATC', 'GATC', True, False, 'phosphate', 'phosphate')

    # Define steps
    step = Ligate(['poly'], 'product')

    # Define ConstructionFile
    cf = ConstructionFile([step], {'poly': poly})

    # Simulate ConstructionFile
    result = simulate_CF(cf)

    # Assert that the product is correct
    expected = Polynucleotide('GATCGTATACCCA', '', '', True, True, '', '')
    assert result['product'] == expected
