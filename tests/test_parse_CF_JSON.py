import pytest
import json

from pydna_cf_simulator.parse_CF_JSON import parse_CF_JSON
from pydna_cf_simulator.construction_file import ConstructionFile, PCR, Digest, Ligate, GoldenGate, Gibson, Transform
from pydna_cf_simulator.polynucleotide import Polynucleotide

def test_pcr():
    json_string = '''{
        "steps": [
            {
                "operation": "PCR",
                "output": "P6",
                "forward_oligo": "P6libF",
                "reverse_oligo": "P6libR",
                "template": "pTP1"
            }
        ],
        "sequences": {}
    }'''

    cf = parse_CF_JSON(json_string)

    assert len(cf.steps) == 1
    assert cf.steps[0].operation == "PCR"
    assert cf.steps[0].output == "P6"
    assert cf.steps[0].forward_oligo == "P6libF"
    assert cf.steps[0].reverse_oligo == "P6libR"
    assert cf.steps[0].template == "pTP1"
    
def test_pcr_with_product_size():
    json_string = '''{
        "steps": [
            {
                "operation": "PCR",
                "output": "P6",
                "forward_oligo": "P6libF",
                "reverse_oligo": "P6libR",
                "template": "pTP1",
                "product_size": 3583
            }
        ],
        "sequences": {}
    }'''

    cf = parse_CF_JSON(json_string)

    assert len(cf.steps) == 1
    assert cf.steps[0].operation == "PCR"
    assert cf.steps[0].output == "P6"
    assert cf.steps[0].forward_oligo == "P6libF"
    assert cf.steps[0].reverse_oligo == "P6libR"
    assert cf.steps[0].template == "pTP1"
    assert cf.steps[0].product_size == 3583

def test_digest():
    json_string = '''{
        "steps": [
            {
                "operation": "Digest",
                "dna": "pTP1",
                "enzymes": ["EcoRI"],
                "fragSelect": 3583,
                "output": "digested_dna"
            }
        ],
        "sequences": {}
    }'''

    cf = parse_CF_JSON(json_string)

    assert len(cf.steps) == 1
    assert cf.steps[0].operation == "Digest"
    assert cf.steps[0].dna == "pTP1"
    assert cf.steps[0].enzymes == ["EcoRI"]
    assert cf.steps[0].fragSelect == 3583
    assert cf.steps[0].output == "digested_dna"

def test_ligate():
    json_string = '''{
        "steps": [
            {
                "operation": "Ligate",
                "dnas": ["Fragment1", "Fragment2"],
                "output": "LigatedProduct"
            }
        ],
        "sequences": {}
    }'''

    cf = parse_CF_JSON(json_string)

    assert len(cf.steps) == 1
    assert cf.steps[0].operation == "Ligate"
    assert cf.steps[0].dnas == ["Fragment1", "Fragment2"]
    assert cf.steps[0].output == "LigatedProduct"

def test_golden_gate():
    json_string = '''{
        "steps": [
            {
                "operation": "GoldenGate",
                "inputs": ["Fragment1", "Fragment2"],
                "enzyme": "BsaI",
                "output": "ProductName"
            }
        ],
        "sequences": {}
    }'''

    cf = parse_CF_JSON(json_string)

    assert len(cf.steps) == 1
    assert cf.steps[0].operation == "GoldenGate"
    assert cf.steps[0].inputs == ["Fragment1", "Fragment2"]
    assert cf.steps[0].enzyme == "BsaI"
    assert cf.steps[0].output == "ProductName"

def test_gibson():
    json_string = '''{
        "steps": [
            {
                "operation": "Gibson",
                "inputs": ["Fragment1", "Fragment2"],
                "output": "ProductName"
            }
        ],
        "sequences": {}
    }'''

    cf = parse_CF_JSON(json_string)

    assert len(cf.steps) == 1
    assert cf.steps[0].operation == "Gibson"
    assert cf.steps[0].inputs == ["Fragment1", "Fragment2"]
    assert cf.steps[0].output == "ProductName"

def test_transform():
    json_string = '''{
        "steps": [
            {
                "operation": "Transform",
                "dna": "ProductName",
                "strain": "E. coli",
                "antibiotics": ["Amp"],
                "output": "TransformedProduct",
                "temperature": 37
            }
        ],
        "sequences": {}
    }'''

    cf = parse_CF_JSON(json_string)

    assert len(cf.steps) == 1
    assert cf.steps[0].operation == "Transform"
    assert cf.steps[0].dna == "ProductName"
    assert cf.steps[0].strain == "E. coli"
    assert cf.steps[0].antibiotics == ["Amp"]
    assert cf.steps[0].output == "TransformedProduct"
    assert cf.steps[0].temperature == 37

def test_sequences():
    json_string = '''{
        "steps": [],
        "sequences": {
            "pTP1": {"sequence": "ATTACCGCCTTTGAGTGG", "ext5": "", "ext3": "", "is_double_stranded": true, "is_circular": true, "mod_ext5": null, "mod_ext3": null}
        }
    }'''

    cf = parse_CF_JSON(json_string)

    assert len(cf.sequences) == 1
    assert isinstance(cf.sequences["pTP1"], Polynucleotide)
    assert cf.sequences["pTP1"].sequence == "ATTACCGCCTTTGAGTGG"
    assert cf.sequences["pTP1"].is_double_stranded == True
    assert cf.sequences["pTP1"].is_circular == True

def test_unrecognized_operation():
    json_string = '''{
        "steps": [
            {
                "operation": "InvalidOperation",
                "output": "P6",
                "forward_oligo": "P6libF",
                "reverse_oligo": "P6libR",
                "template": "pTP1",
                "product_size": 3583
            }
        ],
        "sequences": {}
    }'''

    with pytest.raises(ValueError):
        cf = parse_CF_JSON(json_string)

def test_unrecognized_enzyme_golden_gate():
    json_string = '''{
        "steps": [
            {
                "operation": "GoldenGate",
                "inputs": ["Fragment1", "Fragment2"],
                "enzyme": "InvalidEnzyme",
                "output": "ProductName"
            }
        ],
        "sequences": {}
    }'''

    with pytest.raises(ValueError):
        cf = parse_CF_JSON(json_string)

def test_unrecognized_enzyme_digest():
    json_string = '''{
        "steps": [
            {
                "operation": "Digest",
                "dna": "pTP1",
                "enzymes": ["InvalidEnzyme"],
                "fragSelect": 3583,
                "output": "digested_dna"
            }
        ],
        "sequences": {}
    }'''

    with pytest.raises(ValueError):
        cf = parse_CF_JSON(json_string)

def test_multiple_steps():
    json_string = '''{
        "steps": [
            {
                "operation": "PCR",
                "output": "P6",
                "forward_oligo": "P6libF",
                "reverse_oligo": "P6libR",
                "template": "pTP1",
                "product_size": 3583
            },
            {
                "operation": "Digest",
                "dna": "pTP1",
                "enzymes": ["EcoRI"],
                "fragSelect": 3583,
                "output": "digested_dna"
            },
            {
                "operation": "Ligate",
                "dnas": ["Fragment1", "Fragment2"],
                "output": "LigatedProduct"
            }
        ],
        "sequences": {
            "pTP1": {"sequence": "ATTACCGCCTTTGAGTGG", "ext5": "", "ext3": "", "is_double_stranded": true, "is_circular": true, "mod_ext5": null, "mod_ext3": null}
        }
    }'''

    cf = parse_CF_JSON(json_string)

    assert len(cf.steps) == 3
    assert cf.steps[0].operation == "PCR"
    assert cf.steps[1].operation == "Digest"
    assert cf.steps[2].operation == "Ligate"

    assert len(cf.sequences) == 1
    assert isinstance(cf.sequences["pTP1"], Polynucleotide)
    assert cf.sequences["pTP1"].sequence == "ATTACCGCCTTTGAGTGG"

def test_sequences_and_steps():
    json_string = '''{
        "steps": [
            {
                "operation": "PCR",
                "output": "P6",
                "forward_oligo": "P6libF",
                "reverse_oligo": "P6libR",
                "template": "pTP1",
                "product_size": 3583
            }
        ],
        "sequences": {
            "pTP1": {"sequence": "ATTACCGCCTTTGAGTGG", "ext5": "", "ext3": "", "is_double_stranded": true, "is_circular": true, "mod_ext5": null, "mod_ext3": null}
        }
    }'''

    cf = parse_CF_JSON(json_string)

    assert len(cf.steps) == 1
    assert cf.steps[0].operation == "PCR"
    assert len(cf.sequences) == 1
    assert isinstance(cf.sequences["pTP1"], Polynucleotide)

def test_no_steps_or_sequences():
    json_string = '''{
        "steps": [],
        "sequences": {}
    }'''

    cf = parse_CF_JSON(json_string)

    assert len(cf.steps) == 0
    assert len(cf.sequences) == 0

def test_invalid_antibiotics():
    json_string = '''{
        "steps": [
            {
                "operation": "Transform",
                "dna": "ProductName",
                "strain": "E. coli",
                "antibiotics": ["InvalidAntibiotic"],
                "output": "TransformedProduct",
                "temperature": 37
            }
        ],
        "sequences": {}
    }'''

    with pytest.raises(ValueError):
        cf = parse_CF_JSON(json_string)

def test_missing_required_fields():
    json_string = '''{
        "steps": [
            {
                "operation": "PCR",
                "output": "P6",
                "forward_oligo": "P6libF",
                "reverse_oligo": "P6libR",
                "product_size": 3583
            }
        ],
        "sequences": {}
    }'''

    with pytest.raises(ValueError) as e:
        cf = parse_CF_JSON(json_string)
