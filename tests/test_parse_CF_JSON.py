import pytest
import json

from pydna_cf_simulator.parse_CF_JSON import parse_CF_JSON
from pydna_cf_simulator.construction_file import ConstructionFile, PCR, Digest, Ligate, GoldenGate, Gibson, Transform
from pydna_cf_simulator.polynucleotide import Polynucleotide

def test_pcr():
    # JSON string for a Construction File with one PCR operation
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
            "P6libF": {"sequence": "ccaaaggtctcATTATANNNNNNNNNNNNNNNNNTGTCAANNNNGAacccaggactcctcgaagtcgttcttaagacaac", "ext5": null, "ext3": null, "is_double_stranded": false, "is_circular": false, "mod_ext5": "hydroxyl", "mod_ext3": null},
            "P6libR": {"sequence": "cagttGGTCTCAATAATNNNNNNANNNNGTtagtatttctcctcgtctacggttaactgatactc", "ext5": null, "ext3": null, "is_double_stranded": false, "is_circular": false, "mod_ext5": "hydroxyl", "mod_ext3": null},
            "pTP1": {"sequence": "ATTACCGCCTTTGAGTGG", "ext5": "", "ext3": "", "is_double_stranded": true, "is_circular": true, "mod_ext5": null, "mod_ext3": null}
        }
    }'''

    # Parse the Construction File
    cf = parse_CF_JSON(json_string)

    # Check the output of the PCR operation
    assert len(cf.steps) == 1
    assert cf.steps[0].operation == "PCR"
    assert cf.steps[0].output == "P6"
    assert cf.steps[0].forward_oligo == "P6libF"
    assert cf.steps[0].reverse_oligo == "P6libR"
    assert cf.steps[0].template == "pTP1"
    assert cf.steps[0].product_size == 3583
