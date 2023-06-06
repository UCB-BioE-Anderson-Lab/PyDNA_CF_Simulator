import json

from .construction_file import ConstructionFile, PCR, Digest, Ligate, GoldenGate, Gibson, Transform
from .polynucleotide import Polynucleotide

ALL_ENZYMES = ['AarI', 'BbsI', 'BsaI', 'BsmBI', 'SapI', 'BseRI', 'BamHI', 'BglII', 'EcoRI', 'XhoI', 'SpeI', 'XbaI', 'PstI', 'HindIII', 'NotI', 'XmaI', 'SmaI', 'KpnI', 'SacI', 'SalI']
TYPE_IIS_ENZYMES = ['AarI', 'BbsI', 'BsaI', 'BsmBI', 'SapI', 'BseRI']
VALID_ANTIBIOTICS = {'Amp', 'Carb', 'Cam', 'Kan', 'Gen', 'Spec', 'Trim'}

def parse_CF_JSON(json_string):
    cf_dict = json.loads(json_string)
    
    # Parse sequences
    sequences = {}
    for name, seq_dict in cf_dict['sequences'].items():
        sequences[name] = Polynucleotide(seq_dict["sequence"], seq_dict["ext5"], seq_dict["ext3"], 
                                         seq_dict["is_double_stranded"], seq_dict["is_circular"], 
                                         seq_dict["mod_ext5"], seq_dict["mod_ext3"])
    
    # Parse steps
    steps = []
    for step_dict in cf_dict['steps']:
        op = step_dict.pop('operation')
        if op == 'PCR':
            steps.append(PCR(**step_dict))
        elif op == 'Digest':
            # Validate enzymes
            enzymes = step_dict.get('enzymes', [])
            unrecognized_enzymes = [enzyme for enzyme in enzymes if enzyme not in ALL_ENZYMES]
            if unrecognized_enzymes:
                raise ValueError(f"Unrecognized enzyme(s): {', '.join(unrecognized_enzymes)}")
            steps.append(Digest(**step_dict))
        elif op == 'Ligate':
            steps.append(Ligate(**step_dict))
        elif op == 'GoldenGate':
            # Validate enzyme
            enzyme = step_dict.get('enzyme')
            if enzyme not in TYPE_IIS_ENZYMES:
                raise ValueError(f"Invalid enzyme {enzyme} for GoldenGate. Must be one of {TYPE_IIS_ENZYMES}.")
            steps.append(GoldenGate(**step_dict))
        elif op == 'Gibson':
            steps.append(Gibson(**step_dict))
        elif op == 'Transform':
            # Validate antibiotics
            antibiotics = step_dict.get('antibiotics', [])
            unrecognized_antibiotics = [antibiotic for antibiotic in antibiotics if antibiotic not in VALID_ANTIBIOTICS]
            if unrecognized_antibiotics:
                raise ValueError(f"Unrecognized antibiotic(s): {', '.join(unrecognized_antibiotics)}")
            steps.append(Transform(**step_dict))
        else:
            raise ValueError(f"Unrecognized operation: {op}")
    
    return ConstructionFile(steps, sequences)

