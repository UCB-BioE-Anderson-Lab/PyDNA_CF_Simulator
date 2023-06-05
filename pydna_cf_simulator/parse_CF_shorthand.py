import re

from .construction_file import ConstructionFile, PCR, Digest, Ligate, GoldenGate, Gibson, Transform
from .polynucleotide import Polynucleotide, oligo, plasmid

ALL_ENZYMES = ['AarI', 'BbsI', 'BsaI', 'BsmBI', 'SapI', 'BseRI', 'BamHI', 'BglII', 'EcoRI', 'XhoI', 'SpeI', 'XbaI', 'PstI', 'HindIII', 'NotI', 'XmaI', 'SmaI', 'KpnI', 'SacI', 'SalI']
TYPE_IIS_ENZYMES = ['AarI', 'BbsI', 'BsaI', 'BsmBI', 'SapI', 'BseRI']
VALID_ANTIBIOTICS = {'Amp', 'Carb', 'Cam', 'Kan', 'Gen', 'Spec', 'Trim'}

def parse_CF_shorthand(cf_shorthand):
    lines = cf_shorthand.split('\n')
    steps = []
    sequences = {}

    for line_num, line in enumerate(lines, start=1):
        elements = line.split()

        try:
            # Check if the line is a sequence or a step
            if len(elements) == 2:
                # It's a sequence, store it
                name, sequence = elements
                if re.match("^[ATCGNRKYSWBVHDM]+$", sequence):
                    if len(sequence) < 100:
                        sequences[name] = oligo(sequence)
                    else:
                        sequences[name] = plasmid(sequence)
                else:
                    raise ValueError(f"Error in line {line_num}: Invalid sequence format. Sequences must only contain characters 'A', 'T', 'C', 'G', 'N', 'R', 'K', 'Y', 'S', 'W', 'B', 'V', 'H', 'D', and be at least one character long.")
            else:
                # It's a step, parse it
                operation = elements[0]

                if operation == 'PCR':
                    try:
                        if len(elements) == 6:
                            forward_primer, reverse_primer, template, product_size, product_name = elements[1:]
                            steps.append(PCR(forward_primer, reverse_primer, template, product_name, int(product_size)))
                        else:
                            forward_primer, reverse_primer, template, product_name = elements[1:]
                            steps.append(PCR(forward_primer, reverse_primer, template, product_name))
                    except ValueError:
                        raise ValueError(f"Error in line {line_num}: Invalid number of arguments for PCR operation.")
                elif operation == 'Digest':
                    try:
                        if len(elements) == 5:
                            dna, enzymes, frag_select, product_name = elements[1:]
                            enzymes = enzymes.split(',')
                            unrecognized_enzymes = [enzyme for enzyme in enzymes if enzyme not in ALL_ENZYMES]
                            if unrecognized_enzymes:
                                raise ValueError(f"Error in line {line_num}: Unrecognized enzyme(s): {', '.join(unrecognized_enzymes)}")
                            steps.append(Digest(dna, enzymes, int(frag_select), product_name))
                        else:
                            raise ValueError(f"Error in line {line_num}: Invalid number of arguments for Digest operation.")
                    except ValueError:
                        raise ValueError(f"Error in line {line_num}: Invalid argument type for Digest operation.")
                elif operation == 'Ligate':
                    try:
                        if len(elements) >= 4:
                            dnas, product_name = elements[1:-1], elements[-1]
                            if len(dnas) < 2:
                                raise ValueError(f"Error in line {line_num}: Ligate operation requires at least 2 DNA inputs.")
                            steps.append(Ligate(dnas, product_name))
                        else:
                            raise ValueError(f"Error in line {line_num}: Invalid number of arguments for Ligate operation.")
                    except ValueError:
                        raise ValueError(f"Error in line {line_num}: Invalid argument type for Ligate operation.")
                elif operation == 'GoldenGate':
                    try:
                        if len(elements) >= 5:
                            inputs, enzyme, product_name = elements[1:-1], elements[-2], elements[-1]
                            if len(inputs) < 2:
                                raise ValueError(f"Error in line {line_num}: GoldenGate operation requires at least 2 inputs.")
                            if enzyme not in TYPE_IIS_ENZYMES:
                                raise ValueError(f"Invalid enzyme {enzyme} for GoldenGate. Must be one of {TYPE_IIS_ENZYMES}.")
                            steps.append(GoldenGate(inputs, enzyme, product_name))
                        else:
                            raise ValueError(f"Error in line {line_num}: Invalid number of arguments for GoldenGate operation.")
                    except ValueError:
                        raise ValueError(f"Error in line {line_num}: Invalid argument type for GoldenGate operation.")
                elif operation == 'Gibson':
                    try:
                        if len(elements) >= 4:
                            inputs, product_name = elements[1:-1], elements[-1]
                            if len(inputs) < 2:
                                raise ValueError(f"Error in line {line_num}: Gibson operation requires at least 2 inputs.")
                            steps.append(Gibson(inputs, product_name))
                        else:
                            raise ValueError(f"Error in line {line_num}: Invalid number of arguments for Gibson operation.")
                    except ValueError:
                        raise ValueError(f"Error in line {line_num}: Invalid argument type for Gibson operation.")
                elif operation == 'Transform':
                    try:
                        if len(elements) == 6:
                            dna, strain, antibiotics, temperature, output = elements[1:]
                            antibiotics = antibiotics.split(',')
                            validate_antibiotics(antibiotics, line_num)
                            try:
                                temperature = int(temperature)
                            except ValueError:
                                raise ValueError(f"Error in line {line_num}: Invalid temperature '{temperature}'. Must be an integer.")
                        else:
                            dna, strain, antibiotics, output = elements[1:]
                            antibiotics = antibiotics.split(',')
                            validate_antibiotics(antibiotics, line_num)
                            temperature = None
                        steps.append(Transform(dna, strain, antibiotics, output, temperature))
                    except ValueError:
                        raise ValueError(f"Error in line {line_num}: Invalid argument type for Transform operation.")
                    except IndexError:
                        raise ValueError(f"Error in line {line_num}: Invalid number of arguments for Transform operation.")
                else:
                    raise ValueError(f"Error in line {line_num}: Invalid operation {operation}")
        except (ValueError, IndexError) as e:
            raise ValueError(f"Error in line {line_num}: {str(e)}")

    return ConstructionFile(steps, sequences)

def validate_antibiotics(antibiotics, line_num):
    unrecognized_antibiotics = [antibiotic for antibiotic in antibiotics if antibiotic not in VALID_ANTIBIOTICS]
    if unrecognized_antibiotics:
        raise ValueError(f"Error in line {line_num}: Unrecognized antibiotic(s): {', '.join(unrecognized_antibiotics)}")
