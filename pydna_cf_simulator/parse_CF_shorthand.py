import re

from .construction_file import ConstructionFile, PCR, Digest, Ligate, GoldenGate, Gibson, Transform
from .polynucleotide import Polynucleotide, oligo, plasmid, dsDNA

ALL_ENZYMES = ['AarI', 'BbsI', 'BsaI', 'BsmBI', 'SapI', 'BseRI', 'BamHI', 'BglII', 'EcoRI', 'XhoI', 'SpeI', 'XbaI', 'PstI', 'HindIII', 'NotI', 'XmaI', 'SmaI', 'KpnI', 'SacI', 'SalI']
TYPE_IIS_ENZYMES = ['AarI', 'BbsI', 'BsaI', 'BsmBI', 'SapI', 'BseRI']
VALID_ANTIBIOTICS = {'Amp', 'Carb', 'Cam', 'Kan', 'Gen', 'Spec', 'Trim'}

def parse_CF_shorthand(cf_shorthand):
    # Remove /*comment*/ style comments
    cf_shorthand = re.sub(r'/\*.*?\*/', '', cf_shorthand, flags=re.DOTALL)

    lines = cf_shorthand.split('\n')
    steps = []
    sequences = {}

    for line_num, line in enumerate(lines, start=1):
        # Remove # and // style comments
        line = re.split(r'#|//', line)[0]

        # Skip empty or whitespace-only lines
        if not line.strip():
            continue
        
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
                operation = elements[0].lower()

                if operation == 'pcr':
                    try:
                        if len(elements) == 6:
                            forward_primer, reverse_primer, template, product_size, product_name = elements[1:]
                            validate_int(product_size, line_num)
                            steps.append(PCR(forward_primer, reverse_primer, template, product_name, int(product_size)))
                        else:
                            forward_primer, reverse_primer, template, product_name = elements[1:]
                            steps.append(PCR(forward_primer, reverse_primer, template, product_name))
                    except ValueError:
                        raise ValueError(f"Error in line {line_num}: Invalid argument type for PCR operation.")
                    except IndexError:
                        raise ValueError(f"Error in line {line_num}: Invalid number of arguments for PCR operation.")
                elif operation == 'digest':
                    try:
                        if len(elements) == 5:
                            dna, enzymes, frag_select, product_name = elements[1:]
                            validate_int(frag_select, line_num)
                            frag_select = int(frag_select)
                            if frag_select < 0:
                                raise ValueError(f"Error in line {line_num}: Invalid sequence number '{frag_select}'. Must be a non-negative integer.")
                            enzymes = enzymes.split(',')
                            unrecognized_enzymes = [enzyme for enzyme in enzymes if enzyme not in ALL_ENZYMES]
                            if unrecognized_enzymes:
                                raise ValueError(f"Error in line {line_num}: Unrecognized enzyme(s): {', '.join(unrecognized_enzymes)}")
                            steps.append(Digest(dna, enzymes, frag_select, product_name))
                        else:
                            raise ValueError(f"Error in line {line_num}: Invalid number of arguments for Digest operation.")
                    except ValueError:
                        raise ValueError(f"Error in line {line_num}: Invalid argument type for Digest operation.")
                elif operation == 'ligate':
                    try:
                        if len(elements) >= 4:
                            dnas, product_name = elements[1:-1], elements[-1]
                            if len(dnas) < 1:
                                raise ValueError(f"Error in line {line_num}: Ligate operation requires at least 1 DNA inputs.")
                            steps.append(Ligate(dnas, product_name))
                        else:
                            raise ValueError(f"Error in line {line_num}: Invalid number of arguments for Ligate operation.")
                    except ValueError:
                        raise ValueError(f"Error in line {line_num}: Invalid argument type for Ligate operation.")
                elif operation == 'goldengate':
                    try:
                        if len(elements) >= 5:
                            inputs, enzyme, product_name = elements[1:-1], elements[-2], elements[-1]
                            if len(inputs) < 1:
                                raise ValueError(f"Error in line {line_num}: GoldenGate operation requires at least 1 inputs.")
                            if enzyme not in TYPE_IIS_ENZYMES:
                                raise ValueError(f"Invalid enzyme {enzyme} for GoldenGate. Must be one of {TYPE_IIS_ENZYMES}.")
                            steps.append(GoldenGate(inputs, enzyme, product_name))
                        else:
                            raise ValueError(f"Error in line {line_num}: Invalid number of arguments for GoldenGate operation.")
                    except ValueError:
                        raise ValueError(f"Error in line {line_num}: Invalid argument type for GoldenGate operation.")
                elif operation == 'gibson':
                    try:
                        if len(elements) >= 4:
                            inputs, product_name = elements[1:-1], elements[-1]
                            if len(inputs) < 1:
                                raise ValueError(f"Error in line {line_num}: Gibson operation requires at least 1 inputs.")
                            steps.append(Gibson(inputs, product_name))
                        else:
                            raise ValueError(f"Error in line {line_num}: Invalid number of arguments for Gibson operation.")
                    except ValueError:
                        raise ValueError(f"Error in line {line_num}: Invalid argument type for Gibson operation.")
                elif operation == 'transform':
                    try:
                        if len(elements) == 6:
                            dna, strain, antibiotics, temperature, output = elements[1:]
                            antibiotics = antibiotics.split(',')
                            validate_antibiotics(antibiotics, line_num)
                            validate_int(temperature, line_num)
                            steps.append(Transform(dna, strain, antibiotics, output, int(temperature)))
                        else:
                            dna, strain, antibiotics, output = elements[1:]
                            antibiotics = antibiotics.split(',')
                            validate_antibiotics(antibiotics, line_num)
                            steps.append(Transform(dna, strain, antibiotics, output))
                    except ValueError:
                        raise ValueError(f"Error in line {line_num}: Invalid argument type for Transform operation.")
                    except IndexError:
                        raise ValueError(f"Error in line {line_num}: Invalid number of arguments for Transform operation.")

                # Or it's a sequence with a biological type
                elif operation == 'oligo':
                    if len(elements) == 3:
                        # It's a sequence, store it
                        seqname, sequence = elements[1:]
                        if re.match("^[ATCGNRKYSWBVHDM]+$", sequence):
                            sequences[seqname] = oligo(sequence)
                        else:
                            raise ValueError(f"Error in line {line_num}: Invalid sequence format. Sequences must only contain characters 'A', 'T', 'C', 'G', 'N', 'R', 'K', 'Y', 'S', 'W', 'B', 'V', 'H', 'D', and be at least one character long.")
                    else:
                        raise ValueError(f"Error in line {line_num}: Invalid number of arguments for oligo operation.")
                elif operation == 'plasmid':
                    if len(elements) == 3:
                        # It's a sequence, store it
                        seqname, sequence = elements[1:]
                        if re.match("^[ATCGNRKYSWBVHDM]+$", sequence):
                            sequences[seqname] = plasmid(sequence)
                        else:
                            raise ValueError(f"Error in line {line_num}: Invalid sequence format. Sequences must only contain characters 'A', 'T', 'C', 'G', 'N', 'R', 'K', 'Y', 'S', 'W', 'B', 'V', 'H', 'D', and be at least one character long.")
                    else:
                        raise ValueError(f"Error in line {line_num}: Invalid number of arguments for plasmid operation.")
                elif operation == 'dsdna':
                    if len(elements) == 3:
                        # It's a sequence, store it
                        seqname, sequence = elements[1:]
                        if re.match("^[ATCGNRKYSWBVHDM]+$", sequence):
                            sequences[seqname] = dsDNA(sequence)
                        else:
                            raise ValueError(f"Error in line {line_num}: Invalid sequence format. Sequences must only contain characters 'A', 'T', 'C', 'G', 'N', 'R', 'K', 'Y', 'S', 'W', 'B', 'V', 'H', 'D', and be at least one character long.")
                    else:
                        raise ValueError(f"Error in line {line_num}: Invalid number of arguments for dsdna operation.")

                else:
                    raise ValueError(f"Error in line {line_num}: Invalid operation {operation}")
        except (ValueError, IndexError) as e:
            raise ValueError(f"Error in line {line_num}: {str(e)}")

    return ConstructionFile(steps, sequences)

def validate_int(value, line_num):
    try:
        int(value)
    except ValueError:
        raise ValueError(f"Error in line {line_num}: Invalid argument type. Must be an integer.")

def validate_antibiotics(antibiotics, line_num):
    unrecognized_antibiotics = [antibiotic for antibiotic in antibiotics if antibiotic not in VALID_ANTIBIOTICS]
    if unrecognized_antibiotics:
        raise ValueError(f"Error in line {line_num}: Unrecognized antibiotic(s): {', '.join(unrecognized_antibiotics)}")
