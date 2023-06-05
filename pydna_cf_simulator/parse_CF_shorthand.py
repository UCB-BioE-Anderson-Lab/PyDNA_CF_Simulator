from .construction_file import ConstructionFile, PCR, Digest, Ligate, GoldenGate, Gibson, Transform
from .polynucleotide import Polynucleotide, oligo, plasmid

ALL_ENZYMES = ['AarI', 'BbsI', 'BsaI', 'BsmBI', 'SapI', 'BseRI', 'BamHI', 'BglII', 'EcoRI', 'XhoI', 'SpeI', 'XbaI', 'PstI', 'HindIII', 'NotI', 'XmaI', 'SmaI', 'KpnI', 'SacI', 'SalI']
TYPE_IIS_ENZYMES = ['AarI', 'BbsI', 'BsaI', 'BsmBI', 'SapI', 'BseRI']

def parse_CF_shorthand(cf_shorthand):
    lines = cf_shorthand.split('\n')
    steps = []
    sequences = {}

    for line in lines:
        elements = line.split()

        # Check if the line is a sequence or a step
        if len(elements) == 2:
            # It's a sequence, store it
            name, sequence = elements
            if len(sequence) < 100:
                sequences[name] = oligo(sequence)
            else:
                sequences[name] = plasmid(sequence)
        else:
            # It's a step, parse it
            operation = elements[0]

            if operation == 'PCR':
                if len(elements) == 6:
                    forward_primer, reverse_primer, template, product_size, product_name = elements[1:]
                    steps.append(PCR(forward_primer, reverse_primer, template, product_name, int(product_size)))
                else:
                    forward_primer, reverse_primer, template, product_name = elements[1:]
                    steps.append(PCR(forward_primer, reverse_primer, template, product_name))
            elif operation == 'Digest':
                dna, enzymes, frag_select, product_name = elements[1:]
                enzymes = [enzyme for enzyme in enzymes.split(',') if enzyme in ALL_ENZYMES]
                steps.append(Digest(dna, enzymes, int(frag_select), product_name))
            elif operation == 'Ligate':
                dnas, product_name = elements[1:-1], elements[-1]
                steps.append(Ligate(dnas, product_name))
            elif operation == 'GoldenGate':
                inputs, enzyme, product_name = elements[1:-1], elements[-2], elements[-1]
                if enzyme not in TYPE_IIS_ENZYMES:
                    raise ValueError(f"Invalid enzyme {enzyme} for GoldenGate. Must be one of {TYPE_IIS_ENZYMES}.")
                steps.append(GoldenGate(inputs, enzyme, product_name))
            elif operation == 'Gibson':
                inputs, product_name = elements[1:-1], elements[-1]
                steps.append(Gibson(inputs, product_name))
            elif operation == 'Transform':
                if len(elements) == 6:
                    dna, strain, antibiotics, temperature, output = elements[1:]
                    antibiotics = antibiotics.split(',')
                    temperature = int(temperature)
                else:
                    dna, strain, antibiotics, output = elements[1:]
                    antibiotics = antibiotics.split(',')
                    temperature = None
                steps.append(Transform(dna, strain, antibiotics, output, temperature))

    return ConstructionFile(steps, sequences)
