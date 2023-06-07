from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer
from pydna.dseq import Dseq
from pydna.assembly import Assembly
from pydna.amplify import pcr
from Bio.Restriction import *
from .polynucleotide_to_dseqrecord import polynucleotide_to_dseqrecord
from .dseqrecord_to_polynucleotide import dseqrecord_to_polynucleotide

from .polynucleotide import oligo


def simulate_CF(construction_file):
    # Convert all Polynucleotide sequences in the CF to Dseqrecord or Primer
    polyDictionary = construction_file.sequences
    dseqDictionary = {}
    for name, poly in construction_file.sequences.items():
        if poly.is_double_stranded:
            dseqDictionary[name] = polynucleotide_to_dseqrecord(poly)
        else:
            dseqDictionary[name] = Primer(poly.sequence)

    # Iterate through the steps
    for step in construction_file.steps:
        operation = step.operation
        product_name = step.output

        # Switch based on the operation
        if operation == 'PCR':
            # Get inputs
            forward = dseqDictionary[step.forward_oligo]
            reverse = dseqDictionary[step.reverse_oligo]
            template = dseqDictionary[step.template]
            # Simulate PCR
            amplicon = pcr(forward, reverse, template)
            # Store product
            dseqDictionary[product_name] = amplicon
            product_poly = dseqrecord_to_polynucleotide(amplicon, polyDictionary[step.forward_oligo].mod_ext5, polyDictionary[step.reverse_oligo].mod_ext5) 
            polyDictionary[product_name] = product_poly;
            
        elif operation == 'Digest':
            # Get inputs
            sequence = dseqDictionary[step.dna]
            enzymes = [getattr(Restriction, name) for name in step.enzymes.split(',')]
            fragselect = step.fragSelect
            # Simulate Digest
            fragments = sequence.cut(enzymes)
            for i, fragment in enumerate(fragments):
                print(f"Fragment {i}:")
                print("Watson: ", fragment.seq.watson)
                print("Crick: ", fragment.seq.crick)
                print("Overhang: ", fragment.seq.ovhg)
                print("\n")

            # Adjust fragment selection if the sequence was circular
            if sequence.circular:
                if not sequence.seq.watson.startswith(fragments[0].seq.watson):
                    fragselect = (fragselect + 1) % len(fragments)

            product = fragments[fragselect]
            # Store product
            dseqDictionary[product_name] = product
            product_poly = dseqrecord_to_polynucleotide(product, 'phosphate', 'phosphate')
            polyDictionary[product_name] = product_poly;
            
        elif operation == 'Ligate':
            raise NotImplementedError('Ligate operation is not implemented')
        elif operation == 'GoldenGate':
            raise NotImplementedError('GoldenGate operation is not implemented')
        elif operation == 'Gibson':
            raise NotImplementedError('Gibson operation is not implemented')
        elif operation == 'Transform':
            raise NotImplementedError('Transform operation is not implemented')

    return polyDictionary
