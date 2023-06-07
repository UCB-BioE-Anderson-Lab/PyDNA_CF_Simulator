from pydna.dseqrecord import Dseqrecord
from pydna.primer import Primer
from pydna.dseq import Dseq
from pydna.assembly import Assembly
from pydna.amplicon import Amplicon
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
            amplicon = Amplicon(dseqDictionary[step.forward_oligo], dseqDictionary[step.reverse_oligo], dseqDictionary[step.template])
            dseqDictionary[product_name] = amplicon
            product_poly = dseqrecord_to_polynucleotide(amplicon, polyDictionary[step.forward_oligo].mod_ext5, polyDictionary[step.reverse_oligo].mod_ext5) 
            polyDictionary[step.template] = product_poly;
        elif operation == 'Digest':
            raise NotImplementedError('Digest operation is not implemented')
        elif operation == 'Ligate':
            raise NotImplementedError('Ligate operation is not implemented')
        elif operation == 'GoldenGate':
            raise NotImplementedError('GoldenGate operation is not implemented')
        elif operation == 'Gibson':
            raise NotImplementedError('Gibson operation is not implemented')
        elif operation == 'Transform':
            raise NotImplementedError('Transform operation is not implemented')

    return polyDictionary
