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
            enzymes = [Restriction.__dict__[name] for name in step.enzymes]

            fragselect = step.fragSelect
            # Simulate Digest and select fragment
            fragments = sequence.cut(enzymes)
            product = fragments[fragselect]
            # Store product
            dseqDictionary[product_name] = product
            product_poly = dseqrecord_to_polynucleotide(product, 'phosphate', 'phosphate')
            polyDictionary[product_name] = product_poly;
            
        elif operation == 'Ligate':
            # Get inputs
            fragments = [dseqDictionary[dna] for dna in step.dnas]
            # Simulate Ligate
            product = fragments[0]
            for fragment in fragments[1:]:
                product += fragment
            # Convert to Polynucleotide and check for circularization
            product_poly = dseqrecord_to_polynucleotide(product, polyDictionary[step.dnas[0]].mod_ext5, polyDictionary[step.dnas[-1]].mod_ext3)
            if product_poly.ext5 == product_poly.ext3 and product_poly.mod_ext5 == 'phosphate' and product_poly.mod_ext3 == 'phosphate':
                product_poly.sequence = product_poly.ext5 + product_poly.sequence
                product_poly.ext5 = ""
                product_poly.ext3 = ""
                product_poly.mod_ext5 = ""
                product_poly.mod_ext3 = ""
                product_poly.is_circular = True
            # Convert back to Dseqrecord for storage
            product = polynucleotide_to_dseqrecord(product_poly)
            # Store product
            dseqDictionary[product_name] = product
            polyDictionary[product_name] = product_poly;

        elif operation == 'Gibson':
            # Get inputs
            fragments = [dseqDictionary[dna] for dna in step.dnas]
            # Simulate Gibson Assembly
            assembly = Assembly(fragments)
            assemblies = assembly.assemble_circular()
            product = assemblies[0]
            # Store product
            dseqDictionary[product_name] = product
            product_poly = dseqrecord_to_polynucleotide(product, polyDictionary[step.dnas[0]].mod_ext5, polyDictionary[step.dnas[-1]].mod_ext3)
            polyDictionary[product_name] = product_poly;

        elif operation == 'GoldenGate':
            raise NotImplementedError('GoldenGate operation is not implemented')
 
        elif operation == 'Transform':
            # Get the input DNA
            transPoly = polyDictionary[step.dna]

            # Check if DNA is circular
            if not transPoly.is_circular:
                raise ValueError(f"The DNA in the Transform operation must be circular, but '{step.dna}' is not.")
            polyDictionary[product_name] = transPoly;

    return polyDictionary
