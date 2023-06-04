from pydna.dseq import Dseq
from pydna.dseqrecord import Dseqrecord
from pydna.assembly import Assembly
from pydna.amplify import pcr
from pydna.primer import Primer
from Bio.Restriction import Restriction

from construction_file import PCR, Digest, Ligate, GoldenGate, Gibson, Transform


def simulate_CF(CF):
    # Store the result DNA sequences
    result_sequences = CF.sequences

    # Perform the operations
    for step in CF.steps:
        # Perform the operation
        if type(step) == PCR:
            # Get the input sequences
            forward_oligo = result_sequences[step.forward_oligo]
            reverse_oligo = result_sequences[step.reverse_oligo]
            template = result_sequences[step.template]

            # Perform the PCR operation
            product_sequence = pcr(forward_oligo, reverse_oligo, template)

        elif type(step) == Digest:
            # Get the input sequence
            dna = Dseq(result_sequences[step.dna])

            # Create the enzyme objects
            enzymes = [getattr(Restriction, enzyme) for enzyme in step.enzymes]

            # Digest the DNA with the enzymes
            product_sequences = dna.cut(*enzymes)
            # Get the index of the fragment to carry forward
            fragment_index = step.fragSelect
            # Check if the index is valid
            if fragment_index < 0 or fragment_index >= len(product_sequences):
                raise ValueError('Invalid fragment index: ' + str(fragment_index))
            # Carry forward the specified fragment
            product_sequence = product_sequences[fragment_index]


        elif type(step) == Ligate:
            # Assume the first two sequences are fragments
            # Ligate the fragments
            product_sequence = input_sequences[0] + input_sequences[1]

        elif type(step) == GoldenGate:
            # Assume the first two sequences are fragments and the third one is the enzyme
            # Digest the fragments with the enzyme and then ligate them
            product_sequence = input_sequences[0].cut(step.enzyme) + input_sequences[1].cut(step.enzyme)

        elif type(step) == Gibson:
            # Assume the first two sequences are fragments
            # Perform Gibson assembly
            product_sequence = Assembly(input_sequences).assemble_circular()[0]

        elif type(step) == Transform:
            # Not sure how to simulate transformation, just return the plasmid sequence
            product_sequence = input_sequences[0]

        else:
            raise ValueError('Invalid operation: ' + step.operation)

        # Check if the sequence is an oligo or a plasmid
        if len(product_sequence) < 100:
            # It's an oligo, make it single stranded and linear
            product_sequence = Dseqrecord(product_sequence.seq, linear=True)
        else:
            # It's a plasmid, make it double stranded and circular
            product_sequence = Dseqrecord(product_sequence.seq, circular=True)

        # Store the product sequence
        result_sequences[step.output] = product_sequence

    # Return the resulting sequences
    return result_sequences
