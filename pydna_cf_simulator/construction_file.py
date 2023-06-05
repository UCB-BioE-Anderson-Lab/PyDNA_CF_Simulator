class Step:
    def __init__(self, output):
        self.output = output


class PCR(Step):
    def __init__(self, forward_oligo, reverse_oligo, template, output, product_size=None):
        super().__init__(output)
        self.forward_oligo = forward_oligo
        self.reverse_oligo = reverse_oligo
        self.template = template
        self.product_size = product_size


class Digest(Step):
    def __init__(self, dna, enzymes, fragSelect, output, product_size=None):
        super().__init__(output)
        self.dna = dna
        self.enzymes = enzymes
        self.fragSelect = fragSelect
        self.product_size = product_size


class Ligate(Step):
    def __init__(self, dnas, output):
        super().__init__(output)
        self.dnas = dnas


class GoldenGate(Step):
    def __init__(self, inputs, enzyme, output):
        super().__init__(output)
        self.inputs = inputs
        self.enzyme = enzyme


class Gibson(Step):
    def __init__(self, inputs, output):
        super().__init__(output)
        self.inputs = inputs


class Transform(Step):
    def __init__(self, dna, strain, antibiotics, output, temperature=None):
        super().__init__(output)
        self.dna = dna
        self.strain = strain
        self.antibiotics = antibiotics
        self.temperature = temperature

class ConstructionFile:
    def __init__(self, steps, sequences):
        self.steps = steps
        self.sequences = sequences
