class ConstructionFile:
    def __init__(self, steps, sequences):
        self.steps = steps
        self.sequences = sequences

    def __repr__(self):
        return f'ConstructionFile(steps={self.steps}, sequences={self.sequences})'


class PCR:
    def __init__(self, operation, output, forward_oligo, reverse_oligo, template, product_size):
        self.operation = operation
        self.output = output
        self.forward_oligo = forward_oligo
        self.reverse_oligo = reverse_oligo
        self.template = template
        self.product_size = product_size

    def __repr__(self):
        return f'PCR(operation={self.operation}, output={self.output}, forward_oligo={self.forward_oligo}, reverse_oligo={self.reverse_oligo}, template={self.template}, product_size={self.product_size})'


class GoldenGate:
    def __init__(self, operation, output, dnas, enzyme):
        self.operation = operation
        self.output = output
        self.dnas = dnas
        self.enzyme = enzyme

    def __repr__(self):
        return f'GoldenGate(operation={self.operation}, output={self.output}, dnas={self.dnas}, enzyme={self.enzyme})'


class Transform:
    def __init__(self, operation, output, dna, strain, antibiotics, temperature=None):
        self.operation = operation
        self.output = output
        self.dna = dna
        self.strain = strain
        self.antibiotics = antibiotics
        self.temperature = temperature

    def __repr__(self):
        return f'Transform(operation={self.operation}, output={self.output}, dna={self.dna}, strain={self.strain}, antibiotics={self.antibiotics}, temperature={self.temperature})'


class Digest:
    def __init__(self, operation, output, dna, fragSelect, enzymes):
        self.operation = operation
        self.output = output
        self.dna = dna
        self.fragSelect = fragSelect
        self.enzymes = enzymes

    def __repr__(self):
        return f'Digest(operation={self.operation}, output={self.output}, dna={self.dna}, fragSelect={self.fragSelect}, enzymes={self.enzymes})'


class Ligate:
    def __init__(self, operation, output, dnas):
        self.operation = operation
        self.output = output
        self.dnas = dnas

    def __repr__(self):
        return f'Ligate(operation={self.operation}, output={self.output}, dnas={self.dnas})'


class Gibson:
    def __init__(self, operation, output, dnas):
        self.operation = operation
        self.output = output
        self.dnas = dnas

    def __repr__(self):
        return f'Gibson(operation={self.operation}, output={self.output}, dnas={self.dnas})'