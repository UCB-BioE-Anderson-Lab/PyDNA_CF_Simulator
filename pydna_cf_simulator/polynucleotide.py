class Polynucleotide:

    def __eq__(self, other):
        if not isinstance(other, Polynucleotide):
            return NotImplemented
        return (self.sequence == other.sequence and
                self.ext5 == other.ext5 and
                self.ext3 == other.ext3 and
                self.is_double_stranded == other.is_double_stranded and
                self.is_circular == other.is_circular and
                self.mod_ext5 == other.mod_ext5 and
                self.mod_ext3 == other.mod_ext3)
    def __init__(self, sequence, ext5, ext3, is_double_stranded, is_circular, mod_ext5, mod_ext3):
        self.sequence = sequence.upper() if sequence else sequence
        self.ext5 = ext5.upper() if ext5 else ext5
        self.ext3 = ext3.upper() if ext3 else ext3
        self.is_double_stranded = is_double_stranded
        self.is_circular = is_circular
        self.mod_ext5 = mod_ext5
        self.mod_ext3 = mod_ext3

    def __str__(self):
        return f'Polynucleotide(sequence={self.sequence}, ext5={self.ext5}, ext3={self.ext3}, is_double_stranded={self.is_double_stranded}, is_circular={self.is_circular}, mod_ext5={self.mod_ext5}, mod_ext3={self.mod_ext3})'

    def __repr__(self):
        return self.__str__()

def dsDNA(sequence):
    return Polynucleotide(sequence, '', '', True, False, 'hydroxyl', 'hydroxyl')


def oligo(sequence):
    return Polynucleotide(sequence, None, None, False, False, 'hydroxyl', None)


def plasmid(sequence):
    return Polynucleotide(sequence, '', '', True, True, None, None)