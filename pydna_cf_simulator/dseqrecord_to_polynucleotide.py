from pydna.dseqrecord import Dseqrecord
from pydna_cf_simulator.polynucleotide import Polynucleotide
from Bio.Seq import Seq


def dseqrecord_to_polynucleotide(dseqrecord: Dseqrecord, mod_ext5: str, mod_ext3: str) -> Polynucleotide:
    """
    Convert a Dseqrecord to a Polynucleotide.
    """
    # Get sequence and extensions
    watson = dseqrecord.seq.watson
    crick = dseqrecord.seq.crick

    # Determine the type of overhangs
    lefty = 'prime3' if dseqrecord.seq.ovhg < 0 else 'blunt' if dseqrecord.seq.ovhg == 0 else 'prime5'
    lefty_overhang = dseqrecord.seq.ovhg
    print(f'lefty: {lefty}, lefty_overhang: {lefty_overhang}')
    
    righty_overhang = -1 * (len(watson) - len(crick) - lefty_overhang)
    righty = 'prime3' if righty_overhang < 0 else 'blunt' if righty_overhang == 0 else 'prime5'
    print(f'righty: {righty}, righty_overhang: {righty_overhang}')
    
    # Compute ext5 and ext3
    if lefty == 'blunt':
        ext_5 = ''
    elif lefty == 'prime5':
        ext_5 = watson[:lefty_overhang]
    else:  # lefty == 'prime3'
        ext_5 = str(Seq(crick[-abs(lefty_overhang):]).reverse_complement())
    
    if righty == 'blunt':
        ext_3 = ''
    elif righty == 'prime5':
        ext_3 = str(Seq(crick[:abs(righty_overhang)]).reverse_complement())
    else:  # righty == 'prime3'
        ext_3 = watson[-abs(righty_overhang):]
    
    print(f'Watson: {watson}, Calculated ext_3: {ext_3}')
    
    # Add hyphen if overhang is negative
    ext_5 = '-' + ext_5 if lefty_overhang < 0 else ext_5
    ext_3 = '-' + ext_3 if righty_overhang < 0 else ext_3
    
    print(f'ext_5: {ext_5}, ext_3: {ext_3}')
    
    # Compute sequence
    sequence = watson
    if lefty_overhang > 0:
        sequence = sequence[lefty_overhang:]
    if righty_overhang < 0:
        sequence = sequence[:len(sequence) - abs(righty_overhang)]
    
    print(f'sequence: {sequence}')
    
    # Create and return Polynucleotide
    return Polynucleotide(sequence, ext_5, ext_3, True, dseqrecord.circular, mod_ext5, mod_ext3)
