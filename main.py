import sys
from parse_CF_shorthand import parse_CF_shorthand
from simulate_CF import simulate_CF

# Get the CF shorthand string from the command line arguments
CF_shorthand_string = sys.argv[1]

# Parse the CF shorthand string
CF = parse_CF_shorthand(CF_shorthand_string)

# Simulate the CF
result_sequences = simulate_CF(CF)

# Print the result sequences
for name, sequence in result_sequences.items():
    print(f'{name}: {sequence}')