import re
from construction_file import ConstructionFile, PCR, GoldenGate, Transform, Digest, Ligate, Gibson

def parse_CF_shorthand(CF_string):
    # Remove /*comment*/ style comments from CF_string
    CF_string = re.sub(r'/\*.*?\*/', '', CF_string, flags=re.DOTALL)

    # Split CF_string into lines
    CF_lines = CF_string.split('\n')

    # Store the sequences and steps
    sequences = {}
    steps = []

    # Parse the lines
    for line in CF_lines:
        # Remove # and // style comments
        line = re.split(r'#|//', line)[0]

        # Skip empty lines
        if not line.strip():
            continue

        # Split the line into elements
        elements = line.split()

        # Check if the line is a sequence or a step
        if len(elements) == 2:
            # It's a sequence, store it
            sequences[elements[0]] = elements[1]
        else:
            # It's a step, parse it
            # Get the operation and the product name
            operation = elements[0].lower()
            product_name = elements[-1]

            # Get the input sequences
            input_sequences = elements[1:-1]

            # Create the operation object

            if operation == 'pcr':
                if len(elements) != 5 and len(elements) != 6:
                    raise ValueError('Invalid number of parameters for PCR operation. Expected 5 or 6, got ' + str(len(elements)))
                product_size = None
                if len(elements) == 6:
                    if elements[4].isdigit():
                        product_size = int(elements[4])
                    else:
                        raise ValueError('Invalid product size for PCR operation. Expected an integer, got ' + elements[4])
                steps.append(PCR(operation, product_name, input_sequences[0], input_sequences[1], input_sequences[2], product_size))

            elif operation == 'digest':
                if len(elements) < 4:
                    raise ValueError('Invalid number of parameters for Digest operation. Expected at least 4, got ' + str(len(elements)))
                if not elements[3].isdigit():
                    raise ValueError('Invalid fragment selection for Digest operation. Expected an integer, got ' + elements[3])
                steps.append(Digest(operation, product_name, input_sequences[0], int(elements[3]), elements[2].split(',')))

            elif operation == 'ligate':
                if len(elements) < 3:
                    raise ValueError('Invalid number of parameters for Ligate operation. Expected at least 3, got ' + str(len(elements)))
                steps.append(Ligate(operation, product_name, input_sequences))

            elif operation == 'goldengate':
                if len(elements) != 4:
                    raise ValueError('Invalid number of parameters for GoldenGate operation. Expected 4, got ' + str(len(elements)))
                steps.append(GoldenGate(operation, product_name, input_sequences, elements[2]))

            elif operation == 'gibson':
                if len(elements) < 3:
                    raise ValueError('Invalid number of parameters for Gibson operation. Expected at least 3, got ' + str(len(elements)))
                steps.append(Gibson(operation, product_name, input_sequences))

            elif operation == 'transform':
                if len(elements) != 5 and len(elements) != 6:
                    raise ValueError('Invalid number of parameters for Transform operation. Expected 5 or 6, got ' + str(len(elements)))
                temperature = None
                if len(elements) == 6:
                    if elements[4].isdigit():
                        temperature = int(elements[4])
                    else:
                        raise ValueError('Invalid temperature for Transform operation. Expected an integer, got ' + elements[4])
                steps.append(Transform(operation, product_name, input_sequences[0], elements[2], elements[3], temperature))

            else:
                raise ValueError('Invalid operation: ' + operation)

    # Return the ConstructionFile object
    return ConstructionFile(steps, sequences)
