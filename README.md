# PyDNA Construction File Simulator

This is a ChatGPT plugin that wraps the PyDNA library to simulate construction files (CFs). The plugin is designed to be used with OpenAI's ChatGPT model, and it allows the model to simulate molecular biology experiments described in CF format.

## Construction Files (CFs)

A Construction File (CF) is a list of steps, each representing a specific operation in a molecular biology experiment. Each step is written on a new line. Parameters are separated by whitespace, preferably TSV. The step includes the names of input DNA sequence(s), non-sequence parameters, and ends with the name of the product DNA sequence. The input sequences can refer to products of previous steps. The product is the output of the operation.

In addition, a CF can include sequences in the form of 'name sequence' lines.

Comment lines follow '#', '//', or '/*comment*/' syntax.

This wrapper supports two versions of CFs: JSON and shorthand, but only shorthand is hosted as an API. The specifications for both versions can be found in the `docs` directory.  It can simulate PCR, restriction enzyme digestion, ligation, and Gibson assembly.

## Limitations

This plugin is an experiment exploring the ability of ChatGPT to design and simulate genetic engineering experiments. It is not intended for production use and has not been tested in real-world scenarios. The following limitations should be noted:

- It cannot handle 5' phosphates entirely.
- It cannot perform Golden Gate or Transform operations.
- Due to token limits, it is constrained to short sequences.

## Usage

To use this plugin, you need to have a running instance of the ChatGPT model with this plugin installed. The plugin provides a `simulate` function that takes a CF as input and returns the result of the simulation.

Here is an example of how to use the `simulate` function with CF shorthand:

```
# Define sequences
forward CCGCAACACACTTAACCTTG
reverse GTGGTTGTGGCCGGTCAAATC
template CCGCAACACACTTAACCTTGGCGTCGGGATACGTACATTGGAGAACGGTTGGCTGTACGGACTTAATACTTTTTATGATAATGATTTGACCGGCCACAACCACCG

# Define steps
PCR forward reverse template product
```

This will simulate a PCR operation using the defined sequences and return the result.  

## Installation

This project is not available via pip. To install it, you need to clone the repository and set it up manually. Here is the clone URL:

```
https://github.com/UCB-BioE-Anderson-Lab/PyDNA_CF_Simulator.git
```

After cloning the repository, navigate to the project directory and install the dependencies using pip:

```
pip install -r requirements.txt
```

Then, you can start the Flask server by running:

```
python app.py
```

The server will start on `http://localhost:8234`.