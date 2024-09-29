import subprocess
import re
import numpy as np
import os
import webbrowser


def get_entropy(sequence: str, view_struct: bool):
    script_dir = os.path.dirname(__file__)  # Directory of the current script

    # Full path to RNAfold and RNAplot executables
    # For some reason, I cant import it so I just downloaded it
    rnafold_path = os.path.join(script_dir, '../ViennaRNA Package', 'RNAfold.exe')
    rnaplot_path = os.path.join(script_dir, '../ViennaRNA Package', 'RNAplot.exe')
    rna2d_path = os.path.join(script_dir, '../ViennaRNA Package', "RNA2Dfold.exe")

    # Run RNAfold to get the structure
    fold_result = subprocess.run(
        [rnafold_path, "-p"],
        input=sequence,
        capture_output=True,
        text=True
    )

    # Extract the structure line from the RNAfold output
    dot_bracket_pattern = r'^[^\s]+'
    structure = fold_result.stdout.splitlines()[1]
    structure = re.search(dot_bracket_pattern, structure)[0]

    # Save the sequence and structure to a file in dot-bracket format
    with open("output/rna_structure.txt", "w") as f:
        f.write(f"{sequence}\n{structure}\n")

    svg_out = os.path.join(script_dir, 'output')
    # Run RNAplot, creates an SVG file of the output
    result_plot = subprocess.run(
        [rnaplot_path, "-o", 'svg', "output/rna_structure.txt"],
        capture_output=True,
        text=True
    )

    if view_struct:
        # Create FORNA URL with parameters
        # Forna displays the 2D struct
        forna_url = f"http://nibiru.tbi.univie.ac.at/forna/forna.html?id=fasta&file=%3Eheader%5Cn{sequence}%5Cn{structure}"

        # Open the FORNA visualization in the web browser
        webbrowser.open(forna_url)

    result_2d = subprocess.run(
        [rna2d_path],
        input=sequence + "\n" + structure + "\n" + structure,
        capture_output=True,
        text=True
    )

    # Split lines and extract probabilities
    lines = result_2d.stdout.strip().split("\n")[5:]  # Skip headers
    # print(result_2d.stdout)
    probabilities = [float(line.split()[2]) for line in lines]

    # Calculate Shannon entropy
    def shannon_entropy(probabilities):
        return -np.sum([p * np.log2(p) for p in probabilities if p > 0])

    entropy_value = shannon_entropy(probabilities)
    print("Entropy:", entropy_value)


if __name__ == "__main__":
    # RNA sequence
    sequence = "AUGUAGAUUUUUGUGGGGGUUUUCAGUGGGGUAGUAGUGUUGAGUGGAGUGAG"
    get_entropy(sequence, view_struct=False)
