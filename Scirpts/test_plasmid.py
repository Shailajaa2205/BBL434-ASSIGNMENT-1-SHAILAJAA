import os
from plasmid_designer import read_fasta

def test_ecori_removed():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_path = os.path.join(base_dir, "Output", "Output.Fa")

    seq = read_fasta(output_path)

    assert "GAATTC" not in seq, "EcoRI (GAATTC) still present!"
    assert "CTTAAG" not in seq, "EcoRI reverse complement still present!"

    print("Test Passed: EcoRI successfully removed.")


def test_output_exists():
    base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    output_path = os.path.join(base_dir, "Output", "Output.Fa")

    assert os.path.exists(output_path), "Output file not created."

    print("Test Passed: Output file exists.")


if __name__ == "__main__":
    test_output_exists()
    test_ecori_removed()
