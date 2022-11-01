import argparse
from pathlib import Path

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--transcripts", type=str)
    parser.add_argument("--annotation", type=str)
    parser.add_argument("--prob_inclusion", type=float)
    args = parser.parse_args()

    input_transcripts_file = args.transcripts
    input_annotations_file = args.annotation
    prob_inclusion = args.prob_inclusion
    input_transcripts_path = Path(input_transcripts_file)
    input_annotations_path = Path(input_annotations_file)
    output_transcripts_file = "generated_" + input_transcripts_path.stem + ".csv"
    output_annotations_file = "generated_" + input_annotations_path.name
