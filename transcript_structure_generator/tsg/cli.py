import argparse
import logging
from pathlib import Path

from tsg.main import sample_transcripts


def setup_logging(loglevel: str = None) -> None:
    """Set up logging. Loglevel can be one of ["CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"].

    Args:
        loglevel (str, optional): Level of log output. Defaults to None.

    Raises:
        ValueError: If string that is not a log level is passed, raise error.

    Returns:
        None
    """
    if loglevel:
        numeric_level = getattr(logging, loglevel.upper())
        if not isinstance(numeric_level, int):
            raise ValueError("Invalid log level: %s" % loglevel)
    else:
        numeric_level = getattr(logging, "INFO")

    logging.basicConfig(
        format='[%(asctime)s: %(levelname)s] %(message)s (module "%(module)s")',
        level=numeric_level,
    )


def build_arg_parser() -> argparse.ArgumentParser:
    """ Builds the argument parser.

    Args:
        1) path to the csv-file with the number of transcripts (str)
        2) path to the gtf-file with the annotations for each transcript (str)
        3) a value for the probability of intron inclusion (float)
        4) a log message (str)
    
    Raises:
        None  

    Returns:
        parser  
    """
    parser = argparse.ArgumentParser()
    parser.add_argument("--transcripts", type=str)
    parser.add_argument("--annotation", type=str)
    parser.add_argument("--prob_inclusion", type=float)
    parser.add_argument("--log", type=str)

    return parser


def get_args() -> argparse.Namespace:
    """Builds a parser and returns its arguments.

    Args:
        None
    
    Raises:
        None 

    Returns:
        arguments for parser   
    """
    parser = build_arg_parser()
    args = parser.parse_args()

    return args


def output_filename(filename: str) -> str:
    """Generate output filename for given input filename.

    Args:
        filename (str): Input filename

    Raises:
        NotImplementedError: Only accept filetypes .csv, .tsv and .gtf.
        FileExistsError: If the output file exists, raise error.

    Returns:
        str: Output filename
    """
    filepath = Path(filename)
    if filepath.suffix == ".csv" or filepath.suffix == ".tsv":
        outfile = "generated_" + filepath.stem + ".csv"
    elif filepath.suffix == ".gtf":
        outfile = "generated_" + filepath.name
    else:
        raise NotImplementedError()

    if Path(outfile).exists():
        raise FileExistsError(f"The output file {outfile} already exists.")
        
    return outfile


def app():
    """Gets the args, sets up the logging and starts the programm with the provided parameters.

    Args: 
        1) path to the csv-file with the number of transcripts (str)
        2) path to the gtf-file with the annotations for each transcript (str)
        3) a value for the probability of intron inclusion (float)
        4) a log message (str)
    
    Raises:
        None  

    Returns:
        None  
    """
    args = get_args()

    setup_logging(args.log)
    sample_transcripts(
        args.transcripts,
        args.annotation,
        args.prob_inclusion,
        output_filename(args.transcripts),
        output_filename(args.annotation),
    )