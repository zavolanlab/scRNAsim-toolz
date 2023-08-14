# -*- coding: utf-8 -*-
"""
Changed on Fr Dec 23 16:49:50 2022

@author: RobinC
"""
import argparse
import logging
import main

class CLI():
    def create_parser(self):
        """This function creates the parser"""

        parser = argparse.ArgumentParser(
            prog = 'PrimingSitePredictor',
            description = 'Takes a cutoff energy and the predicts location of priming sites of transcripts',
            epilog = 'To predict or not to predict')
        parser.add_argument('--float', type=float, required=True, help='A energy-cutoff float number')
        parsed_args = parser.parse_args()
        energy_cutoff = parsed_args.float
        return energy_cutoff

    def letsgo(self):
        """This function creates a parser and prints the energycutoff"""
        energy_cutoff = CLI.create_parser()
        print(f"Your energy cutoff is {energy_cutoff}")

if __name__ == '__main__':
    logging.basicConfig(
        format='[%(asctime)s: %(levelname)s] %(message)s (module "%(module)s")',
        level=logging.INFO,
        )
    LOG = logging.getLogger(__name__)
    main()

