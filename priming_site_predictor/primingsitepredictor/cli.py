# -*- coding: utf-8 -*-
"""
Created on Mon Nov 14 14:49:50 2022

@author: baerma
"""
import argparse
import logging

def create_parser():
    """This function creates the parser"""
    parser = argparse.ArgumentParser(
    prog = 'Priming site predictor',
    description = 'Takes a cutoff energy and the predicts location of priming sites of transcripts',
    epilog = 'To predict or not to predict')
    parser.add_argument('energycutoff', type=float, help='a float as energy Cutoff')
    #parser.add_argument('transcripts', help='fastafile containing transcripts') #What type is that? fasta? Actually doesn't make sense here
    args = parser.parse_args()
    energy_cutoff = args.energycutoff
    return energy_cutoff
#possibly make a class out of this although I think it's an overkill

def letsgo():
    energy_cutoff = create_parser()
    print(f"Your energy cutoff is {energy_cutoff}")

if __name__ == '__main__':
    logging.basicConfig(
        format='[%(asctime)s: %(levelname)s] %(message)s (module "%(module)s")',
        level=logging.INFO,
    )
    LOG = logging.getLogger(__name__)
    letsgo()
    #here we would point to the main module and parse the energy cutoff
