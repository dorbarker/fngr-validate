#!/usr/bin/env python3

from Bio import SeqIO
import argparse
import os
import random
import sys

def arguments():

    parser = argparse.ArgumentParser()

    return parser.parse_args()

def load_genome(handle: str) -> dict:
    """Parse FASTA formatted string"""
    pass

def validate_ingroup(ingroup: list):
    """Use a set of high quality genomes of the target species
    to look for any spurious identification of 'foreign' sequence
    """
    pass

def validate_outgroup(ingroup: list):
    """Use a set of high quality genomes known not to be of the target
    species to validate that fngr correctly identifies foreign sequence
    """
    pass

def contigify(mean: float, stdev: float) -> dict:
    """Cut single-sequence genomes into artificial contigs based on 
    empirical distribution contig sizes in draft assemblies
    """
    pass

def select_subsequence(mean: float, stdev: float) -> str:
    """Return a gene-like subsequence from a source genome"""
    pass

def integrate(transposon: str, contig: str) -> str:
    """Take a gene-like subsequence and integrate it into a target contig"""
    pass

def iterate():
    """Repeat a treatment with a random variable"""
    pass

def main():
    pass

if __name__ == '__main__':
    main()
