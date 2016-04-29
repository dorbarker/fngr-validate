#!/usr/bin/env python3

from Bio import SeqIO
from itertools import dropwhile
import argparse
import os
import random
import sys

def arguments():

    parser = argparse.ArgumentParser()

    return parser.parse_args()

def handle_input(filepath: str or sys.stdin) -> 'file_handle':

    o = sys.stdin if filepath == '-' else open(filepath, 'r')

    with o as f:
        return f

def load_genome(handle: 'file_handle') -> dict:
    """Parse FASTA formatted string"""

    g = {contig.id: str(contig.seq) for contig in SeqIO.parse(handle, 'fasta')}
    return g

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

def contigify(sequence: str, mean: float, stdev: float) -> dict:
    """Cut single-sequence genomes into artificial contigs based on
    empirical distribution contig sizes in draft assemblies
    """

    def chunk_genome():

        processed = 0
        counter = 0

        while processed < len(sequence):

            counter += 1

            n = int(random.gauss(mean, stdev))

            chunk = sequence[processed:processed + n]

            processed += n

            yield 'contig_{:04d}'.format(counter), chunk

    return dict(chunk_genome())

def select_subsequence(mean: float, stdev: float) -> str:
    """Return a gene-like subsequence from a source genome"""
    pass

def integrate(transposon: str, contig: str, breakpoint: int) -> str:
    """Take a gene-like subsequence and integrate it into a target contig"""

    first, last = contig[:breakpoint], contig[breakpoint:]
    return ''.join((first, transposon, last))

def contaminate(contaminant: str, genome: dict) -> dict:
    """Add a contamination contig to genome"""

    def suffix_max(dictionary):

        k = dictionary.keys()

        o = int(''.join(dropwhile(lambda x: x not in map(str, range(10)), k)))
        return o

    # not strictly necessary, but makes the program easier to follow
    contaminated_genome = dict(genome.items())

    name = 'contamination_{}'.format(suffix_max(contaminated_genome) + 1)
    contaminated_genome[name] = contaminant

    return contaminated_genome

def iterate():
    """Repeat a treatment with a random variable"""
    pass

def main():
    pass

if __name__ == '__main__':
    main()
