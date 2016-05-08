#!/usr/bin/env python3

import subprocess
import collections
from multiprocessing import cpu_count

Metadata = collections.namedtuple('Metadata', ['contig', 'start', 'length'])

Result = collections.namedtuple('Results', ['result', 'metadata'])

def fasta(f):
    if '.f' not in f:
        msg = 'Requires FASTA files (*.fasta, *.f, *.fna, etc)'
        raise argparse.ArgumentTypeError(msg)
    return f

class Fngr(object):

    def __init__(self, prog, organism, kraken_db, nt_db,
                 threshold = 100, fragment = 250, cores = None):

        self.prog = prog
        self.organism = organism
        self.kraken_db = kraken_db
        self.nt_db = nt_db
        self.fragment = fragment
        self.threshold = threshold
        self.cores = cores or cpu_count()

    def fngr(self, assembly):

        assem = self._format_assembly(assembly)

        cmd = ('python3', self.prog,
               '--organism', self.organism,
               '--kraken-database', self.kraken_db,
               '--nt-database' if self.nt_db else '', self.nt_db or '',
               '--threshold', str(self.threshold),
               '--fragment', str(self.fragment),
               '--cores', str(self.cores),
               '-')

        out = subprocess.check_output([x for x in cmd if x], input = assem,
                                      universal_newlines = True)

        return out  # JSON

    def _format_assembly(self, assembly):

        def dict_to_fasta(assem_dict):

            out = []
            for name, seq in assem_dict.items():

                out.extend([">{}".format(name), seq])

            return '\n'.join(out)

        if type(assembly) is dict:
            assem = dict_to_fasta(assembly)
        elif type(assembly) is str:
            assem = assembly
        else:
            raise TypeError('assembly is not of type str or dict')

        return assem
