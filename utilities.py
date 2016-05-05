#!/usr/bin/env python3

import subprocess
from multiprocessing import cpu_count

def fasta(f):
    if '.f' not in f:
        msg = 'Requires FASTA files (*.fasta, *.f, *.fna, etc)'
        raise argparse.ArgumentTypeError(msg)
    return f

class Fngr(object):

    def __init__(self, fngr, organism, kraken_db, nt_db,
                 threshold = 100, fragment = 250, cores = None):

        self.fngr = fngr
        self.organism = organism
        self.kraken_db = kraken_db
        self.nt_db = blast_db
        self.fragment = fragment
        self.threshold = threshold
        self.cores = cores or cpu_count()

    def fngr(self, assembly):

        assem = self._format_assembly(assembly)

        cmd = ('python3', self.fngr,
               '--organism', self.organism,
               '--kraken-database', self.kraken_db,
               '--nt-database', self.blast_db,
               '--threshold', self.threshold,
               '--fragment', self.fragment,
               '--cores', self.cores,
               '-')

        out = subprocess.check_output(cmd, input = assem,
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
