import os
import logging
import pandas

from staramr.blast.AbstractBlastDatabase import AbstractBlastDatabase

logger = logging.getLogger('ResfinderBlastDatabase')

class ResfinderBlastDatabase(AbstractBlastDatabase):

    def __init__(self, database_dir):
        super().__init__(database_dir)
        self._parse_resfinder_info(os.path.join(database_dir, "notes.txt"))

    def _parse_resfinder_info(self, database_info_file):
        self._resfinder_info = pandas.read_csv(database_info_file, sep=':', comment="#", header=None, index_col=False,
                                               names=("gene", "phenotype"))

    def get_database_names(self):
        return [f[:-len(self.fasta_suffix)] for f in os.listdir(self.database_dir) if
                (os.path.isfile(os.path.join(self.database_dir, f)) and f.endswith(self.fasta_suffix))]

    def get_path(self, database_name):
        return os.path.join(self.database_dir, database_name + self.fasta_suffix)

    def get_phenotype(self, gene):
        table = self._resfinder_info
        phenotype = table[table['gene'] == gene]['phenotype']
        if phenotype.size == 0:
            logger.warning("No phenotype matches in resfinder info for gene [" + gene + "]")
            return 'None'
        elif phenotype.size == 1:
            return phenotype.iloc[0]
        else:
            raise Exception(
                "Invalid number of matches in resfinder info for gene [" + gene + "], got=" + str(phenotype))
