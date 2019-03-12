import logging
from os import path, listdir

import pandas as pd

from staramr.blast.AbstractBlastDatabase import AbstractBlastDatabase

logger = logging.getLogger('PlasmidfinderBlastDatabase')

"""
A Class for pulling information from the PlasmidFinder database.
"""


class PlasmidfinderBlastDatabase(AbstractBlastDatabase):

    def __init__(self, database_dir: str, bacteria=None) -> None:
        """
        Creates a new PlasmidfinderBlastDatabase.
        :param database_dir: The specific PlasmidFinder database directory.
        """

        super().__init__(database_dir)
        self.plasmidfinder_database_dir = path.join(self.database_dir)
        self.bacteria = bacteria

        if bacteria is not None:
            self.plasmidfinder_database_file = path.join(self.database_dir, bacteria)
            self.plasmidfinder_database_file += ".fsa"

            if (not path.isfile(self.plasmidfinder_database_file)):
                raise Exception(
                    "Error, plasmidfinder bacteria [" + bacteria + "] is either incorrect or plasmidfinder database not installed properly")
            elif bacteria not in PlasmidfinderBlastDatabase.get_bacterias(database_dir):
                raise Exception("Plasmidfinder bacteria [" + bacteria + "] is not valid")

    def get_name(self) -> str:
        return 'plasmidfinder'

    def get_database_names(self):
        if self.bacteria is not None:
            names = [self.bacteria]
        else:
            names = [f[:-len(self.fasta_suffix)] for f in listdir(self.plasmidfinder_database_dir) if
                (path.isfile(path.join(self.plasmidfinder_database_dir, f)) and f.endswith(self.fasta_suffix))]

        return names

    @classmethod
    def get_available_bacteria(cls):
        """
        A Class Method to get a list of bacteria that are currently supported by staramr.
        :return: The list of bacteria currently supported by staramr.
        """
        return ['gram_positive', 'enterobacteriaceae']

    @classmethod
    def get_bacterias(cls, database_dir):
        """
        A Class Method to get the list of bacteria from the PlasmidFinder database root directory.
        :param database_dir: The PlasmidFinder database root directory.
        :return: A list of bacteria.
        """
        config = pd.read_csv(path.join(database_dir, 'config'), sep='\t', comment='#', header=None,
                               names=['db_prefix', 'name', 'description'])

        return config['db_prefix'].tolist()