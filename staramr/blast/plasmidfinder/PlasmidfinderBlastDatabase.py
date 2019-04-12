import logging
from os import path
from typing import List

import pandas as pd

from staramr.blast.AbstractBlastDatabase import AbstractBlastDatabase

logger = logging.getLogger('PlasmidfinderBlastDatabase')

"""
A Class for pulling information from the PlasmidFinder database.
"""


class PlasmidfinderBlastDatabase(AbstractBlastDatabase):

    def __init__(self, database_dir: str, database_type: str = None) -> None:
        """
        Creates a new PlasmidfinderBlastDatabase.
        :param database_dir: The specific PlasmidFinder database directory.
        """

        super().__init__(database_dir)
        self.database_type = database_type

        if database_type is not None:
            self.plasmidfinder_database_file = path.join(self.database_dir, database_type)
            self.plasmidfinder_database_file += self.fasta_suffix

            if (not path.isfile(self.plasmidfinder_database_file)):
                raise Exception(
                    "Error, plasmidfinder database type [{}] is either incorrect or plasmidfinder database not installed properly",
                    database_type)
            elif database_type not in PlasmidfinderBlastDatabase.get_database_types(database_dir):
                raise Exception("Plasmidfinder database type [{}] is not valid", database_type)

    def get_name(self) -> str:
        return 'plasmidfinder'

    def get_database_names(self) -> List[str]:
        if self.database_type is not None:
            names = [self.database_type]
        else:
            names = super().get_database_names()

        return names

    @classmethod
    def get_available_databases(self) -> list:
        """
        A Class Method to get a list of plasmidfinder databases that are currently supported by staramr.
        :return: The list of database_type currently supported by staramr.
        """
        return ['gram_positive', 'enterobacteriaceae']

    @classmethod
    def get_database_types(self, database_dir: str) -> list:
        """
        A Class Method to get the list of databases from the PlasmidFinder database root directory.
        :param database_dir: The PlasmidFinder database root directory.
        :return: A list of databases in Plasmidfinder.
        """
        config = pd.read_csv(path.join(database_dir, 'config'), sep='\t', comment='#', header=None,
                             names=['db_prefix', 'name', 'description'])

        return config['db_prefix'].tolist()
