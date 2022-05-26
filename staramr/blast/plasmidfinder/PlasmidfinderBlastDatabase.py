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
        self.database_names = None

        existing_database_types = PlasmidfinderBlastDatabase.get_database_types(database_dir)

        if database_type is not None:
            # gram positive is special since at one point it corresponded to a single fasta file
            # but later it corresponds to multiple fasta files
            if database_type == 'gram_positive' and 'gram_positive' not in existing_database_types:
                config_df = self.get_config_table(database_dir)
                self.database_names = config_df.loc[config_df["name"] == "Gram Positive", "db_prefix"].tolist()
            else:
                if (database_type == 'enterobacteriaceae') and (database_type not in existing_database_types):
                    logger.info(f'PlasmidFinder database type [{database_type}] is out of date, switching to [enterobacteriales]')
                    database_type = 'enterobacteriales'

                self.plasmidfinder_database_file = path.join(self.database_dir, database_type)
                self.plasmidfinder_database_file += self.fasta_suffix
    
                if (not path.isfile(self.plasmidfinder_database_file)):
                    raise Exception(
                        f"Error, could not find file [{self.plasmidfinder_database_file}]. "
                        f"plasmidfinder database type [{database_type}] is either incorrect or plasmidfinder database not installed properly")
                elif database_type not in existing_database_types:
                    raise Exception(f"Plasmidfinder database type [{database_type}] is not valid")
                else:
                    self.database_names = [database_type]

    def get_name(self) -> str:
        return 'plasmidfinder'

    def get_database_names(self) -> List[str]:
        if self.database_names is not None:
            names = self.database_names
        else:
            names = super().get_database_names()

        return names

    @classmethod
    def get_available_databases(cls) -> list:
        """
        A Class Method to get a list of plasmidfinder databases that are currently supported by staramr.
        :return: The list of database_type currently supported by staramr.
        """
        return ['gram_positive', 'enterobacteriaceae', 'enterobacteriales']

    @classmethod
    def get_config_table(cls, database_dir: str) -> pd.DataFrame:
        """
        A Class Method to get the config table from the plasmidfinder database root directory.
        :return: A DataFrame containing the config table.
        """
        config = pd.read_csv(path.join(database_dir, 'config'), sep='\t', comment='#', header=None,
                             names=['db_prefix', 'name', 'description'])
        return config

    @classmethod
    def get_database_types(cls, database_dir: str) -> list:
        """
        A Class Method to get the list of databases from the PlasmidFinder database root directory.
        :param database_dir: The PlasmidFinder database root directory.
        :return: A list of databases in Plasmidfinder.
        """
        config = cls.get_config_table(database_dir)

        return config['db_prefix'].tolist()
