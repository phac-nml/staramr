import logging
import os

from staramr.blast.AbstractBlastDatabase import AbstractBlastDatabase

logger = logging.getLogger('ResfinderBlastDatabase')

"""
A Class for pulling information from the ResFinder database.
"""


class ResfinderBlastDatabase(AbstractBlastDatabase):

    def __init__(self, database_dir):
        """
        Creates a new ResfinderBlastDatabase.
        :param database_dir: The specific ResFinder database (drug class) directory.
        """
        super().__init__(database_dir)

    def get_database_names(self):
        return [f[:-len(self.fasta_suffix)] for f in os.listdir(self.database_dir) if
                (os.path.isfile(os.path.join(self.database_dir, f)) and f.endswith(self.fasta_suffix))]

    def get_path(self, database_name):
        return os.path.join(self.database_dir, database_name + self.fasta_suffix)
