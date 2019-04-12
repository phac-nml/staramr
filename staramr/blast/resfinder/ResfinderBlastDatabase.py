import logging

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

    def get_name(self):
        return 'resfinder'
