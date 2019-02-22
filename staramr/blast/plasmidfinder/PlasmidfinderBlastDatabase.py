import logging
import os

from staramr.blast.AbstractBlastDatabase import AbstractBlastDatabase

logger = logging.getLogger('PlasmidfinderBlastDatabase')

"""
A Class for pulling information from the PlasmidFinder database.
"""


class PlasmidfinderBlastDatabase(AbstractBlastDatabase):

    def __init__(self, database_dir: str) -> None:
        """
        Creates a new PlasmidfinderBlastDatabase.
        :param database_dir: The specific ResFinder database (drug class) directory.
        """
        super().__init__(database_dir)

    def get_name(self):
        return 'plasmidfinder'
