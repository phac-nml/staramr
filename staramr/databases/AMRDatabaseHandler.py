import logging
import shutil

from staramr.databases.BlastDatabaseRepository import BlastDatabaseRepository

logger = logging.getLogger('AMRDatabaseHandler')

"""
A Class used to handle interactions with the ResFinder/PointFinder database files.
"""


class AMRDatabaseHandler:
    TIME_FORMAT = "%a, %d %b %Y %H:%M"

    def __init__(self, database_dir):
        """
        Creates a new AMRDatabaseHandler.
        :param database_dir: The root directory for both the ResFinder/PointFinder databases.
        """
        self._resfinder_repository = BlastDatabaseRepository(database_dir, 'resfinder',
                                                             'https://bitbucket.org/genomicepidemiology/resfinder_db.git')
        self._pointfinder_repository = BlastDatabaseRepository(database_dir, 'pointfinder',
                                                               'https://bitbucket.org/genomicepidemiology/pointfinder_db.git')

    def build(self, resfinder_commit=None, pointfinder_commit=None):
        """
        Downloads and builds a new ResFinder/PointFinder database.
        :param resfinder_commit: The specific git commit for ResFinder.
        :param pointfinder_commit: The specific git commit for PointFinder.
        :return: None
        """

        self._resfinder_repository.build(resfinder_commit)
        self._pointfinder_repository.build(pointfinder_commit)

    def update(self, resfinder_commit=None, pointfinder_commit=None):
        """
        Updates an existing ResFinder/PointFinder database to the latest revisions (or passed specific revisions).
        :param resfinder_commit: The specific git commit for ResFinder.
        :param pointfinder_commit: The specific git commit for PointFinder.
        :return: None
        """

        self._resfinder_repository.update(resfinder_commit)
        self._pointfinder_repository.update(pointfinder_commit)

    def remove(self):
        """
        Removes the databases stored in this directory.
        :return: None
        """
        self._resfinder_repository.remove()
        self._pointfinder_repository.remove()

    def info(self):
        """
        Gets information on the ResFinder/PointFinder databases.
        :return: Database information as a OrderedDict of key/value pairs.
        """
        info = self._resfinder_repository.info()
        info.update(self._pointfinder_repository.info())

        return info

    def get_database_dir(self):
        """
        Gets the root database dir.
        :return: The root database dir.
        """
        return self._resfinder_repository.get_database_dir()

    def get_resfinder_dir(self):
        """
        Gets the ResFinder database directory.
        :return: The ResFinder database directory.
        """
        return self._resfinder_repository.get_git_dir()

    def get_pointfinder_dir(self):
        """
        Gets the PointFinder database directory.
        :return: The PointFinder database directory.
        """
        return self._pointfinder_repository.get_git_dir()
