import logging

from staramr.databases.AMRDatabaseHandler import AMRDatabaseHandler
from staramr.databases.BlastDatabaseRepositoryStripGitDir import BlastDatabaseRepositoryStripGitDir

logger = logging.getLogger('AMRDatabaseHandlerStripGitDir')

"""
A Class used to handle interactions with the ResFinder/PointFinder database files, stripping out the .git directory.
"""


class AMRDatabaseHandlerStripGitDir(AMRDatabaseHandler):
    GIT_INFO_SECTION = 'GitInfo'

    def __init__(self, database_dir):
        """
        Creates a new AMRDatabaseHandlerStripGitDir.
        :param database_dir: The root directory for both the ResFinder/PointFinder databases.
        """
        super().__init__(database_dir)

        self._resfinder_repository = BlastDatabaseRepositoryStripGitDir(database_dir, 'resfinder',
                                                                        'https://bitbucket.org/genomicepidemiology/resfinder_db.git')
        self._pointfinder_repository = BlastDatabaseRepositoryStripGitDir(database_dir, 'pointfinder',
                                                                          'https://bitbucket.org/genomicepidemiology/pointfinder_db.git')

    def update(self, resfinder_commit=None, pointfinder_commit=None):
        """
        Updates an existing ResFinder/PointFinder database to the latest revisions (or passed specific revisions).
        :param resfinder_commit: The specific git commit for ResFinder.
        :param pointfinder_commit: The specific git commit for PointFinder.
        :return: None
        """
        raise Exception("Cannot update when .git directory has been removed")
