import logging
import shutil
import time
from collections import OrderedDict
from os import path

import git

from staramr.exceptions.DatabaseErrorException import DatabaseErrorException
from staramr.exceptions.DatabaseNotFoundException import DatabaseNotFoundException

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
        self._database_dir = database_dir
        self._resfinder_dir = path.join(database_dir, 'resfinder')
        self._pointfinder_dir = path.join(database_dir, 'pointfinder')

        self._resfinder_url = "https://bitbucket.org/genomicepidemiology/resfinder_db.git"
        self._pointfinder_url = "https://bitbucket.org/genomicepidemiology/pointfinder_db.git"

    def build(self, resfinder_commit=None, pointfinder_commit=None):
        """
        Downloads and builds a new ResFinder/PointFinder database.
        :param resfinder_commit: The specific git commit for ResFinder.
        :param pointfinder_commit: The specific git commit for PointFinder.
        :return: None
        """

        try:
            logger.info("Cloning resfinder db [%s] to [%s]", self._resfinder_url, self._resfinder_dir)
            resfinder_repo = git.repo.base.Repo.clone_from(self._resfinder_url, self._resfinder_dir)

            if resfinder_commit is not None:
                logger.info("Checking out resfinder commit %s", resfinder_commit)
                resfinder_repo.git.checkout(resfinder_commit)

            logger.info("Cloning pointfinder db [%s] to [%s]", self._pointfinder_url, self._pointfinder_dir)
            pointfinder_repo = git.repo.base.Repo.clone_from(self._pointfinder_url, self._pointfinder_dir)

            if pointfinder_commit is not None:
                logger.info("Checking out pointfinder commit %s", pointfinder_commit)
                pointfinder_repo.git.checkout(pointfinder_commit)
        except Exception as e:
            raise DatabaseErrorException("Could not build database in [" + self._database_dir + "]") from e

    def update(self, resfinder_commit=None, pointfinder_commit=None):
        """
        Updates an existing ResFinder/PointFinder database to the latest revisions (or passed specific revisions).
        :param resfinder_commit: The specific git commit for ResFinder.
        :param pointfinder_commit: The specific git commit for PointFinder.
        :return: None
        """

        if not path.exists(self._database_dir):
            self.build(resfinder_commit=resfinder_commit, pointfinder_commit=pointfinder_commit)
        else:
            try:
                resfinder_repo = git.Repo(self._resfinder_dir)
                pointfinder_repo = git.Repo(self._pointfinder_dir)

                logger.info("Updating %s", self._resfinder_dir)
                resfinder_repo.heads.master.checkout()
                resfinder_repo.remotes.origin.pull()

                if resfinder_commit is not None:
                    logger.info("Checking out resfinder commit %s", resfinder_commit)
                    resfinder_repo.git.checkout(resfinder_commit)

                resfinder_repo.git.reset('--hard')

                logger.info("Updating %s", self._pointfinder_dir)
                pointfinder_repo.heads.master.checkout()
                pointfinder_repo.remotes.origin.pull()

                if pointfinder_commit is not None:
                    logger.info("Checking out pointfinder commit %s", pointfinder_commit)
                    pointfinder_repo.git.checkout(pointfinder_commit)

                resfinder_repo.git.reset('--hard')
            except Exception as e:
                raise DatabaseErrorException("Could not build database in [" + self._database_dir + "]") from e

    def remove(self):
        """
        Removes the databases stored in this directory.
        :return: None
        """
        shutil.rmtree(self._database_dir)

    def info(self):
        """
        Gets information on the ResFinder/PointFinder databases.
        :return: Database information as a OrderedDict of key/value pairs.
        """
        data = OrderedDict()

        try:
            resfinder_repo = git.Repo(self._resfinder_dir)
            resfinder_repo_head = resfinder_repo.commit('HEAD')

            data['resfinder_db_dir'] = self._resfinder_dir
            data['resfinder_db_url'] = self._resfinder_url
            data['resfinder_db_commit'] = str(resfinder_repo_head)
            data['resfinder_db_date'] = time.strftime(self.TIME_FORMAT, time.gmtime(resfinder_repo_head.committed_date))

            pointfinder_repo = git.Repo(self._pointfinder_dir)
            pointfinder_repo_head = pointfinder_repo.commit('HEAD')

            data['pointfinder_db_dir'] = self._pointfinder_dir
            data['pointfinder_db_url'] = self._pointfinder_url
            data['pointfinder_db_commit'] = str(pointfinder_repo_head)
            data['pointfinder_db_date'] = time.strftime(self.TIME_FORMAT,
                                                        time.gmtime(pointfinder_repo_head.committed_date))

        except git.exc.NoSuchPathError as e:
            raise DatabaseNotFoundException('Invalid database in [' + self._database_dir + ']') from e

        return data

    def get_database_dir(self):
        """
        Gets the root database dir.
        :return: The root database dir.
        """
        return self._database_dir

    def get_resfinder_dir(self):
        """
        Gets the ResFinder database directory.
        :return: The ResFinder database directory.
        """
        return self._resfinder_dir

    def get_pointfinder_dir(self):
        """
        Gets the PointFinder database directory.
        :return: The PointFinder database directory.
        """
        return self._pointfinder_dir
