import logging
import subprocess
import time
from os import path

import git

from staramr.blast.pointfinder.PointfinderBlastDatabase import PointfinderBlastDatabase
from staramr.blast.resfinder.ResfinderBlastDatabase import ResfinderBlastDatabase

logger = logging.getLogger('AMRDatabaseHandler')

"""
A Class used to handle interactions with the ResFinder/PointFinder database files.
"""


class AMRDatabaseHandler:

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
        logger.info("Cloning resfinder db [" + self._resfinder_url + "] to [" + self._resfinder_dir + "]")
        resfinder_repo = git.repo.base.Repo.clone_from(self._resfinder_url, self._resfinder_dir)

        if resfinder_commit is not None:
            logger.info("Checking out resfinder commit " + resfinder_commit)
            resfinder_repo.git.checkout(resfinder_commit)

        logger.info("Cloning pointfinder db [" + self._pointfinder_url + "] to [" + self._pointfinder_dir + "]")
        pointfinder_repo = git.repo.base.Repo.clone_from(self._pointfinder_url, self._pointfinder_dir)

        if pointfinder_commit is not None:
            logger.info("Checking out pointfinder commit " + pointfinder_commit)
            pointfinder_repo.git.checkout(pointfinder_commit)

        self._blast_format()

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
            resfinder_repo = git.Repo(self._resfinder_dir)
            pointfinder_repo = git.Repo(self._pointfinder_dir)

            logger.info("Updating " + self._resfinder_dir)
            resfinder_repo.heads.master.checkout()
            resfinder_repo.remotes.origin.pull()

            if resfinder_commit is not None:
                logger.info("Checking out resfinder commit " + resfinder_commit)
                resfinder_repo.git.checkout(resfinder_commit)

            resfinder_repo.git.reset('--hard')

            logger.info("Updating " + self._pointfinder_dir)
            pointfinder_repo.heads.master.checkout()
            pointfinder_repo.remotes.origin.pull()

            if pointfinder_commit is not None:
                logger.info("Checking out pointfinder commit " + pointfinder_commit)
                pointfinder_repo.git.checkout(pointfinder_commit)

            resfinder_repo.git.reset('--hard')

            self._blast_format()

    def info(self):
        """
        Gets information on the ResFinder/PointFinder databases.
        :return: Database information as a list containing key/value pairs.
        """
        data = []

        resfinder_repo = git.Repo(self._resfinder_dir)
        resfinder_repo_head = resfinder_repo.commit('HEAD')

        data.append(['resfinder_db_dir', self._resfinder_dir])
        data.append(['resfinder_db_url', self._resfinder_url])
        data.append(['resfinder_db_commit', str(resfinder_repo_head)])
        data.append(
            ['resfinder_db_date', time.strftime("%a, %d %b %Y %H:%M", time.gmtime(resfinder_repo_head.committed_date))])

        pointfinder_repo = git.Repo(self._pointfinder_dir)
        pointfinder_repo_head = pointfinder_repo.commit('HEAD')
        data.append(['pointfinder_db_dir', self._pointfinder_dir])
        data.append(['pointfinder_db_url', self._pointfinder_url])
        data.append(['pointfinder_db_commit', str(pointfinder_repo_head)])
        data.append(['pointfinder_db_date',
                     time.strftime("%a, %d %b %Y %H:%M", time.gmtime(pointfinder_repo_head.committed_date))])

        return data

    def _blast_format(self):

        logger.info("Formatting resfinder db")
        resfinder_db = ResfinderBlastDatabase(self._resfinder_dir)
        for path in resfinder_db.get_database_paths():
            self._make_blast_db(path)

        logger.info("Formatting pointfinder db")
        for organism_db in PointfinderBlastDatabase.build_databases(self._pointfinder_dir):
            for path in organism_db.get_database_paths():
                self._make_blast_db(path)

    def _make_blast_db(self, path):
        command = ['makeblastdb', '-in', path, '-dbtype', 'nucl', '-parse_seqids']
        logger.debug(' '.join(command))
        subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE).check_returncode()

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
