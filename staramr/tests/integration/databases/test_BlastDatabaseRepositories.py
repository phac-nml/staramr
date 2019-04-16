import tempfile
import unittest
from os import path

import git

from staramr.databases.BlastDatabaseRepositories import BlastDatabaseRepositories


class BlastDatabaseRepositoriesIT(unittest.TestCase):
    RESFINDER_VALID_COMMIT = 'dc33e2f9ec2c420f99f77c5c33ae3faa79c999f2'
    RESFINDER_VALID_COMMIT2 = 'a4a699f3d13974477c7120b98fb0c63a1b70bd16'
    POINTFINDER_VALID_COMMIT = 'ba65c4d175decdc841a0bef9f9be1c1589c0070a'
    POINTFINDER_VALID_COMMIT2 = '0de22bff78214208171aef70461c639227e62e5d'
    PLASMIDFINDER_VALID_COMMIT = '21d154b1ccf877348e27da7782e25323760662d1'
    PLASMIDFINDER_VALID_COMMIT2 = 'a5c40dc39bc5f2ba7f5dbe2f8d4c341ce24f8dfe'

    def setUp(self):
        self.databases_dir = tempfile.TemporaryDirectory()
        self.database_repositories = BlastDatabaseRepositories.create_default_repositories(self.databases_dir.name)

    def tearDown(self):
        self.databases_dir.cleanup()

    def testBuild(self):
        # Verify that databases don't exist beforehand
        self.assertFalse(path.exists(self.database_repositories.get_repo_dir('resfinder')),
                         'resfinder path exists before creation of database')
        self.assertFalse(path.exists(self.database_repositories.get_repo_dir('pointfinder')),
                         'pointfinder path exists before creation of database')
        self.assertFalse(path.exists(self.database_repositories.get_repo_dir('plasmidfinder')),
                         'plasmidfinder path exists before creation of database')

        # Build database
        self.database_repositories.build(
            {'resfinder': self.RESFINDER_VALID_COMMIT, 'pointfinder': self.POINTFINDER_VALID_COMMIT,
             'plasmidfinder': self.PLASMIDFINDER_VALID_COMMIT})

        # Verify database is built properly
        self.assertTrue(path.exists(self.database_repositories.get_repo_dir('resfinder')),
                        'No resfinder dir')
        self.assertTrue(path.exists(self.database_repositories.get_repo_dir('pointfinder')),
                        'No pointfinder dir')
        self.assertTrue(path.exists(self.database_repositories.get_repo_dir('plasmidfinder')),
                        'No plasmidfinder dir')

        # Verify correct commits
        resfinder_repo_head = git.Repo(self.database_repositories.get_repo_dir('resfinder')).commit('HEAD')
        self.assertEqual(str(resfinder_repo_head), self.RESFINDER_VALID_COMMIT,
                         'Resfinder commits invalid')
        pointfinder_repo_head = git.Repo(self.database_repositories.get_repo_dir('pointfinder')).commit('HEAD')
        self.assertEqual(str(pointfinder_repo_head), self.POINTFINDER_VALID_COMMIT,
                         'Pointfinder commits invalid')
        plasmidfinder_repo_head = git.Repo(self.database_repositories.get_repo_dir('plasmidfinder')).commit('HEAD')
        self.assertEqual(str(plasmidfinder_repo_head), self.PLASMIDFINDER_VALID_COMMIT,
                         'Plasmidfinder commits invalid')

    def testUpdate(self):
        # Build database
        self.database_repositories.build(
            {'resfinder': self.RESFINDER_VALID_COMMIT, 'pointfinder': self.POINTFINDER_VALID_COMMIT,
             'plasmidfinder': self.PLASMIDFINDER_VALID_COMMIT})

        # Update database
        self.database_repositories.update(
            {'resfinder': self.RESFINDER_VALID_COMMIT2, 'pointfinder': self.POINTFINDER_VALID_COMMIT2,
             'plasmidfinder': self.PLASMIDFINDER_VALID_COMMIT2})

        # Verify correct commits
        resfinder_repo_head = git.Repo(self.database_repositories.get_repo_dir('resfinder')).commit('HEAD')
        self.assertEqual(str(resfinder_repo_head), self.RESFINDER_VALID_COMMIT2,
                         'Resfinder commits invalid')
        pointfinder_repo_head = git.Repo(self.database_repositories.get_repo_dir('pointfinder')).commit('HEAD')
        self.assertEqual(str(pointfinder_repo_head), self.POINTFINDER_VALID_COMMIT2,
                         'Pointfinder commits invalid')
        plasmidfinder_repo_head = git.Repo(self.database_repositories.get_repo_dir('plasmidfinder')).commit('HEAD')
        self.assertEqual(str(plasmidfinder_repo_head), self.PLASMIDFINDER_VALID_COMMIT2,
                         'Plasmidfinder commits invalid')

    def testInfo(self):
        # Build database
        self.database_repositories.build(
            {'resfinder': self.RESFINDER_VALID_COMMIT, 'pointfinder': self.POINTFINDER_VALID_COMMIT,
             'plasmidfinder': self.PLASMIDFINDER_VALID_COMMIT})

        database_info = self.database_repositories.info()

        # Verify correct commits in info
        self.assertEqual(database_info['resfinder_db_commit'], self.RESFINDER_VALID_COMMIT,
                         'Resfinder commits invalid')
        self.assertEqual(database_info['pointfinder_db_commit'], self.POINTFINDER_VALID_COMMIT,
                         'Pointfinder commits invalid')
        self.assertEqual(database_info['plasmidfinder_db_commit'], self.PLASMIDFINDER_VALID_COMMIT,
                         'Plasmidfinder commits invalid')
