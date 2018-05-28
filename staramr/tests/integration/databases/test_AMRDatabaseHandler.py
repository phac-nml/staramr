import tempfile
import unittest
from os import path

import git

from staramr.databases.AMRDatabaseHandler import AMRDatabaseHandler


class AMRDatabaseHandlerIT(unittest.TestCase):
    RESFINDER_VALID_COMMIT = 'dc33e2f9ec2c420f99f77c5c33ae3faa79c999f2'
    RESFINDER_VALID_COMMIT2 = 'a4a699f3d13974477c7120b98fb0c63a1b70bd16'
    POINTFINDER_VALID_COMMIT = 'ba65c4d175decdc841a0bef9f9be1c1589c0070a'
    POINTFINDER_VALID_COMMIT2 = '0de22bff78214208171aef70461c639227e62e5d'

    def setUp(self):
        self.databases_dir = tempfile.TemporaryDirectory()
        self.database_handler = AMRDatabaseHandler(database_dir=self.databases_dir.name)

    def tearDown(self):
        self.databases_dir.cleanup()

    def testBuild(self):
        # Verify that databases don't exist beforehand
        self.assertFalse(path.exists(self.database_handler.get_resfinder_dir()),
                         'resfinder path exists before creation of database')
        self.assertFalse(path.exists(self.database_handler.get_pointfinder_dir()),
                         'pointfinder path exists before creation of database')

        # Build database
        self.database_handler.build(resfinder_commit=self.RESFINDER_VALID_COMMIT,
                                    pointfinder_commit=self.POINTFINDER_VALID_COMMIT)

        # Verify database is built properly
        self.assertTrue(path.exists(self.database_handler.get_resfinder_dir()),
                        'No resfinder dir')
        self.assertTrue(path.exists(self.database_handler.get_pointfinder_dir()),
                        'No pointfinder dir')

        # Verify correct commits
        resfinder_repo_head = git.Repo(self.database_handler.get_resfinder_dir()).commit('HEAD')
        self.assertEqual(str(resfinder_repo_head), self.RESFINDER_VALID_COMMIT,
                         'Resfinder commits invalid')
        pointfinder_repo_head = git.Repo(self.database_handler.get_pointfinder_dir()).commit('HEAD')
        self.assertEqual(str(pointfinder_repo_head), self.POINTFINDER_VALID_COMMIT,
                         'Pointfinder commits invalid')

    def testUpdate(self):
        # Build database
        self.database_handler.build(resfinder_commit=self.RESFINDER_VALID_COMMIT,
                                    pointfinder_commit=self.POINTFINDER_VALID_COMMIT)

        # Update database
        self.database_handler.update(resfinder_commit=self.RESFINDER_VALID_COMMIT2,
                                     pointfinder_commit=self.POINTFINDER_VALID_COMMIT2)

        # Verify correct commits
        resfinder_repo_head = git.Repo(self.database_handler.get_resfinder_dir()).commit('HEAD')
        self.assertEqual(str(resfinder_repo_head), self.RESFINDER_VALID_COMMIT2,
                         'Resfinder commits invalid')
        pointfinder_repo_head = git.Repo(self.database_handler.get_pointfinder_dir()).commit('HEAD')
        self.assertEqual(str(pointfinder_repo_head), self.POINTFINDER_VALID_COMMIT2,
                         'Pointfinder commits invalid')

    def testInfo(self):
        # Build database
        self.database_handler.build(resfinder_commit=self.RESFINDER_VALID_COMMIT,
                                    pointfinder_commit=self.POINTFINDER_VALID_COMMIT)

        database_info = self.database_handler.info()

        # Verify correct commits in info
        self.assertEqual(database_info['resfinder_db_commit'], self.RESFINDER_VALID_COMMIT,
                         'Resfinder commits invalid')
        self.assertEqual(database_info['pointfinder_db_commit'], self.POINTFINDER_VALID_COMMIT,
                         'Pointfinder commits invalid')
