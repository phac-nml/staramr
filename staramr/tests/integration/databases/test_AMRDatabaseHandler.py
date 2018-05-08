import os
import tempfile
import unittest
from os import path

import git

from staramr.databases.AMRDatabaseHandler import AMRDatabaseHandler
from staramr.exceptions.DatabaseErrorException import DatabaseErrorException


class AMRDatabaseHandlerIT(unittest.TestCase):
    RESFINDER_VALID_COMMIT = 'dc33e2f9ec2c420f99f77c5c33ae3faa79c999f2'
    RESFINDER_VALID_COMMIT2 = 'a4a699f3d13974477c7120b98fb0c63a1b70bd16'
    RESFINDER_ERROR_COMMIT = 'fcc4259ef58ebcd85af38509102f8a9d3a36ad18'
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

        # Verify files from `makeblastdb` command got created
        makeblastdb_resfinder_file_count = len(
            [x for x in os.listdir(self.database_handler.get_resfinder_dir()) if x.endswith('.nhr')])
        self.assertEqual(makeblastdb_resfinder_file_count, 15,
                         "Not all 'makeblastdb' commands were run properly for resfinder")
        makeblastdb_pointfinder_file_count = len(
            [x for x in os.listdir(path.join(self.database_handler.get_pointfinder_dir(), 'salmonella')) if
             x.endswith('.nhr')])
        self.assertEqual(makeblastdb_pointfinder_file_count, 8,
                         "Not all 'makeblastdb' commands were run properly for pointfinder salmonella")

        # Verify no error files
        self.assertFalse(path.exists(path.join(self.database_handler.get_database_dir(), '.error')),
                         'A .error file exists')

    def testBuildWithError(self):
        # Build database
        self.assertRaises(DatabaseErrorException, self.database_handler.build, self.RESFINDER_ERROR_COMMIT,
                          self.POINTFINDER_VALID_COMMIT)

        # Verify error file exists
        self.assertTrue(path.exists(path.join(self.database_handler.get_database_dir(), '.error')),
                        'A .error file does not exist')

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

        # Verify no error files
        self.assertFalse(path.exists(path.join(self.database_handler.get_database_dir(), '.error')),
                         'A .error file exists')

    def testUpdateWithError(self):
        # Build database
        self.database_handler.build(resfinder_commit=self.RESFINDER_VALID_COMMIT,
                                    pointfinder_commit=self.POINTFINDER_VALID_COMMIT)

        # Update database
        self.assertRaises(DatabaseErrorException, self.database_handler.build, self.RESFINDER_ERROR_COMMIT,
                          self.POINTFINDER_VALID_COMMIT)

        # Verify error file exists
        self.assertTrue(path.exists(path.join(self.database_handler.get_database_dir(), '.error')),
                        'A .error file does not exist')

    def testInfo(self):
        # Build database
        self.database_handler.build(resfinder_commit=self.RESFINDER_VALID_COMMIT,
                                    pointfinder_commit=self.POINTFINDER_VALID_COMMIT)

        database_info = self.database_handler.info()
        resfinder_commit_info = [x[1] for x in database_info if x[0] == 'resfinder_db_commit'][0]
        pointfinder_commit_info = [x[1] for x in database_info if x[0] == 'pointfinder_db_commit'][0]

        # Verify correct commits in info
        self.assertEqual(resfinder_commit_info, self.RESFINDER_VALID_COMMIT,
                         'Resfinder commits invalid')
        self.assertEqual(pointfinder_commit_info, self.POINTFINDER_VALID_COMMIT,
                         'Pointfinder commits invalid')

    def testInfoWithError(self):
        # Build database
        self.assertRaises(DatabaseErrorException, self.database_handler.build, self.RESFINDER_ERROR_COMMIT,
                          self.POINTFINDER_VALID_COMMIT)

        # Verify error file exists
        self.assertRaises(DatabaseErrorException, self.database_handler.info)
