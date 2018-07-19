import configparser
import tempfile
import unittest
from os import path

from staramr.databases.AMRDatabaseHandler import AMRDatabaseHandler
from staramr.databases.AMRDatabaseHandlerStripGitDir import AMRDatabaseHandlerStripGitDir
from staramr.databases.AMRDatabasesManager import AMRDatabasesManager


class AMRDatabasesManagerIT(unittest.TestCase):
    RESFINDER_DEFAULT_COMMIT = 'e8f1eb2585cd9610c4034a54ce7fc4f93aa95535'
    POINTFINDER_DEFAULT_COMMIT = '8706a6363bb29e47e0e398c53043b037c24b99a7'

    def setUp(self):
        self.databases_dir = tempfile.TemporaryDirectory()
        self.databases_manager = AMRDatabasesManager(database_dir=self.databases_dir.name, sub_dirs=True)

    def tearDown(self):
        self.databases_dir.cleanup()

    def testGetHandlerGitStripDir(self):
        self.assertIsInstance(self.databases_manager.get_database_handler(), AMRDatabaseHandlerStripGitDir,
                              'Invalid instance returned')

    def testGetHandlerGit(self):
        self.assertIsInstance(self.databases_manager.get_database_handler(force_use_git=True), AMRDatabaseHandler,
                              'Invalid instance returned')

    def testSetupDefault(self):
        database_handler = self.databases_manager.get_database_handler()

        # Verify that databases don't exist beforehand
        self.assertFalse(path.exists(database_handler.get_resfinder_dir()),
                         'resfinder path exists before creation of database')
        self.assertFalse(path.exists(database_handler.get_pointfinder_dir()),
                         'pointfinder path exists before creation of database')

        # Setup default database
        self.databases_manager.setup_default()

        # Verify that resfinder/pointfinder paths exist
        self.assertTrue(path.exists(database_handler.get_resfinder_dir()), 'resfinder path does not exist')
        self.assertTrue(path.exists(database_handler.get_resfinder_dir()), 'pointfinder path does not exist')
        self.assertTrue(path.exists(path.join(database_handler.get_database_dir(), 'info.ini')),
                        'info file does not exist')

        # Verify we've removed the .git directories
        self.assertFalse(path.exists(path.join(database_handler.get_resfinder_dir(), '.git')),
                         'resfinder .git directory was not removed')
        self.assertFalse(path.exists(path.join(database_handler.get_pointfinder_dir(), '.git')),
                         'pointfinder .git directory was not removed')

        config = configparser.ConfigParser()
        config.read(path.join(database_handler.get_database_dir(), 'info.ini'))

        # Verify that the info.ini file has correct git commits for default database
        self.assertEqual(config['GitInfo']['resfinder_db_commit'], self.RESFINDER_DEFAULT_COMMIT,
                         'invalid resfinder commit')
        self.assertEqual(config['GitInfo']['pointfinder_db_commit'], self.POINTFINDER_DEFAULT_COMMIT,
                         'invalid pointfinder commit')

    def testRestoreDefault(self):
        # Build initial default database
        self.databases_manager.setup_default()

        # Build updated database
        database_handler_git = self.databases_manager.get_database_handler(force_use_git=True)
        database_handler_git.build(resfinder_commit=self.RESFINDER_DEFAULT_COMMIT,
                                   pointfinder_commit=self.POINTFINDER_DEFAULT_COMMIT)

        # Verify that updated database is the one that gets returned by get_database_handler()
        database_handler = self.databases_manager.get_database_handler()
        self.assertIsInstance(database_handler, AMRDatabaseHandler, 'Invalid instance returned')
        self.assertTrue(path.exists(path.join(database_handler.get_resfinder_dir(), '.git')),
                        'Not using git version (updated version) of resfinder database')
        self.assertTrue(path.exists(path.join(database_handler.get_pointfinder_dir(), '.git')),
                        'Not using git version (updated version) of pointfinder database')

        # Restore default database
        self.databases_manager.restore_default()

        # Verify that default database (git stripped version) is the one that gets returned by get_database_handler()
        database_handler = self.databases_manager.get_database_handler()
        self.assertIsInstance(database_handler, AMRDatabaseHandlerStripGitDir, 'Invalid instance returned')
        self.assertFalse(path.exists(path.join(database_handler.get_resfinder_dir(), '.git')),
                         'resfinder .git directory was not removed')
        self.assertFalse(path.exists(path.join(database_handler.get_pointfinder_dir(), '.git')),
                         'pointfinder .git directory was not removed')

    def testIsHandlerDefaultCommitsTrue(self):
        # Setup default database
        self.databases_manager.setup_default()

        database_handler = self.databases_manager.get_database_handler()

        self.assertTrue(AMRDatabasesManager.is_handler_default_commits(database_handler), "Database is not default")

    def testIsHandlerDefaultCommitsFalse(self):
        # Setup database
        database_handler = self.databases_manager.get_database_handler(force_use_git=True)
        database_handler.update(resfinder_commit='dc33e2f9ec2c420f99f77c5c33ae3faa79c999f2')

        self.assertFalse(AMRDatabasesManager.is_handler_default_commits(database_handler), "Database is default")
