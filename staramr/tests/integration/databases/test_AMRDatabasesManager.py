import configparser
import tempfile
import unittest
from os import path

from staramr.databases.AMRDatabasesManager import AMRDatabasesManager


class AMRDatabasesManagerIT(unittest.TestCase):
    RESFINDER_DEFAULT_COMMIT = 'e8f1eb2585cd9610c4034a54ce7fc4f93aa95535'
    POINTFINDER_DEFAULT_COMMIT = '8706a6363bb29e47e0e398c53043b037c24b99a7'
    PLASMIDFINDER_DEFAULT_COMMIT = '81919954cbedaff39056610ab584ab4c06011ed8'

    def setUp(self):
        self.databases_dir = tempfile.TemporaryDirectory()
        self.databases_manager = AMRDatabasesManager(database_dir=self.databases_dir.name, sub_dirs=True)

    def tearDown(self):
        self.databases_dir.cleanup()

    def testSetupDefault(self):
        blast_database_repos = self.databases_manager.get_database_repos()

        # Verify that databases don't exist beforehand
        self.assertFalse(path.exists(blast_database_repos.get_repo_dir('resfinder')),
                         'resfinder path exists before creation of database')
        self.assertFalse(path.exists(blast_database_repos.get_repo_dir('pointfinder')),
                         'pointfinder path exists before creation of database')
        self.assertFalse(path.exists(blast_database_repos.get_repo_dir('plasmidfinder')),
                         'plasmidfinder path exists before creation of database')

        # Setup default database
        self.databases_manager.setup_default()

        # Verify that resfinder/pointfinder paths exist
        self.assertTrue(path.exists(blast_database_repos.get_repo_dir('resfinder')), 'resfinder path does not exist')
        self.assertTrue(path.exists(blast_database_repos.get_repo_dir('pointfinder')),
                        'pointfinder path does not exist')
        self.assertTrue(path.exists(blast_database_repos.get_repo_dir('plasmidfinder')),
                        'plasmidfinder path does not exist')

        self.assertTrue(path.exists(path.join(blast_database_repos.get_database_dir(), 'resfinder-info.ini')),
                        'resfinder info file does not exist')
        self.assertTrue(path.exists(path.join(blast_database_repos.get_database_dir(), 'pointfinder-info.ini')),
                        'pointfinder info file does not exist')
        self.assertTrue(path.exists(path.join(blast_database_repos.get_database_dir(), 'plasmidfinder-info.ini')),
                        'plasmidfinder info file does not exist')

        # Verify we've removed the .git directories
        self.assertFalse(path.exists(path.join(blast_database_repos.get_repo_dir('resfinder'), '.git')),
                         'resfinder .git directory was not removed')
        self.assertFalse(path.exists(path.join(blast_database_repos.get_repo_dir('pointfinder'), '.git')),
                         'pointfinder .git directory was not removed')
        self.assertFalse(path.exists(path.join(blast_database_repos.get_repo_dir('plasmidfinder'), '.git')),
                         'plasmidfinder .git directory was not removed')

        config = configparser.ConfigParser()

        # Verify that the info.ini file has correct git commits for default database
        config.read(path.join(blast_database_repos.get_database_dir(), 'resfinder-info.ini'))
        self.assertEqual(config['GitInfo']['resfinder_db_commit'], self.RESFINDER_DEFAULT_COMMIT,
                         'invalid resfinder commit')

        config.read(path.join(blast_database_repos.get_database_dir(), 'pointfinder-info.ini'))
        self.assertEqual(config['GitInfo']['pointfinder_db_commit'], self.POINTFINDER_DEFAULT_COMMIT,
                         'invalid pointfinder commit')

        config.read(path.join(blast_database_repos.get_database_dir(), 'plasmidfinder-info.ini'))
        self.assertEqual(config['GitInfo']['plasmidfinder_db_commit'], self.PLASMIDFINDER_DEFAULT_COMMIT,
                         'invalid plasmidfinder commit')

    def testRestoreDefault(self):
        # Build initial default database
        self.databases_manager.setup_default()

        # Build updated database
        blast_database_repos_git = self.databases_manager.get_database_repos(force_use_git=True)
        blast_database_repos_git.build(
            {'resfinder': self.RESFINDER_DEFAULT_COMMIT, 'pointfinder': self.POINTFINDER_DEFAULT_COMMIT,
             'plasmidfinder': self.PLASMIDFINDER_DEFAULT_COMMIT})

        # Verify that updated database is the one that gets returned by get_database_handler()
        blast_database_repos = self.databases_manager.get_database_repos()
        self.assertFalse(blast_database_repos.is_dist(), 'Invalid is_dist')
        self.assertEqual(blast_database_repos.get_database_dir(), path.join(self.databases_dir.name, 'update'),
                         'Invalid database directory')
        self.assertTrue(path.exists(path.join(blast_database_repos.get_repo_dir('resfinder'), '.git')),
                        'Not using git version (updated version) of resfinder database')
        self.assertTrue(path.exists(path.join(blast_database_repos.get_repo_dir('pointfinder'), '.git')),
                        'Not using git version (updated version) of pointfinder database')
        self.assertTrue(path.exists(path.join(blast_database_repos.get_repo_dir('plasmidfinder'), '.git')),
                        'Not using git version (updated version) of plasmidfinder database')

        # Restore default database
        self.databases_manager.restore_default()

        # Verify that default database (git stripped version) is the one that gets returned by get_database_handler()
        blast_database_repos = self.databases_manager.get_database_repos()
        self.assertTrue(blast_database_repos.is_dist(), 'Invalid is_dist')
        self.assertEqual(blast_database_repos.get_database_dir(), path.join(self.databases_dir.name, 'dist'),
                         'Invalid database directory')
        self.assertFalse(path.exists(path.join(blast_database_repos.get_repo_dir('resfinder'), '.git')),
                         'resfinder .git directory was not removed')
        self.assertFalse(path.exists(path.join(blast_database_repos.get_repo_dir('pointfinder'), '.git')),
                         'pointfinder .git directory was not removed')
        self.assertFalse(path.exists(path.join(blast_database_repos.get_repo_dir('plasmidfinder'), '.git')),
                         'plasmidfinder .git directory was not removed')

    def testIsHandlerDefaultCommitsTrue(self):
        # Setup default database
        self.databases_manager.setup_default()

        blast_database_repos = self.databases_manager.get_database_repos()

        self.assertTrue(AMRDatabasesManager.is_database_repos_default_commits(blast_database_repos),
                        "Database is not default")

    def testIsHandlerDefaultCommitsFalse(self):
        # Setup database
        blast_database_repos = self.databases_manager.get_database_repos(force_use_git=True)
        blast_database_repos.update({'resfinder': 'dc33e2f9ec2c420f99f77c5c33ae3faa79c999f2'})

        self.assertFalse(AMRDatabasesManager.is_database_repos_default_commits(blast_database_repos),
                         "Database is default")
