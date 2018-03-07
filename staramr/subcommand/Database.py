"""
Classes for interacting with the (ResFinder/PointFinder) databases used to detect AMR genes.
"""
import argparse
from os import path, mkdir

from staramr.SubCommand import SubCommand
from staramr.databases.AMRDatabaseHandler import AMRDatabaseHandler
from staramr.exceptions.CommandParseException import CommandParseException

"""
Base class for interacting with a database.
"""


class Database(SubCommand):

    def __init__(self, arg_parser, script_dir):
        """
        Builds a SubCommand for interacting with databases.
        :param arg_parser: The argparse.ArgumentParser to use.
        :param script_dir: The directory containing the main application script.
        """
        super().__init__(arg_parser, script_dir)

    def _setup_args(self, arg_parser):
        subparsers = arg_parser.add_subparsers(dest='db_command',
                                               help='Subcommand for ResFinder/PointFinder databases.')

        Build(subparsers.add_parser('build', help='Downloads and builds databases in the given directory.'),
              self._script_dir)
        Update(subparsers.add_parser('update', help='Updates databases in the given directories.'), self._script_dir)
        Info(subparsers.add_parser('info', help='Prints information on databases in the given directories.'),
             self._script_dir)

    def run(self, args):
        if args.db_command is None:
            self._root_arg_parser.print_help()


"""
Class for building a new database.
"""


class Build(Database):

    def __init__(self, arg_parser, script_dir):
        """
        Creates a SubCommand for building a new database.
        :param arg_parser: The argparse.ArgumentParser to use.
        :param script_dir: The directory containing the main application script.
        """
        super().__init__(arg_parser, script_dir)

    def _setup_args(self, arg_parser):
        default_dir = AMRDatabaseHandler.get_default_database_directory(self._script_dir)
        arg_parser.add_argument('--dir', action='store', dest='destination', type=str,
                                help='The directory to download the databases into [' + default_dir + '].',
                                default=default_dir, required=False)

    def run(self, args):
        super(Build, self).run(args)

        if path.exists(args.destination):
            raise CommandParseException("Error, destination [" + args.destination + "] already exists",
                                        self._root_arg_parser)
        else:
            mkdir(args.destination)

        database_handler = AMRDatabaseHandler(args.destination)
        database_handler.build()


"""
Class for updating an existing database.
"""


class Update(Database):

    def __init__(self, arg_parser, script_dir):
        """
        Creates a SubCommand for updating an existing database.
        :param arg_parser: The argparse.ArgumentParser to use.
        :param script_dir: The directory containing the main application script.
        """
        super().__init__(arg_parser, script_dir)

    def _setup_args(self, arg_parser):
        default_dir = AMRDatabaseHandler.get_default_database_directory(self._script_dir)
        arg_parser.add_argument('-d', '--update-default', action='store_true', dest='update_default',
                                help='Updates default database directory (' + default_dir + ').', required=False)
        arg_parser.add_argument('directories', nargs=argparse.REMAINDER)

    def run(self, args):
        super(Update, self).run(args)

        if len(args.directories) == 0:
            if not args.update_default:
                raise CommandParseException("Must pass at least one directory to update", self._root_arg_parser)
            else:
                database_handler = AMRDatabaseHandler.create_default_handler(self._script_dir)
                database_handler.update()
        else:
            for directory in args.directories:
                database_handler = AMRDatabaseHandler(directory)
                database_handler.update()


"""
Class for getting information from an existing database.
"""


class Info(Database):

    def __init__(self, arg_parser, script_dir):
        """
        Creates a SubCommand for printing information about a database.
        :param arg_parser: The argparse.ArgumentParser to use.
        :param script_dir: The directory containing the main application script.
        """
        super().__init__(arg_parser, script_dir)

    def _setup_args(self, arg_parser):
        arg_parser.add_argument('directories', nargs=argparse.REMAINDER)

    def run(self, args):
        super(Info, self).run(args)

        if len(args.directories) == 0:
            database_handler = AMRDatabaseHandler.create_default_handler(self._script_dir)
            database_handler.info()
        elif len(args.directories) == 1:
            database_handler = AMRDatabaseHandler(args.directories[0])
            database_handler.info()
        else:
            for directory in args.directories:
                database_handler = AMRDatabaseHandler(directory)
                database_handler.info()
                print()
