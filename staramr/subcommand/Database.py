"""
Classes for interacting with the (ResFinder/PointFinder) databases used to detect AMR genes.
"""
import argparse
import sys
from os import path, mkdir

from staramr.SubCommand import SubCommand
from staramr.Utils import get_string_with_spacing
from staramr.databases.AMRDatabaseHandlerFactory import AMRDatabaseHandlerFactory
from staramr.exceptions.CommandParseException import CommandParseException

"""
Base class for interacting with a database.
"""


class Database(SubCommand):

    def __init__(self, subparser, script_name):
        """
        Builds a SubCommand for interacting with databases.
        :param subparser: The subparser to use.  Generated from argparse.ArgumentParser.add_subparsers().
        :param script_name: The name of the script being run.
        """
        super().__init__(subparser, script_name)

    def _setup_args(self, arg_parser):
        arg_parser = self._subparser.add_parser('db', help='Download ResFinder/PointFinder databases')
        subparser = arg_parser.add_subparsers(dest='db_command',
                                              help='Subcommand for ResFinder/PointFinder databases.')
        arg_parser.add_argument('--version', action='store_true', dest='version',
                                help='Prints version information.', required=False)

        Build(subparser, self._script_name + " db")
        Update(subparser, self._script_name + " db")
        Info(subparser, self._script_name + " db")

        return arg_parser

    def run(self, args):
        super(Database, self).run(args)

        if args.db_command is None:
            self._root_arg_parser.print_help()


"""
Class for building a new database.
"""


class Build(Database):

    def __init__(self, subparser, script_name):
        """
        Creates a SubCommand for building a new database.
        :param subparser: The subparser to use.  Generated from argparse.ArgumentParser.add_subparsers().
        :param script_name: The name of the script being run.
        """
        super().__init__(subparser, script_name)

    def _setup_args(self, arg_parser):
        name = self._script_name
        default_dir = AMRDatabaseHandlerFactory.get_default_database_directory()
        epilog = ("Example:\n"
                  "\t" + name + " build\n"
                                "\t\tBuilds a new ResFinder/PointFinder database under " + default_dir + " if it does not exist\n\n" +
                  "\t" + name + " build --dir databases\n" +
                  "\t\tBuilds a new ResFinder/PointFinder database under databases/")

        arg_parser = self._subparser.add_parser('build',
                                                epilog=epilog,
                                                formatter_class=argparse.RawTextHelpFormatter,
                                                help='Downloads and builds databases in the given directory.')
        arg_parser.add_argument('--dir', action='store', dest='destination', type=str,
                                help='The directory to download the databases into [' + default_dir + '].',
                                default=default_dir, required=False)
        arg_parser.add_argument('--resfinder-commit', action='store', dest='resfinder_commit', type=str,
                                help='The specific git commit for the resfinder database [latest].', required=False)
        arg_parser.add_argument('--pointfinder-commit', action='store', dest='pointfinder_commit', type=str,
                                help='The specific git commit for the pointfinder database [latest].', required=False)
        return arg_parser

    def run(self, args):
        super(Build, self).run(args)

        if path.exists(args.destination):
            raise CommandParseException("Error, destination [" + args.destination + "] already exists",
                                        self._root_arg_parser)
        else:
            mkdir(args.destination)

        if args.destination == AMRDatabaseHandlerFactory.get_default_database_directory():
            database_handler = AMRDatabaseHandlerFactory.create_default_factory().get_database_handler()
        else:
            database_handler = AMRDatabaseHandlerFactory(args.destination).get_database_handler()
        database_handler.build(resfinder_commit=args.resfinder_commit, pointfinder_commit=args.pointfinder_commit)


"""
Class for updating an existing database.
"""


class Update(Database):

    def __init__(self, subparser, script_name):
        """
        Creates a SubCommand for updating an existing database.
        :param subparser: The subparser to use.  Generated from argparse.ArgumentParser.add_subparsers().
        :param script_name: The name of the script being run.
        """
        super().__init__(subparser, script_name)

    def _setup_args(self, arg_parser):
        default_dir = AMRDatabaseHandlerFactory.get_default_database_directory()
        name = self._script_name
        epilog = ("Example:\n"
                  "\t" + name + " update databases/\n"
                                "\t\tUpdates the ResFinder/PointFinder database under databases/\n\n" +
                  "\t" + name + " update -d\n" +
                  "\t\tUpdates the default ResFinder/PointFinder database under " + default_dir)
        arg_parser = self._subparser.add_parser('update',
                                                epilog=epilog,
                                                formatter_class=argparse.RawTextHelpFormatter,
                                                help='Updates databases in the given directories.')

        arg_parser.add_argument('-d', '--update-default', action='store_true', dest='update_default',
                                help='Updates default database directory (' + default_dir + ').', required=False)
        arg_parser.add_argument('--resfinder-commit', action='store', dest='resfinder_commit', type=str,
                                help='The specific git commit for the resfinder database [latest].', required=False)
        arg_parser.add_argument('--pointfinder-commit', action='store', dest='pointfinder_commit', type=str,
                                help='The specific git commit for the pointfinder database [latest].', required=False)
        arg_parser.add_argument('directories', nargs=argparse.REMAINDER)

        return arg_parser

    def run(self, args):
        super(Update, self).run(args)

        if len(args.directories) == 0:
            if not args.update_default:
                raise CommandParseException("Must pass at least one directory to update", self._root_arg_parser)
            else:
                database_handler = AMRDatabaseHandlerFactory.create_default_factory().get_database_handler(
                    force_use_git=True)
                database_handler.update(resfinder_commit=args.resfinder_commit,
                                        pointfinder_commit=args.pointfinder_commit)
        else:
            for directory in args.directories:
                database_handler = AMRDatabaseHandlerFactory(directory).get_database_handler()
                database_handler.update(resfinder_commit=args.resfinder_commit,
                                        pointfinder_commit=args.pointfinder_commit)


"""
Class for getting information from an existing database.
"""


class Info(Database):

    def __init__(self, subparser, script_name):
        """
        Creates a SubCommand for printing information about a database.
        :param subparser: The subparser to use.  Generated from argparse.ArgumentParser.add_subparsers().
        :param script_dir: The directory containing the main application script.
        :param script_name: The name of the script being run.
        """
        super().__init__(subparser, script_name)

    def _setup_args(self, arg_parser):
        name = self._script_name
        default_dir = AMRDatabaseHandlerFactory.get_default_database_directory()
        epilog = ("Example:\n"
                  "\t" + name + " info\n"
                                "\t\tPrints information about the default database in " + default_dir + "\n\n" +
                  "\t" + name + " info databases\n" +
                  "\t\tPrints information on the database stored in databases/")
        arg_parser = self._subparser.add_parser('info',
                                                epilog=epilog,
                                                formatter_class=argparse.RawTextHelpFormatter,
                                                help='Prints information on databases in the given directories.')
        arg_parser.add_argument('directories', nargs=argparse.REMAINDER)

        return arg_parser

    def run(self, args):
        super(Info, self).run(args)

        if len(args.directories) == 0:
            database_handler = AMRDatabaseHandlerFactory.create_default_factory().get_database_handler()
            sys.stdout.write(get_string_with_spacing(database_handler.info()))
        elif len(args.directories) == 1:
            database_handler = AMRDatabaseHandlerFactory(args.directories[0]).get_database_handler()
            sys.stdout.write(get_string_with_spacing(database_handler.info()))
        else:
            for directory in args.directories:
                database_handler = AMRDatabaseHandlerFactory(directory).get_database_handler()
                sys.stdout.write(get_string_with_spacing(database_handler.info()))
