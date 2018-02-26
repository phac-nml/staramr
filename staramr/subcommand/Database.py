import argparse
from os import path,mkdir

from staramr.SubCommand import SubCommand
from staramr.exceptions.CommandParseException import CommandParseException
from staramr.databases.AMRDatabaseHandler import AMRDatabaseHandler

class Database(SubCommand):

    def __init__(self, arg_parser):
        super().__init__(arg_parser)

    def _setup_args(self, arg_parser):
        subparsers = arg_parser.add_subparsers(dest='db_command', help='Subcommand for ResFinder/PointFinder databases.')

        Build(subparsers.add_parser('build', help='Downloads and builds databases in the given directory.'))
        Update(subparsers.add_parser('update', help='Updates databases in the given directories.'))
        Info(subparsers.add_parser('info', help='Prints information on databases in the given directories.'))

    def run(self, args):
        if args.db_command is None:
            self._root_arg_parser.print_help()

class Build(Database):

    def __init__(self, arg_parser):
        super().__init__(arg_parser)

    def _setup_args(self, arg_parser):
        arg_parser.add_argument('--dir', action='store', dest='destination', type=str, help='The directory to download the databases into [databases].',
                            default='databases', required=False)

    def run(self, args):
        super(Build, self).run(args)

        if path.exists(args.destination):
            raise CommandParseException("Error, destination [" + args.destination + "] already exists", self._root_arg_parser)
        else:
            mkdir(args.destination)

        database_handler = AMRDatabaseHandler(args.destination)
        database_handler.build()

class Update(Database):

    def __init__(self, arg_parser):
        super().__init__(arg_parser)

    def _setup_args(self, arg_parser):
        arg_parser.add_argument('directories', nargs=argparse.REMAINDER)

    def run(self, args):
        super(Update, self).run(args)

        if len(args.directories) == 0:
            raise CommandParseException("Must pass at least one directory to update", self._root_arg_parser)
        else:
            for directory in args.directories:
                database_handler = AMRDatabaseHandler(directory)
                database_handler.update()

class Info(Database):

    def __init__(self, arg_parser):
        super().__init__(arg_parser)

    def _setup_args(self, arg_parser):
        arg_parser.add_argument('directories', nargs=argparse.REMAINDER)

    def run(self, args):
        super(Info, self).run(args)

        if len(args.directories) == 0:
            raise CommandParseException("Must pass at least one directory", self._root_arg_parser)
        if len(args.directories) == 1:
            database_handler = AMRDatabaseHandler(args.directories[0])
            database_handler.info()
        else:
            for directory in args.directories:
                database_handler = AMRDatabaseHandler(directory)
                print()
                database_handler.info()