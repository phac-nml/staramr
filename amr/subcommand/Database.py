import argparse
from os import path,mkdir

from amr.SubCommand import SubCommand
from amr.databases.AMRDatabaseHandler import AMRDatabaseHandler

class Database(SubCommand):

    def __init__(self, arg_parser):
        super().__init__(arg_parser)

    def _setup_args(self, arg_parser):
        subparsers = arg_parser.add_subparsers(dest='db_command', help='db --help')

        Build(subparsers.add_parser('build', help='Builds ResFinder/PointFinder databases'))
        Update(subparsers.add_parser('update', help='Updates ResFinder/PointFinder databases'))
        Info(subparsers.add_parser('info', help='Prints information on ResFinder/PointFinder databases'))

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
            raise Exception("Error, destination [" + args.destination + "] already exists")
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

        if len(args.directories) == 1:
            database_handler = AMRDatabaseHandler(args.directories[0])
            database_handler.info()
        else:
            for directory in args.directories:
                database_handler = AMRDatabaseHandler(directory)
                print()
                database_handler.info()