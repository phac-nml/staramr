import argparse
from os import path,mkdir

from amr.SubCommand import SubCommand
from amr.databases.AMRDatabaseHandler import AMRDatabaseHandler

class Build(SubCommand):

    def __init__(self, arg_parser):
        super().__init__(arg_parser)

    def _setup_args(self, arg_parser):
        arg_parser.add_argument('--dir', action='store', dest='destination', type=str, help='The directory to download the databases into [databases].',
                            default='databases', required=False)

    def run(self, args):
        if path.exists(args.destination):
            raise Exception("Error, destination [" + args.destination + "] already exists")
        else:
            mkdir(args.destination)

        database_handler = AMRDatabaseHandler(args.destination)
        database_handler.build()


