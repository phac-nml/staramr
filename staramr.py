#!/usr/bin/env python
import argparse
import logging
import sys
from os import path

from staramr.exceptions.CommandParseException import CommandParseException
from staramr.subcommand.Database import Database
from staramr.subcommand.Search import Search

logger = logging.getLogger("staramr-detection")
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

script_dir = path.dirname(path.realpath(sys.argv[0]))

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Do AMR detection for genes and point mutations')
    subparsers = parser.add_subparsers(dest='command', help='Subcommand for AMR detection.')

    Search(subparsers.add_parser('search', help='Search for AMR genes'), script_dir)
    Database(subparsers.add_parser('db', help='Download ResFinder/PointFinder databases'), script_dir)

    args = parser.parse_args()
    if args.command is None:
        parser.print_help()
    else:
        try:
            args.run_command(args)
        except CommandParseException as e:
            logger.error(str(e))
            e.get_parser().print_help()
