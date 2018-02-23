#!/usr/bin/env python
import argparse
import logging
from os import path
import sys

from amr.subcommand.Search import Search
from amr.subcommand.Build import Build

logger = logging.getLogger("amr-detection")
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

script_dir = path.dirname(path.realpath(sys.argv[0]))

default_database_dir = path.join(script_dir, "databases")
resfinder_database_dir = path.join(default_database_dir, "resfinder")
pointfinder_database_root_dir = path.join(default_database_dir, "pointfinder")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Do AMR detection for genes and point mutations')
    subparsers = parser.add_subparsers(dest='command', help='sub-command --help')

    search = Search(subparsers.add_parser('search', help='Search for AMR genes'), resfinder_database_dir, pointfinder_database_root_dir)
    db_download = Build(subparsers.add_parser('build', help='Download ResFinder/PointFinder databases'))

    args = parser.parse_args()
    if args.command is None:
        parser.print_help()
    else:
        args.run_command(args)
