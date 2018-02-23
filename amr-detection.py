#!/usr/bin/env python
import argparse
import logging
from os import path

from amr.subcommand.Search import Search

logger = logging.getLogger("amr-detection")
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

resfinder_database_dir = path.join("databases", "resfinder")
pointfinder_database_root_dir = path.join("databases", "pointfinder")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Do AMR detection for genes and point mutations')
    subparsers = parser.add_subparsers(help='sub-command help')
    search = Search(subparsers.add_parser('search', help='Search for AMR genes'), resfinder_database_dir, pointfinder_database_root_dir)

    args = parser.parse_args()
    args.func(args)
