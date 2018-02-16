#!/usr/bin/env python
import argparse
import logging
import sys
from os import path, mkdir

from amr.blast.BlastHandler import BlastHandler
from amr.blast.pointfinder.PointfinderBlastDatabase import PointfinderBlastDatabase
from amr.blast.resfinder.ResfinderBlastDatabase import ResfinderBlastDatabase
from amr.blast.results.AMRDetectionSummary import AMRDetectionSummary
from amr.blast.results.pointfinder.BlastResultsParserPointfinder import BlastResultsParserPointfinder
from amr.blast.results.resfinder.BlastResultsParserResfinder import BlastResultsParserResfinder

logger = logging.getLogger("amr-detection")
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

resfinder_database_dir = path.join("databases", "resfinder")
pointfinder_database_root_dir = path.join("databases", "pointfinder")
arg_drug_table_resfinder_file = path.join("data", "ARG_drug_key_resfinder.tsv")
fasta_suffix = ".fsa"


def print_to_file(dataframe, file=None):
    file_handle = sys.stdout

    if file:
        file_handle = open(file, 'w')

    dataframe.to_csv(file_handle, sep="\t", index=False, float_format="%0.2f")

    if file:
        file_handle.close()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Do AMR detection for genes and point mutations')
    parser.add_argument('--threads', action='store', dest='threads', type=int, help='The number of threads to use [1].',
                        default=1, required=False)
    parser.add_argument('--pid-threshold', action='store', dest='pid_threshold', type=float,
                        help='The % identity threshold [98.0].', default=98.0, required=False)
    parser.add_argument('--percent-length-overlap', action='store', dest='plength_threshold', type=float,
                        help='The % length overlap [60.0].', default=60.0, required=False)
    parser.add_argument('--pointfinder-organism', action='store', dest='pointfinder_organism', type=str,
                        help='The organism to use for pointfinder [None].', default=None, required=False)
    parser.add_argument('--output-dir', action='store', dest='output_dir', type=str,
                        help="The output directory for results.  If unset prints all results to stdout.", default=None,
                        required=False)
    parser.add_argument('files', nargs=argparse.REMAINDER)
    args = parser.parse_args()

    if (len(args.files) == 0):
        raise Exception("Must pass a fasta file to process")

    if args.output_dir:
        if path.exists(args.output_dir):
            raise Exception("Error, output directory [" + args.output_dir + "] already exists")
        else:
            mkdir(args.output_dir)

    resfinder_database = ResfinderBlastDatabase(resfinder_database_dir)
    if (args.pointfinder_organism):
        pointfinder_database = PointfinderBlastDatabase(pointfinder_database_root_dir, args.pointfinder_organism)
    else:
        pointfinder_database = None
    blast_handler = BlastHandler(resfinder_database, pointfinder_database, threads=args.threads)

    blast_handler.run_blasts(args.files)

    resfinder_blast_map = blast_handler.get_resfinder_outputs()
    resfinder_parser = BlastResultsParserResfinder(resfinder_blast_map, resfinder_database, args.pid_threshold,
                                                   args.plength_threshold)
    resfinder_dataframe = resfinder_parser.parse_results()

    if args.output_dir:
        print_to_file(resfinder_dataframe, path.join(args.output_dir, "results_tab.tsv"))
    else:
        print_to_file(resfinder_dataframe)

    summary = None
    if (blast_handler.is_pointfinder_configured()):
        pointfinder_blast_map = blast_handler.get_pointfinder_outputs()
        pointfinder_parser = BlastResultsParserPointfinder(pointfinder_blast_map, pointfinder_database,
                                                           args.pid_threshold, args.plength_threshold)
        pointfinder_dataframe = pointfinder_parser.parse_results()

        if args.output_dir:
            print_to_file(pointfinder_dataframe, path.join(args.output_dir, "results_tab.pointfinder.tsv"))
        else:
            print_to_file(pointfinder_dataframe)

        summary = AMRDetectionSummary(resfinder_dataframe, pointfinder_dataframe)
    else:
        summary = AMRDetectionSummary(resfinder_dataframe)

    if args.output_dir:
        print_to_file(summary.create_summary(), path.join(args.output_dir, "summary.tsv"))
    else:
        print_to_file(summary.create_summary())
