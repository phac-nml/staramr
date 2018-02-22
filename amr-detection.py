#!/usr/bin/env python
import argparse
import logging
import sys
from os import path, mkdir

from amr.AMRDetection import AMRDetection
from amr.blast.BlastHandler import BlastHandler
from amr.blast.pointfinder.PointfinderBlastDatabase import PointfinderBlastDatabase
from amr.blast.resfinder.ResfinderBlastDatabase import ResfinderBlastDatabase

logger = logging.getLogger("amr-detection")
logging.basicConfig(level=logging.INFO, format='%(asctime)s %(levelname)s: %(message)s')

resfinder_database_dir = path.join("databases", "resfinder")
pointfinder_database_root_dir = path.join("databases", "pointfinder")


def print_to_file(dataframe, file=None):
    file_handle = sys.stdout

    if dataframe is not None:
        if file:
            file_handle = open(file, 'w')

        dataframe.to_csv(file_handle, sep="\t", float_format="%0.2f")

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

    amr_detection = AMRDetection(resfinder_database, blast_handler, pointfinder_database)
    amr_detection.run_amr_detection(args.files, args.pid_threshold, args.plength_threshold)

    if args.output_dir:
        print_to_file(amr_detection.get_resfinder_results(), path.join(args.output_dir, "results_tab.tsv"))
        print_to_file(amr_detection.get_pointfinder_results(), path.join(args.output_dir, "results_tab.pointfinder.tsv"))
        print_to_file(amr_detection.get_summary_results(), path.join(args.output_dir, "summary.tsv"))
    else:
        print_to_file(amr_detection.get_resfinder_results())
        print_to_file(amr_detection.get_pointfinder_results())
        print_to_file(amr_detection.get_summary_results())
