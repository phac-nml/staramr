import argparse
import logging
import sys
from os import path, mkdir

from amr.detection.AMRDetection import AMRDetection
from amr.blast.BlastHandler import BlastHandler
from amr.blast.pointfinder.PointfinderBlastDatabase import PointfinderBlastDatabase
from amr.blast.resfinder.ResfinderBlastDatabase import ResfinderBlastDatabase
from amr.SubCommand import SubCommand
from amr.exceptions.CommandParseException import CommandParseException

class Search(SubCommand):

    def __init__(self, arg_parser, default_resfinder_dir, default_pointfinder_dir):
        super().__init__(arg_parser)
        self._resfinder_database_dir = default_resfinder_dir
        self._pointfinder_database_root_dir = default_pointfinder_dir

    def _setup_args(self, arg_parser):
        arg_parser.add_argument('--threads', action='store', dest='threads', type=int, help='The number of threads to use [1].',
                            default=1, required=False)
        arg_parser.add_argument('--pid-threshold', action='store', dest='pid_threshold', type=float,
                            help='The percent identity threshold [98.0].', default=98.0, required=False)
        arg_parser.add_argument('--percent-length-overlap', action='store', dest='plength_threshold', type=float,
                            help='The percent length overlap [60.0].', default=60.0, required=False)
        arg_parser.add_argument('--pointfinder-organism', action='store', dest='pointfinder_organism', type=str,
                            help='The organism to use for pointfinder [None].', default=None, required=False)
        arg_parser.add_argument('--output-dir', action='store', dest='output_dir', type=str,
                            help="The output directory for results.  If unset prints all results to stdout.", default=None,
                            required=False)
        arg_parser.add_argument('files', nargs=argparse.REMAINDER)

    def _print_dataframe_to_file(self, dataframe, file=None):
        file_handle = sys.stdout

        if dataframe is not None:
            if file:
                file_handle = open(file, 'w')

            dataframe.to_csv(file_handle, sep="\t", float_format="%0.2f")

            if file:
                file_handle.close()

    def run(self, args):
        if (len(args.files) == 0):
            raise CommandParseException("Must pass a fasta file to process", self._root_arg_parser)

        if args.output_dir:
            if path.exists(args.output_dir):
                raise CommandParseException("Error, output directory [" + args.output_dir + "] already exists", self._root_arg_parser)
            else:
                mkdir(args.output_dir)

        resfinder_database = ResfinderBlastDatabase(self._resfinder_database_dir)
        if (args.pointfinder_organism):
            pointfinder_database = PointfinderBlastDatabase(self._pointfinder_database_root_dir, args.pointfinder_organism)
        else:
            pointfinder_database = None
        blast_handler = BlastHandler(resfinder_database, pointfinder_database, threads=args.threads)

        amr_detection = AMRDetection(resfinder_database, blast_handler, pointfinder_database)
        amr_detection.run_amr_detection(args.files, args.pid_threshold, args.plength_threshold)

        if args.output_dir:
            self._print_dataframe_to_file(amr_detection.get_resfinder_results(), path.join(args.output_dir, "results_tab.tsv"))
            self._print_dataframe_to_file(amr_detection.get_pointfinder_results(),
                          path.join(args.output_dir, "results_tab.pointfinder.tsv"))
            self._print_dataframe_to_file(amr_detection.get_summary_results(), path.join(args.output_dir, "summary.tsv"))
        else:
            self._print_dataframe_to_file(amr_detection.get_resfinder_results())
            self._print_dataframe_to_file(amr_detection.get_pointfinder_results())
            self._print_dataframe_to_file(amr_detection.get_summary_results())