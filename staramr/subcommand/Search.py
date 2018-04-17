import argparse
import datetime
import logging
import multiprocessing
import sys
from os import path, mkdir

import pandas

from staramr.SubCommand import SubCommand
from staramr.Utils import get_string_with_spacing
from staramr.blast.BlastHandler import BlastHandler
from staramr.blast.pointfinder.PointfinderBlastDatabase import PointfinderBlastDatabase
from staramr.blast.resfinder.ResfinderBlastDatabase import ResfinderBlastDatabase
from staramr.databases.AMRDatabaseHandlerFactory import AMRDatabaseHandlerFactory
from staramr.exceptions.CommandParseException import CommandParseException

logger = logging.getLogger("Search")

"""
Class for searching for AMR resistance genes.
"""


class Search(SubCommand):
    blank = '-'

    def __init__(self, amr_detection_factory, subparser, script_name, version):
        """
        Creates a new Search sub-command instance.
        :param amr_detection_factory: A factory of type staramr.detection.AMRDetectionFactory for building necessary objects for AMR detection.
        :param subparser: The subparser to use.  Generated from argparse.ArgumentParser.add_subparsers().
        :param script_name: The name of the script being run.
        :param version: The version of this software.
        """
        super().__init__(subparser, script_name)
        self._amr_detection_factory = amr_detection_factory
        self._version = version

    def _setup_args(self, arg_parser):
        name = self._script_name
        epilog = ("Example:\n"
                  "\t" + name + " search --output-dir out *.fasta\n"
                                "\t\tSearches the files *.fasta for AMR genes using only the ResFinder database, storing results in the out/ directory.\n\n" +
                  "\t" + name + " search --pointfinder-organism salmonella --output-dir out *.fasta\n" +
                  "\t\tSearches *.fasta for AMR genes using ResFinder and PointFinder database with the passed organism, storing results in out/.")

        arg_parser = self._subparser.add_parser('search',
                                                epilog=epilog,
                                                formatter_class=argparse.RawTextHelpFormatter,
                                                help='Search for AMR genes')

        self._default_database_dir = AMRDatabaseHandlerFactory.get_default_database_directory()
        cpu_count = multiprocessing.cpu_count()

        arg_parser.add_argument('-n', '--nprocs', action='store', dest='nprocs', type=int,
                                help='The number of processing cores to use [' + str(cpu_count) + '].',
                                default=cpu_count, required=False)
        arg_parser.add_argument('--pid-threshold', action='store', dest='pid_threshold', type=float,
                                help='The percent identity threshold [98.0].', default=98.0, required=False)
        arg_parser.add_argument('--percent-length-overlap-resfinder', action='store',
                                dest='plength_threshold_resfinder', type=float,
                                help='The percent length overlap for resfinder results [60.0].', default=60.0,
                                required=False)
        arg_parser.add_argument('--percent-length-overlap-pointfinder', action='store',
                                dest='plength_threshold_pointfinder', type=float,
                                help='The percent length overlap for pointfinder results [95.0].', default=95.0,
                                required=False)
        arg_parser.add_argument('--pointfinder-organism', action='store', dest='pointfinder_organism', type=str,
                                help='The organism to use for pointfinder {' + ', '.join(
                                    PointfinderBlastDatabase.get_available_organisms()) + '} [None].', default=None,
                                required=False)
        arg_parser.add_argument('--include-negatives', action='store_true', dest='include_negatives',
                                help='Inclue negative results (those sensitive to antimicrobials) [False].',
                                required=False)
        arg_parser.add_argument('--report-all-blast', action='store_true', dest='report_all_blast',
                                help='Report all blast hits (vs. only top blast hits) [False].',
                                required=False)
        arg_parser.add_argument('-d', '--database', action='store', dest='database', type=str,
                                help='The directory containing the resfinder/pointfinder databases [' + self._default_database_dir + '].',
                                default=self._default_database_dir, required=False)
        arg_parser.add_argument('-o', '--output-dir', action='store', dest='output_dir', type=str,
                                help="The output directory for results.  If unset prints all results to stdout.",
                                default=None, required=False)
        arg_parser.add_argument('--version', action='store_true', dest='version',
                                help='Prints version information.', required=False)
        arg_parser.add_argument('files', nargs=argparse.REMAINDER)

        return arg_parser

    def _print_dataframes_to_excel(self, outfile_path, summary_dataframe, resfinder_dataframe, pointfinder_dataframe,
                                   settings_dataframe):
        writer = pandas.ExcelWriter(outfile_path, engine='xlsxwriter')

        summary_dataframe.to_excel(writer, 'Summary', freeze_panes=[1, 1], na_rep=self.blank)
        resfinder_dataframe.to_excel(writer, 'ResFinder', freeze_panes=[1, 1], na_rep=self.blank)
        if pointfinder_dataframe is not None:
            pointfinder_dataframe.to_excel(writer, 'PointFinder', freeze_panes=[1, 1], na_rep=self.blank)
        settings_dataframe.to_excel(writer, 'Settings')

        writer.save()

    def _print_dataframe_to_text_file(self, dataframe, file=None):
        file_handle = sys.stdout

        if dataframe is not None:
            if file:
                file_handle = open(file, 'w')

            dataframe.to_csv(file_handle, sep="\t", float_format="%0.2f", na_rep=self.blank)

            if file:
                file_handle.close()

    def _print_settings_to_file(self, settings, file):
        file_handle = open(file, 'w')
        file_handle.write(get_string_with_spacing(settings))
        file_handle.close()

    def run(self, args):
        super(Search, self).run(args)

        start_time = datetime.datetime.now()

        if (len(args.files) == 0):
            raise CommandParseException("Must pass a fasta file to process", self._root_arg_parser)

        hits_output_dir = None
        if args.output_dir:
            if path.exists(args.output_dir):
                raise CommandParseException("Output directory [" + args.output_dir + "] already exists",
                                            self._root_arg_parser)
            else:
                hits_output_dir = path.join(args.output_dir, 'hits')
                mkdir(args.output_dir)
                mkdir(hits_output_dir)

        if not path.isdir(args.database):
            raise CommandParseException("Database directory [" + args.database + "] does not exist")

        if args.database == AMRDatabaseHandlerFactory.get_default_database_directory():
            database_handler = AMRDatabaseHandlerFactory.create_default_factory().get_database_handler()
        else:
            database_handler = AMRDatabaseHandlerFactory(args.database)

        resfinder_database_dir = database_handler.get_resfinder_dir()
        pointfinder_database_dir = database_handler.get_pointfinder_dir()

        resfinder_database = ResfinderBlastDatabase(resfinder_database_dir)
        if (args.pointfinder_organism):
            if args.pointfinder_organism not in PointfinderBlastDatabase.get_available_organisms():
                raise CommandParseException("The only Pointfinder organism(s) currently supported are " + str(
                    PointfinderBlastDatabase.get_available_organisms()), self._root_arg_parser)
            pointfinder_database = PointfinderBlastDatabase(pointfinder_database_dir,
                                                            args.pointfinder_organism)
        else:
            pointfinder_database = None
        blast_handler = BlastHandler(resfinder_database, args.nprocs, pointfinder_database)

        amr_detection = self._amr_detection_factory.build(resfinder_database, blast_handler, pointfinder_database,
                                                          args.include_negatives, output_dir=hits_output_dir)
        amr_detection.run_amr_detection(args.files, args.pid_threshold, args.plength_threshold_resfinder,
                                        args.plength_threshold_pointfinder, args.report_all_blast)

        end_time = datetime.datetime.now()
        time_difference = end_time - start_time
        time_difference_minutes = "%0.2f" % (time_difference.total_seconds() / 60)

        logger.info("Finished. Took " + str(time_difference_minutes) + " minutes.")

        if args.output_dir:
            self._print_dataframe_to_text_file(amr_detection.get_resfinder_results(),
                                               path.join(args.output_dir, "resfinder.tsv"))
            self._print_dataframe_to_text_file(amr_detection.get_pointfinder_results(),
                                               path.join(args.output_dir, "pointfinder.tsv"))
            self._print_dataframe_to_text_file(amr_detection.get_summary_results(),
                                               path.join(args.output_dir, "summary.tsv"))

            settings = database_handler.info()
            settings.insert(0, ['command_line', ' '.join(sys.argv)])
            settings.insert(1, ['version', self._version])
            settings.insert(2, ['start_time', start_time.strftime("%Y-%m-%d %H:%M:%S")])
            settings.insert(3, ['end_time', end_time.strftime("%Y-%m-%d %H:%M:%S")])
            settings.insert(4, ['total_minutes', time_difference_minutes])
            self._print_settings_to_file(settings, path.join(args.output_dir, "settings.txt"))

            settings_dataframe = pandas.DataFrame(settings, columns=('Key', 'Value')).set_index('Key')

            self._print_dataframes_to_excel(path.join(args.output_dir, 'results.xlsx'),
                                            amr_detection.get_summary_results(), amr_detection.get_resfinder_results(),
                                            amr_detection.get_pointfinder_results(),
                                            settings_dataframe)

            logger.info("Output files in " + args.output_dir)
        else:
            self._print_dataframe_to_text_file(amr_detection.get_summary_results())
