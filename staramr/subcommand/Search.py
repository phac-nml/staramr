import argparse
import datetime
import logging
import multiprocessing
import sys
import tempfile
from os import path, mkdir

import numpy as np
import pandas as pd

from staramr.SubCommand import SubCommand
from staramr.Utils import get_string_with_spacing
from staramr.blast.BlastHandler import BlastHandler
from staramr.blast.plasmidfinder.PlasmidfinderBlastDatabase import PlasmidfinderBlastDatabase
from staramr.blast.pointfinder.PointfinderBlastDatabase import PointfinderBlastDatabase
from staramr.databases.AMRDatabasesManager import AMRDatabasesManager
from staramr.databases.exclude.ExcludeGenesList import ExcludeGenesList
from staramr.databases.resistance.ARGDrugTable import ARGDrugTable
from staramr.detection.AMRDetectionFactory import AMRDetectionFactory
from staramr.exceptions.CommandParseException import CommandParseException

logger = logging.getLogger("Search")

"""
Class for searching for AMR resistance genes.
"""


class Search(SubCommand):
    BLANK = '-'
    TIME_FORMAT = "%Y-%m-%d %H:%M:%S"

    def __init__(self, subparser, script_name, version):
        """
        Creates a new Search sub-command instance.
        :param subparser: The subparser to use.  Generated from argparse.ArgumentParser.add_subparsers().
        :param script_name: The name of the script being run.
        :param version: The version of this software.
        """
        super().__init__(subparser, script_name)
        self._version = version

    def _setup_args(self, arg_parser):
        name = self._script_name
        epilog = ("Example:\n"
                  "\t" + name + " search -o out *.fasta\n"
                                "\t\tSearches the files *.fasta for AMR genes using only the ResFinder database, storing results in the out/ directory.\n\n" +
                  "\t" + name + " search --pointfinder-organism salmonella --output-excel results.xlsx *.fasta\n" +
                  "\t\tSearches *.fasta for AMR genes using ResFinder and PointFinder database with the passed organism, storing results in results.xlsx.")

        arg_parser = self._subparser.add_parser('search',
                                                epilog=epilog,
                                                formatter_class=argparse.RawTextHelpFormatter,
                                                help='Search for AMR genes')

        self._default_database_dir = AMRDatabasesManager.get_default_database_directory()
        cpu_count = multiprocessing.cpu_count()

        arg_parser.add_argument('--pointfinder-organism', action='store', dest='pointfinder_organism', type=str,
                                help='The organism to use for pointfinder {' + ', '.join(
                                    PointfinderBlastDatabase.get_available_organisms()) + '}. Defaults to disabling search for point mutations. [None].',
                                default=None, required=False)
        arg_parser.add_argument('--plasmidfinder-database-type', action='store', dest='plasmidfinder_database_type',
                                type=str,
                                help='The database type to use for plasmidfinder {' + ', '.join(
                                    PlasmidfinderBlastDatabase.get_available_databases()) + '}. Defaults to using all available database types to search for plasmids. [None].',
                                default=None, required=False)
        arg_parser.add_argument('-d', '--database', action='store', dest='database', type=str,
                                help='The directory containing the resfinder/pointfinder/plasmidfinder databases [' + self._default_database_dir + '].',
                                default=self._default_database_dir, required=False)
        arg_parser.add_argument('-n', '--nprocs', action='store', dest='nprocs', type=int,
                                help='The number of processing cores to use [' + str(cpu_count) + '].',
                                default=cpu_count, required=False)
        arg_parser.add_argument('--ignore-invalid-files', action='store_true', dest='ignore_valid_files',
                                help='Skips over invalid input files', required=False)

        threshold_group = arg_parser.add_argument_group('BLAST Thresholds')
        threshold_group.add_argument('--pid-threshold', action='store', dest='pid_threshold', type=float,
                                     help='The percent identity threshold [98.0].', default=98.0, required=False)
        threshold_group.add_argument('--percent-length-overlap-resfinder', action='store',
                                     dest='plength_threshold_resfinder', type=float,
                                     help='The percent length overlap for resfinder results [60.0].', default=60.0,
                                     required=False)
        threshold_group.add_argument('--percent-length-overlap-pointfinder', action='store',
                                     dest='plength_threshold_pointfinder', type=float,
                                     help='The percent length overlap for pointfinder results [95.0].', default=95.0,
                                     required=False)
        threshold_group.add_argument('--percent-length-overlap-plasmidfinder', action='store',
                                     dest='plength_threshold_plasmidfinder', type=float,
                                     help='The percent length overlap for resfinder results [60.0].', default=60.0,
                                     required=False)

        report_group = arg_parser.add_argument_group('Reporting options')
        report_group.add_argument('--no-exclude-genes', action='store_true', dest='no_exclude_genes',
                                  help='Disable the default exclusion of some genes from ResFinder/PointFinder/PlasmidFinder [False].',
                                  required=False)
        report_group.add_argument('--exclude-genes-file', action='store', dest='exclude_genes_file',
                                  help='A containing a list of ResFinder/PointFinder/PlasmidFinder gene names to exclude from results [{}].'.format(
                                      ExcludeGenesList.get_default_exclude_file()),
                                  default=ExcludeGenesList.get_default_exclude_file(),
                                  required=False)
        report_group.add_argument('--exclude-negatives', action='store_true', dest='exclude_negatives',
                                  help='Exclude negative results (those sensitive to antimicrobials) [False].',
                                  required=False)
        report_group.add_argument('--exclude-resistance-phenotypes', action='store_true',
                                  dest='exclude_resistance_phenotypes',
                                  help='Exclude predicted antimicrobial resistances [False].',
                                  required=False)
        report_group.add_argument('--report-all-blast', action='store_true', dest='report_all_blast',
                                  help='Report all blast hits (vs. only top blast hits) [False].',
                                  required=False)

        output_group = arg_parser.add_argument_group(title='Output',
                                                     description='Use either --output-dir or specify individual output files')
        output_group.add_argument('-o', '--output-dir', action='store', dest='output_dir', type=str,
                                  help="The output directory for results [None].",
                                  default=None, required=False)
        output_group.add_argument('--output-summary', action='store', dest='output_summary', type=str,
                                  help="The name of the output file containing the summary results. Not be be used with '--output-dir'. [None]",
                                  default=None, required=False)
        output_group.add_argument('--output-detailed-summary', action='store', dest='output_detailed_summary', type=str,
                                  help="The name of the output file containing the detailed summary results. Not be be used with '--output-dir'. [None]",
                                  default=None, required=False)
        output_group.add_argument('--output-resfinder', action='store', dest='output_resfinder', type=str,
                                  help="The name of the output file containing the resfinder results. Not be be used with '--output-dir'. [None]",
                                  default=None, required=False)
        output_group.add_argument('--output-pointfinder', action='store', dest='output_pointfinder', type=str,
                                  help="The name of the output file containing the pointfinder results. Not be be used with '--output-dir'. [None]",
                                  default=None, required=False)
        output_group.add_argument('--output-plasmidfinder', action='store', dest='output_plasmidfinder', type=str,
                                  help="The name of the output file containing the plasmidfinder results. Not be be used with '--output-dir'. [None]",
                                  default=None, required=False)
        output_group.add_argument('--output-settings', action='store', dest='output_settings', type=str,
                                  help="The name of the output file containing the settings. Not be be used with '--output-dir'. [None]",
                                  default=None, required=False)
        output_group.add_argument('--output-excel', action='store', dest='output_excel', type=str,
                                  help="The name of the output file containing the excel results. Not be be used with '--output-dir'. [None]",
                                  default=None, required=False)
        output_group.add_argument('--output-hits-dir', action='store', dest='hits_output_dir', type=str,
                                  help="The name of the directory to contain the BLAST hit files. Not be be used with '--output-dir'. [None]",
                                  default=None, required=False)

        arg_parser.add_argument('files', nargs='+')

        return arg_parser

    def _print_dataframes_to_excel(self, outfile_path, summary_dataframe, resfinder_dataframe, pointfinder_dataframe,
                                   plasmidfinder_dataframe, detailed_summary_dataframe,
                                   settings_dataframe):
        writer = pd.ExcelWriter(outfile_path, engine='xlsxwriter')

        sheetname_dataframe = {}
        sheetname_dataframe['Summary'] = summary_dataframe
        sheetname_dataframe['Detailed_Summary'] = detailed_summary_dataframe
        sheetname_dataframe['ResFinder'] = resfinder_dataframe
        sheetname_dataframe['PlasmidFinder'] = plasmidfinder_dataframe
        if pointfinder_dataframe is not None:
            sheetname_dataframe['PointFinder'] = pointfinder_dataframe

        for name in ['Summary', 'Detailed_Summary', 'ResFinder', 'PointFinder', 'PlasmidFinder']:
            if name in sheetname_dataframe:
                sheetname_dataframe[name].to_excel(writer, name, freeze_panes=[1, 1], float_format="%0.2f",
                                                   na_rep=self.BLANK)
        self._resize_columns(sheetname_dataframe, writer, max_width=50)

        settings_dataframe.to_excel(writer, 'Settings')
        self._resize_columns({'Settings': settings_dataframe}, writer, max_width=75, text_wrap=False)

        writer.save()

    def _resize_columns(self, sheetname_dataframe, writer, max_width, text_wrap=True):
        """
        Resizes columns in workbook.
        :param sheetname_dataframe: A map mapping the sheet name to a dataframe.
        :param writer: The ExcelWriter, which the worksheets already added using writer.to_excel
        :param max_width: The maximum width of the columns.
        :param text_wrap: Whether or not to turn on text wrapping if columns surpass max_width.
        :return: None
        """
        workbook = writer.book
        wrap_format = workbook.add_format({'text_wrap': text_wrap})
        for name in sheetname_dataframe:
            for i, width in enumerate(self._get_col_widths(sheetname_dataframe[name])):
                if width > max_width:
                    writer.sheets[name].set_column(i, i, width=max_width, cell_format=wrap_format)
                else:
                    writer.sheets[name].set_column(i, i, width=width)

    def _get_col_widths(self, df):
        """
        Calculate column widths based on column headers and contents
        :param df: The dataframe.
        :return: A generator giving the max width for each column.
        """
        idx_max = max([len(str(s)) for s in df.index.values] + [len(str(df.index.name))])
        yield idx_max

        extra = 2
        for c in df.columns:
            # get max length of column contents and length of column header (plus some extra)
            yield np.max([df[c].astype(str).str.len().max(), len(c)]) + extra

    def _print_dataframe_to_text_file_handle(self, dataframe, file_handle):
        dataframe.to_csv(file_handle, sep="\t", float_format="%0.2f", na_rep=self.BLANK)

    def _print_settings_to_file(self, settings, file):
        file_handle = open(file, 'w')
        file_handle.write(get_string_with_spacing(settings))
        file_handle.close()

    def _generate_results(self, database_repos, resfinder_database, pointfinder_database, plasmidfinder_database,
                          nprocs, include_negatives,
                          include_resistances, hits_output, pid_threshold, plength_threshold_resfinder,
                          plength_threshold_pointfinder, plength_threshold_plasmidfinder, report_all_blast,
                          genes_to_exclude, files, ignore_invalid_files):
        """
        Runs AMR detection and generates results.
        :param database_repos: The database repos object.
        :param resfinder_database: The resfinder database.
        :param pointfinder_database: The pointfinder database.
        :param plasmidfinder_database: The plasmidfinder database.
        :param nprocs: The number of processing cores to use for BLAST.
        :param include_negatives: Whether or not to include negative results in output.
        :param include_resistances: Whether or not to include resistance phenotypes in output.
        :param hits_output: Output directory for hit files.
        :param pid_threshold: The pid threshold.
        :param plength_threshold_resfinder: The plength threshold for resfinder.
        :param plength_threshold_pointfinder: The plength threshold for pointfinder.
        :param plength_threshold_plasmidfinder: The plength threshold for plasmidfinder.
        :param report_all_blast: Whether or not to report all BLAST results.
        :param genes_to_exclude: A list of gene IDs to exclude from the results.
        :param files: The list of files to scan.
        :param ignore_invalid_files: Skips over invalid input files
        :return: A dictionary containing the results as dict['results'] and settings as dict['settings'].
        """
        results = {'results': None, 'settings': None}

        with tempfile.TemporaryDirectory() as blast_out:
            start_time = datetime.datetime.now()

            blast_handler = BlastHandler({'resfinder': resfinder_database, 'pointfinder': pointfinder_database,
                                          'plasmidfinder': plasmidfinder_database}, nprocs, blast_out)

            amr_detection_factory = AMRDetectionFactory()
            amr_detection = amr_detection_factory.build(plasmidfinder_database,
                                                        resfinder_database,
                                                        blast_handler,
                                                        pointfinder_database,
                                                        include_negatives=include_negatives,
                                                        include_resistances=include_resistances,
                                                        output_dir=hits_output,
                                                        genes_to_exclude=genes_to_exclude)
            amr_detection.run_amr_detection(files, pid_threshold, plength_threshold_resfinder,
                                            plength_threshold_pointfinder, plength_threshold_plasmidfinder,
                                            report_all_blast, ignore_invalid_files)

            results['results'] = amr_detection

            end_time = datetime.datetime.now()
            time_difference = end_time - start_time
            time_difference_minutes = "%0.2f" % (time_difference.total_seconds() / 60)

            logger.info("Finished. Took %s minutes.", time_difference_minutes)

            settings = database_repos.info()
            settings['command_line'] = ' '.join(sys.argv)
            settings['version'] = self._version
            settings['start_time'] = start_time.strftime(self.TIME_FORMAT)
            settings['end_time'] = end_time.strftime(self.TIME_FORMAT)
            settings['total_minutes'] = time_difference_minutes
            settings.move_to_end('total_minutes', last=False)
            settings.move_to_end('end_time', last=False)
            settings.move_to_end('start_time', last=False)
            settings.move_to_end('version', last=False)
            settings.move_to_end('command_line', last=False)

            if include_resistances:
                arg_drug_table = ARGDrugTable()
                info = arg_drug_table.get_resistance_table_info()
                settings.update(info)
                logger.info(
                    "Predicting AMR resistance phenotypes is enabled. The predictions are for microbiological " +
                    "resistance and *not* clinical resistance. These results are continually being improved and " +
                    "we welcome any feedback.")

            results['settings'] = settings

        return results

    def run(self, args):
        super(Search, self).run(args)

        if (len(args.files) == 0):
            raise CommandParseException("Must pass a fasta file to process", self._root_arg_parser, print_help=True)

        for file in args.files:
            if not path.exists(file):
                raise CommandParseException('File [' + file + '] does not exist', self._root_arg_parser)

        if not path.isdir(args.database):
            if args.database == self._default_database_dir:
                raise CommandParseException(
                    "Default database does not exist. Perhaps try restoring with 'staramr db restore-default'",
                    self._root_arg_parser)
            else:
                raise CommandParseException(
                    "Database directory [" + args.database + "] does not exist. Perhaps try building with" +
                    "'staramr db build --dir " + args.database + "'",
                    self._root_arg_parser)

        if args.database == AMRDatabasesManager.get_default_database_directory():
            database_repos = AMRDatabasesManager.create_default_manager().get_database_repos()
        else:
            database_repos = AMRDatabasesManager(args.database).get_database_repos()

        if not AMRDatabasesManager.is_database_repos_default_commits(database_repos):
            logger.warning("Using non-default ResFinder/PointFinder. This may lead to differences in the detected " +
                           "AMR genes depending on how the database files are structured.")

        resfinder_database = database_repos.build_blast_database('resfinder')
        if (args.pointfinder_organism):
            if args.pointfinder_organism not in PointfinderBlastDatabase.get_available_organisms():
                raise CommandParseException("The only Pointfinder organism(s) currently supported are " + str(
                    PointfinderBlastDatabase.get_available_organisms()), self._root_arg_parser)
            pointfinder_database = database_repos.build_blast_database('pointfinder',
                                                                       {'organism': args.pointfinder_organism})
        else:
            logger.info("No --pointfinder-organism specified. Will not search the PointFinder databases")
            pointfinder_database = None

        if (args.plasmidfinder_database_type):
            if args.plasmidfinder_database_type not in PlasmidfinderBlastDatabase.get_available_databases():
                raise CommandParseException("The only Plasmidfinder databases that are currently supported are " + str(
                    PlasmidfinderBlastDatabase.get_available_databases()), self._root_arg_parser)
            plasmidfinder_database = database_repos.build_blast_database('plasmidfinder', {
                'database_type': args.plasmidfinder_database_type})
        else:
            logger.info("No --plasmidfinder-database-type specified. Will search the entire PlasmidFinder database")
            plasmidfinder_database = database_repos.build_blast_database('plasmidfinder')

        hits_output_dir = None
        output_summary = None
        output_detailed_summary = None
        output_resfinder = None
        output_pointfinder = None
        output_plasmidfinder = None
        output_excel = None
        output_settings = None
        if args.output_dir:
            if path.exists(args.output_dir):
                raise CommandParseException("Output directory [" + args.output_dir + "] already exists",
                                            self._root_arg_parser)
            elif args.output_summary or args.output_detailed_summary or args.output_resfinder or args.output_pointfinder or args.output_plasmidfinder or args.output_excel or \
                    args.hits_output_dir:
                raise CommandParseException('You cannot use --output-[type] with --output-dir', self._root_arg_parser)
            else:
                mkdir(args.output_dir)

                hits_output_dir = path.join(args.output_dir, 'hits')
                output_resfinder = path.join(args.output_dir, "resfinder.tsv")
                output_pointfinder = path.join(args.output_dir, "pointfinder.tsv")
                output_plasmidfinder = path.join(args.output_dir, "plasmidfinder.tsv")
                output_summary = path.join(args.output_dir, "summary.tsv")
                output_detailed_summary = path.join(args.output_dir, "detailed_summary.tsv")
                output_settings = path.join(args.output_dir, "settings.txt")
                output_excel = path.join(args.output_dir, 'results.xlsx')

                mkdir(hits_output_dir)

                logger.info("--output-dir set. All files will be output to [%s]", args.output_dir)
        elif args.output_summary or args.output_excel or args.output_detailed_summary:
            logger.info('--output-dir not set. Files will be output to the respective --output-[type] setting')
            output_resfinder = args.output_resfinder
            output_pointfinder = args.output_pointfinder
            output_plasmidfinder = args.output_plasmidfinder
            output_summary = args.output_summary
            output_detailed_summary = args.output_detailed_summary
            output_settings = args.output_settings
            output_excel = args.output_excel
            hits_output_dir = args.hits_output_dir

            if hits_output_dir is not None:
                if path.exists(hits_output_dir) and not path.isdir(hits_output_dir):
                    raise CommandParseException(
                        "--output-hits-dir [" + hits_output_dir + "] exists and is not a directory",
                        self._root_arg_parser)
                elif path.exists(hits_output_dir):
                    logger.debug("Found --output-hits-dir [%s] and is a directory. Will write hits here",
                                 hits_output_dir)
                else:
                    logger.debug("Making directory [%s]", hits_output_dir)
                    mkdir(hits_output_dir)
        else:
            raise CommandParseException('You must set one of --output-dir, --output-summary, --output-detailed-summary, or --output-excel',
                                        self._root_arg_parser)

        if args.no_exclude_genes:
            logger.info("--no-exclude-genes enabled. Will not exclude any ResFinder/PointFinder genes.")
            exclude_genes = []
        else:
            if not path.exists(args.exclude_genes_file):
                raise CommandParseException('--exclude-genes-file [{}] does not exist'.format(args.exclude_genes_file),
                                            self._root_arg_parser)
            else:
                logger.info(
                    "Will exclude ResFinder/PointFinder genes listed in [%s]. Use --no-exclude-genes to disable",
                    args.exclude_genes_file)
                exclude_genes = ExcludeGenesList(args.exclude_genes_file).tolist()

        results = self._generate_results(database_repos=database_repos,
                                         resfinder_database=resfinder_database,
                                         pointfinder_database=pointfinder_database,
                                         plasmidfinder_database=plasmidfinder_database,
                                         nprocs=args.nprocs,
                                         include_negatives=not args.exclude_negatives,
                                         include_resistances=not args.exclude_resistance_phenotypes,
                                         hits_output=hits_output_dir,
                                         pid_threshold=args.pid_threshold,
                                         plength_threshold_resfinder=args.plength_threshold_resfinder,
                                         plength_threshold_pointfinder=args.plength_threshold_pointfinder,
                                         plength_threshold_plasmidfinder=args.plength_threshold_plasmidfinder,
                                         report_all_blast=args.report_all_blast,
                                         genes_to_exclude=exclude_genes,
                                         files=args.files,
                                         ignore_invalid_files=args.ignore_valid_files)
        amr_detection = results['results']
        settings = results['settings']

        if output_resfinder:
            logger.info("Writing resfinder to [%s]", output_resfinder)
            with open(output_resfinder, 'w') as fh:
                self._print_dataframe_to_text_file_handle(amr_detection.get_resfinder_results(), fh)
        else:
            logger.info("--output-dir or --output-resfinder unset. No resfinder file will be written")

        if args.pointfinder_organism and output_pointfinder:
            logger.info("Writing pointfinder to [%s]", output_pointfinder)
            with open(output_pointfinder, 'w') as fh:
                self._print_dataframe_to_text_file_handle(amr_detection.get_pointfinder_results(), fh)
        else:
            logger.info("--output-dir or --output-pointfinder unset. No pointfinder file will be written")

        if output_plasmidfinder:
            logger.info("Writing plasmidfinder to [%s]", output_plasmidfinder)
            with open(output_plasmidfinder, 'w') as fh:
                self._print_dataframe_to_text_file_handle(amr_detection.get_plasmidfinder_results(), fh)
        else:
            logger.info("--output-dir or --output-plasmidfinder unset. No plasmidfinder file will be written")

        if output_summary:
            logger.info("Writing summary to [%s]", output_summary)
            with open(output_summary, 'w') as fh:
                self._print_dataframe_to_text_file_handle(amr_detection.get_summary_results(), fh)
        else:
            logger.info("--output-dir or --output-summary unset. No summary file will be written")

        if output_detailed_summary:
            logger.info("Writing detailed summary to [%s]", output_detailed_summary)
            with open(output_detailed_summary, 'w') as fh:
                self._print_dataframe_to_text_file_handle(amr_detection.get_detailed_summary_results(), fh)
        else:
            logger.info("--output-dir or --output-detailed-summary unset. No detailed summary file will be written")

        if output_settings:
            logger.info("Writing settings to [%s]", output_settings)
            self._print_settings_to_file(settings, output_settings)
        else:
            logger.info("--output-dir or --output-settings unset. No settings file will be written")

        if output_excel:
            logger.info("Writing Excel to [%s]", output_excel)
            settings_dataframe = pd.DataFrame.from_dict(settings, orient='index')
            settings_dataframe.index.name = 'Key'
            settings_dataframe.set_axis(['Value'], axis='columns', inplace=True)

            self._print_dataframes_to_excel(output_excel,
                                            amr_detection.get_summary_results(),
                                            amr_detection.get_resfinder_results(),
                                            amr_detection.get_pointfinder_results(),
                                            amr_detection.get_plasmidfinder_results(),
                                            amr_detection.get_detailed_summary_results(),
                                            settings_dataframe)
        else:
            logger.info("--output-dir or --output-excel unset. No excel file will be written")

        if hits_output_dir:
            logger.info("BLAST hits are stored in [%s]", hits_output_dir)
        else:
            logger.info("--output-dir or --output-hits-dir not set. No BLAST hits will be saved.")
