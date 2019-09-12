import logging
import os
import re
import subprocess
import math
from concurrent.futures import ThreadPoolExecutor
from os import path
from typing import Dict, List

from Bio.Blast.Applications import NcbiblastnCommandline

from staramr.blast.AbstractBlastDatabase import AbstractBlastDatabase
from staramr.exceptions.BlastProcessError import BlastProcessError

logger = logging.getLogger('JobHandler')

"""
Class for handling scheduling of BLAST jobs.
"""


class JobHandler:
    BLAST_COLUMNS = [x.strip() for x in '''
    qseqid
    sseqid
    pident
    length
    qstart
    qend
    sstart
    send
    slen
    qlen
    sstrand
    sseq
    qseq
    '''.strip().split('\n')]

    def __init__(self, blast_database_objects_map: Dict[str, AbstractBlastDatabase], threads: int,
                 output_directory: str) -> None:
        """
        Creates a new JobHandler.
        :param blast_database_objects_map: A map containing the blast databases.
        :param threads: The maximum number of threads to use, where one BLAST process gets assigned to one thread.
        :param output_directory: The output directory to store BLAST results.
        """
        if threads is None:
            raise Exception("threads is None")

        self._threads = threads
        self._mlst_version = None

        if output_directory is None:
            raise Exception("output_directory is None")

        self._output_directory = output_directory
        self._input_genomes_tmp_dir = path.join(output_directory, 'input-genomes')

        self._blast_database_objects_map = blast_database_objects_map

        if (self._blast_database_objects_map['pointfinder'] is None):
            self._pointfinder_configured = False  # type: bool
            del self._blast_database_objects_map['pointfinder']
        else:
            self._pointfinder_configured = True  # type: bool

        self._thread_pool_executor = ThreadPoolExecutor(max_workers=self._threads)
        self._max_mlst_columns = 10

        self.reset()

    def reset(self):
        """
        Resets this JobHandler.
        :return: None
        """
        if self._thread_pool_executor is not None:
            self._thread_pool_executor.shutdown()
        self._thread_pool_executor = ThreadPoolExecutor(max_workers=self._threads)
        self._blast_map = {}
        self._future_blasts_map = {}
        self._mlst_data = ""
        self._mlst_map = {}

        if path.exists(self._input_genomes_tmp_dir):
            logger.debug("Directory [%s] already exists", self._input_genomes_tmp_dir)
        else:
            os.mkdir(self._input_genomes_tmp_dir)

    def run_blasts_mlst(self, files, mlst_scheme) -> None:
        """
        Scans all files with BLAST against the ResFinder/PointFinder/Plasmid databases and scans all files with MLST
        :param files: The files to scan.
        :param mlst_scheme: Specifies scheme name for MLST to use.
        :return: None
        """
        db_files = self._make_db_from_input_files(self._input_genomes_tmp_dir, files)
        logger.debug("Done making blast databases for input files")

        future_mlst_db = [] # type: list

        for file in db_files:

            logger.info("Scheduling blasts and MLST for %s", path.basename(file))
            future_mlst_db.append(self._thread_pool_executor.submit(self._schedule_mlst, file, mlst_scheme))

            for name in self._blast_database_objects_map:
                database_object = self._blast_database_objects_map[name]
                self._schedule_blast(file, database_object)

        try:
            for future_mlst in future_mlst_db:
                mlst_result = future_mlst.result()

                self._mlst_data += mlst_result

        except subprocess.CalledProcessError as e:
            err_msg = str(e.stderr.strip())
            raise Exception('Could not run mlst, error {}'.format(err_msg))

    def _make_db_from_input_files(self, db_dir, files):
        logger.info("Making BLAST databases for input files")
        future_makeblastdbs = []
        db_files = []

        for file in files:
            destination = path.join(db_dir, path.basename(file))
            logger.debug("Creating symlink from [%s] to [%s]", file, destination)
            os.symlink(path.abspath(file), destination)
            db_files.append(destination)

            future_makeblastdbs.append(self._thread_pool_executor.submit(self._make_blast_db, destination, file))

        # Blocks until all blast dbs are made. If an exception is raised, will raise same exception
        try:
            for future_blastdb in future_makeblastdbs:
                future_blastdb.result()
        except subprocess.CalledProcessError as e:
            raise BlastProcessError("Error running makeblastdb", e)

        return db_files

    def _schedule_mlst(self, file: str, mlst_scheme: str) -> str:

        command = ['mlst']

        if mlst_scheme is not None:
            command.extend(['--scheme', mlst_scheme])

        command.append(file);

        logger.debug(' '.join(command))
        try:
            output = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

            decoded_output = str(output.stdout, 'utf-8')

        except subprocess.CalledProcessError as e:
            err_msg = str(e.stderr.strip())

            raise Exception('Could not run mlst, error {}'.format(err_msg))

        return decoded_output

    def _schedule_blast(self, file, blast_database):
        database_names = blast_database.get_database_names()
        logger.debug("%s databases: %s", blast_database.get_name(), database_names)
        for database_name in database_names:
            database = blast_database.get_path(database_name)
            file_name = os.path.basename(file)

            blast_out = os.path.join(self._output_directory,
                                     file_name + "." + database_name + "." + blast_database.get_name() + ".blast.tsv")
            if os.path.exists(blast_out):
                raise Exception("Error, blast_out [%s] already exists", blast_out)

            self._get_blast_map(blast_database.get_name()).setdefault(file_name, {})[database_name] = blast_out

            future_blast = self._thread_pool_executor.submit(self._launch_blast, database, file, blast_out)
            self._get_future_blasts_from_map(blast_database.get_name()).append(future_blast)

    def _get_blast_map(self, name: str) -> Dict:
        if name not in self._blast_map:
            self._blast_map[name] = {}

        return self._blast_map[name]

    def _get_mlst_data(self) -> str:

        return self._mlst_data

    def _get_future_blasts_from_map(self, name: str) -> Dict:
        if name not in self._future_blasts_map:
            self._future_blasts_map[name] = []

        return self._future_blasts_map[name]

    def _get_mlst_version(self) -> str:
      command = ['mlst', '--version']

      try:
          output = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)

          mlst_version = str(output.stdout, 'utf-8')

          # Parses out the mlst when the string is given back ex `mlst 2.x.x` and removes new line
          mlst_version = (mlst_version[5:]).rstrip()

      except subprocess.CalledProcessError as e:
          err_msg = str(e.stderr.strip())

          raise Exception('Could not run mlst, error {}'.format(err_msg))

      return mlst_version

    def is_pointfinder_configured(self) -> bool:
        """
        Whether or not PointFinder is being used.
        :return: True if PointFinder is being used, False otherwise.
        """
        return self._pointfinder_configured

    def get_resfinder_outputs(self) -> Dict:
        """
        Gets the ResFinder output files in the form of a dictionary which looks like:
            { 'input_file_name' => 'blast_results_file.xml' }
        :return: A dictionary mapping input file names to ResFinder BLAST output files.
        """

        # Forces any exceptions to be thrown if error with blasts
        for future_blast in self._get_future_blasts_from_map('resfinder'):
            future_blast.result()
        return self._get_blast_map('resfinder')

    def get_plasmidfinder_outputs(self) -> Dict:
        """
        Gets the PlasmidFinder output files in the form of a dictionary which looks like:
            { 'input_file_name' => 'blast_results_file.xml' }
        :return: A dictionary mapping input file names to ResFinder BLAST output files.
        """

        # Forces any exceptions to be thrown if error with blasts
        for future_blast in self._get_future_blasts_from_map('plasmidfinder'):
            future_blast.result()
        return self._get_blast_map('plasmidfinder')

    def get_mlst_outputs(self) -> str:
        """
        Gets the MLST output files from the MLST subprocess
        :return A decoded parsed list that contains all of the found locus in each file
        """

        return self._get_mlst_data();

    def get_mlst_version(self) -> str:

      return self._get_mlst_version(self) # type: ignore

    def get_pointfinder_outputs(self) -> Dict:
        """
        Gets the PointFinder output files in the form of a dictionary which looks like:
            { 'input_file_name' => 'blast_results_file.xml' }
        :return: A dictionary mapping input file names to PointFinder BLAST output files.
        """
        if (self.is_pointfinder_configured()):
            # Forces any exceptions to be thrown if error with blasts
            for future_blast in self._get_future_blasts_from_map('pointfinder'):
                future_blast.result()
            return self._get_blast_map('pointfinder')
        else:
            raise Exception('Error, pointfinder has not been configured')

    def _launch_blast(self, query, db, output) -> None:
        blast_out_format = '"6 ' + ' '.join(self.BLAST_COLUMNS) + '"'

        blastn_command = NcbiblastnCommandline(query=query, db=db, evalue=0.001, outfmt=blast_out_format, out=output)
        stdout, stderr = blastn_command()

        if stderr:
            raise Exception("error with [" + str(blastn_command) + "], stderr=" + stderr)

    def _make_blast_db(self, path: str, file: str) -> None:
        command = ['makeblastdb', '-in', path, '-dbtype', 'nucl', '-parse_seqids']
        logger.debug(' '.join(command))
        try:
            subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
        except subprocess.CalledProcessError as e:
            err_msg = str(e.stderr.strip())
            err_msg = re.findall('REF\|(.*?)\'', err_msg)[0]
            raise Exception('Could not run makeblastdb on file {}, error {}'.format(file, err_msg))
