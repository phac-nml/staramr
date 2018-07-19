import logging
import os
import subprocess
from concurrent.futures import ThreadPoolExecutor
from os import path

from Bio.Blast.Applications import NcbiblastnCommandline

from staramr.exceptions.BlastProcessError import BlastProcessError

logger = logging.getLogger('BlastHandler')

"""
Class for handling scheduling of BLAST jobs.
"""


class BlastHandler:
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

    def __init__(self, resfinder_database, threads, output_directory, pointfinder_database=None):
        """
        Creates a new BlastHandler.
        :param resfinder_database: The staramr.blast.resfinder.ResfinderBlastDatabase for the particular ResFinder database.
        :param threads: The maximum number of threads to use, where one BLAST process gets assigned to one thread.
        :param output_directory: The output directory to store BLAST results.
        :param pointfinder_database: The staramr.blast.pointfinder.PointfinderBlastDatabase to use for the particular PointFinder database.
        """
        self._resfinder_database = resfinder_database

        if threads is None:
            raise Exception("threads is None")

        self._threads = threads

        if output_directory is None:
            raise Exception("output_directory is None")

        self._output_directory = output_directory
        self._input_genomes_tmp_dir = path.join(output_directory, 'input-genomes')

        if (pointfinder_database == None):
            self._pointfinder_configured = False
        else:
            self._pointfinder_database = pointfinder_database
            self._pointfinder_configured = True

        self._thread_pool_executor = None
        self.reset()

    def reset(self):
        """
        Resets this BlastHandler.
        :return: None
        """
        if self._thread_pool_executor is not None:
            self._thread_pool_executor.shutdown()
        self._thread_pool_executor = ThreadPoolExecutor(max_workers=self._threads)
        self._resfinder_blast_map = {}
        self._pointfinder_blast_map = {}
        self._pointfinder_future_blasts = []
        self._resfinder_future_blasts = []

        if path.exists(self._input_genomes_tmp_dir):
            logger.debug("Directory [%s] already exists", self._input_genomes_tmp_dir)
        else:
            os.mkdir(self._input_genomes_tmp_dir)

    def run_blasts(self, files):
        """
        Scans all files with BLAST against the ResFinder/PointFinder databases.
        :param files: The files to scan.
        :return: None
        """
        database_names_resfinder = self._resfinder_database.get_database_names()
        logger.debug("Resfinder Databases: %s", database_names_resfinder)

        if self.is_pointfinder_configured():
            database_names_pointfinder = self._pointfinder_database.get_database_names()
            logger.debug("Pointfinder Databases: %s", database_names_pointfinder)
        else:
            database_names_pointfinder = None

        db_files = self._make_db_from_input_files(self._input_genomes_tmp_dir, files)
        logger.debug("Done making blast databases for input files")

        for file in db_files:
            logger.info("Scheduling blast for %s", path.basename(file))
            self._schedule_resfinder_blast(file, database_names_resfinder)
            if self.is_pointfinder_configured():
                self._schedule_pointfinder_blast(file, database_names_pointfinder)

    def _make_db_from_input_files(self, db_dir, files):
        logger.info("Making BLAST databases for input files")
        future_makeblastdbs = []
        db_files = []

        for file in files:
            destination = path.join(db_dir, path.basename(file))
            logger.debug("Creating symlink from [%s] to [%s]", file, destination)
            os.symlink(path.abspath(file), destination)
            db_files.append(destination)

            future_makeblastdbs.append(self._thread_pool_executor.submit(self._make_blast_db, destination))

        # Blocks until all blast dbs are made. If an exception is raised, will raise same exception
        try:
            for future_blastdb in future_makeblastdbs:
                future_blastdb.result()
        except subprocess.CalledProcessError as e:
            raise BlastProcessError("Error running makeblastdb", e)

        return db_files

    def _schedule_resfinder_blast(self, file, database_names):
        for database_name in database_names:
            database = self._resfinder_database.get_path(database_name)
            file_name = os.path.basename(file)

            blast_out = os.path.join(self._output_directory, file_name + "." + database_name + ".resfinder.blast.xml")
            if os.path.exists(blast_out):
                raise Exception("Error, blast_out [%s] already exists", blast_out)

            self._resfinder_blast_map.setdefault(file_name, {})[database_name] = blast_out

            future_blast = self._thread_pool_executor.submit(self._launch_blast, database, file, blast_out)
            self._resfinder_future_blasts.append(future_blast)

    def _schedule_pointfinder_blast(self, file, database_names):
        for database_name in database_names:
            database = self._pointfinder_database.get_path(database_name)
            file_name = os.path.basename(file)

            blast_out = os.path.join(self._output_directory, file_name + "." + database_name + ".pointfinder.blast.xml")
            if os.path.exists(blast_out):
                raise Exception("Error, blast_out [%s] already exists", blast_out)

            self._pointfinder_blast_map.setdefault(file_name, {})[database_name] = blast_out

            future_blast = self._thread_pool_executor.submit(self._launch_blast, database, file, blast_out)
            self._pointfinder_future_blasts.append(future_blast)

    def is_pointfinder_configured(self):
        """
        Whether or not PointFinder is being used.
        :return: True if PointFinder is being used, False otherwise.
        """
        return self._pointfinder_configured

    def get_resfinder_outputs(self):
        """
        Gets the ResFinder output files in the form of a dictionary which looks like:
            { 'input_file_name' => 'blast_results_file.xml' }
        :return: A dictionary mapping input file names to ResFinder BLAST output files.
        """

        # Forces any exceptions to be thrown if error with blasts
        for future_blast in self._resfinder_future_blasts:
            future_blast.result()
        return self._resfinder_blast_map

    def get_pointfinder_outputs(self):
        """
        Gets the PointFinder output files in the form of a dictionary which looks like:
            { 'input_file_name' => 'blast_results_file.xml' }
        :return: A dictionary mapping input file names to PointFinder BLAST output files.
        """
        if (self.is_pointfinder_configured()):
            # Forces any exceptions to be thrown if error with blasts
            for future_blast in self._pointfinder_future_blasts:
                future_blast.result()
            return self._pointfinder_blast_map
        else:
            raise Exception("Error, pointfinder has not been configured")

    def _launch_blast(self, query, db, output):
        blast_out_format = '"6 ' + ' '.join(self.BLAST_COLUMNS) + '"'
        blastn_command = NcbiblastnCommandline(query=query, db=db, evalue=0.001, outfmt=blast_out_format, out=output)
        logger.debug(blastn_command)
        stdout, stderr = blastn_command()
        if stderr:
            raise Exception("error with [" + str(blastn_command) + "], stderr=" + stderr)

    def _make_blast_db(self, path):
        command = ['makeblastdb', '-in', path, '-dbtype', 'nucl', '-parse_seqids']
        logger.debug(' '.join(command))
        subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=True)
