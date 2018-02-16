from os import path, listdir
import pandas

from amr.blast.AbstractBlastDatabase import AbstractBlastDatabase


class PointfinderBlastDatabase(AbstractBlastDatabase):
    def __init__(self, database_dir, organism):
        super().__init__(database_dir)
        self.organism = organism
        self.pointfinder_database_dir = path.join(self.database_dir, organism)

        if (not path.isdir(self.pointfinder_database_dir)):
            raise Exception(
                "Error, pointfinder organism [" + organism + "] is either incorrect or pointfinder database not installed properly")

        self._parse_pointfinder_info(path.join(self.pointfinder_database_dir, "resistens-overview.txt"))

    def _parse_pointfinder_info(self, database_info_file):
        self._pointfinder_info = pandas.read_csv(database_info_file, sep=':', index_col=False)

    def get_database_names(self):
        return [f[:-len(self.fasta_suffix)] for f in listdir(self.pointfinder_database_dir) if
                (path.isfile(path.join(self.pointfinder_database_dir, f)) and f.endswith(self.fasta_suffix))]

    def get_path(self, database_name):
        return path.join(self.pointfinder_database_dir, database_name + self.fasta_suffix)
