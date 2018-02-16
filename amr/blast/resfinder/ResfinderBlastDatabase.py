import os

from amr.blast.AbstractBlastDatabase import AbstractBlastDatabase


class ResfinderBlastDatabase(AbstractBlastDatabase):

    def __init__(self, database_dir):
        super().__init__(database_dir)

    def get_database_names(self):
        return [f[:-len(self.fasta_suffix)] for f in os.listdir(self.database_dir) if
                (os.path.isfile(os.path.join(self.database_dir, f)) and f.endswith(self.fasta_suffix))]

    def get_path(self, database_name):
        return os.path.join(self.database_dir, database_name + self.fasta_suffix)