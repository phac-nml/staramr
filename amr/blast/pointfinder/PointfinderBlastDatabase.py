from os import path, listdir

from amr.blast.AbstractBlastDatabase import AbstractBlastDatabase
from amr.blast.pointfinder.PointfinderDatabaseInfo import PointfinderDatabaseInfo

class PointfinderBlastDatabase(AbstractBlastDatabase):
    def __init__(self, database_dir, organism):
        super().__init__(database_dir)
        self.organism = organism
        self.pointfinder_database_dir = path.join(self.database_dir, organism)

        if (not path.isdir(self.pointfinder_database_dir)):
            raise Exception(
                "Error, pointfinder organism [" + organism + "] is either incorrect or pointfinder database not installed properly")

        self._pointfinder_info = PointfinderDatabaseInfo.from_file(path.join(self.pointfinder_database_dir, "resistens-overview.txt"))

    def get_database_names(self):
        return [f[:-len(self.fasta_suffix)] for f in listdir(self.pointfinder_database_dir) if
                (path.isfile(path.join(self.pointfinder_database_dir, f)) and f.endswith(self.fasta_suffix))]

    def get_path(self, database_name):
        return path.join(self.pointfinder_database_dir, database_name + self.fasta_suffix)

    def get_resistance_codons(self, gene, codon_mutations):
        return self._pointfinder_info.get_resistance_codons(gene, codon_mutations)

    def get_phenotype(self, gene, codon_mutation):
        return self._pointfinder_info.get_phenotype(gene, codon_mutation)
