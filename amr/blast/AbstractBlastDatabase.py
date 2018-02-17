import abc


class AbstractBlastDatabase:
    fasta_suffix = ".fsa"

    def __init__(self, database_dir):
        __metaclass__ = abc.ABCMeta
        self.database_dir = database_dir

    @abc.abstractmethod
    def get_database_names(self):
        pass

    @abc.abstractmethod
    def get_path(self, database_name):
        pass

    @abc.abstractmethod
    def get_phenotype(self, gene):
        pass