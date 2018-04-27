import pandas
from os import path

class ARGDrugTable:

    DEFAULT_DATA_DIR=path.join(path.dirname(__file__), 'data')

    def __init__(self, file):
        self.file = file
        self.data = pandas.read_csv(file, sep="\t")
