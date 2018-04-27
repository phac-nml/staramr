import pandas


class ARGDrugTable:

    def __init__(self, file):
        self.file = file
        self.data = pandas.read_csv(file, sep="\t")
