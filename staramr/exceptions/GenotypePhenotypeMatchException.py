"""
An Exception to be raised when there are issues matching a genotype to a phenotype.
"""


class GenotypePhenotypeMatchException(Exception):

    def __init__(self, msg):
        """
        Constructs a new GenotypePhenotypeMatchException
        :param msg: The Exception message.
        """
        super().__init__(msg)
