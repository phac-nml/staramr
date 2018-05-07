"""
An Exception to be raised when an error is encountered with the database.
"""


class DatabaseErrorException(Exception):

    def __init__(self, msg):
        """
        Constructs a new DatabaseErrorException
        :param msg: The Exception message.
        """
        super().__init__(msg)
