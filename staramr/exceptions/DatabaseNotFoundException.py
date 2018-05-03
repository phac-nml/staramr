"""
An Exception to be raised when a database could not be found.
"""


class DatabaseNotFoundException(Exception):

    def __init__(self, msg):
        """
        Constructs a new DatabaseNotFoundException
        :param msg: The Exception message.
        """
        super().__init__(msg)
