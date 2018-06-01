"""
An Exception to be raised when invalid blast positions are encountered.
"""


class InvalidPositionException(Exception):

    def __init__(self, msg):
        """
        Constructs a new InvalidPositionException
        :param msg: The Exception message.
        """
        super().__init__(msg)
