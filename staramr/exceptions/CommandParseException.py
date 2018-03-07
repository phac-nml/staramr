"""
An Exception to be raised during command-line parsing.
"""


class CommandParseException(Exception):

    def __init__(self, msg, parser):
        """
        Constructs a new CommandParseException
        :param msg: The Exception message.
        :param parser: The argparse.ArgumentParser for the particular subparser.
        """
        super().__init__(msg)
        self._parser = parser

    def get_parser(self):
        """
        Gets the argparse.ArgumentParser for the particular subparser that threw the Exception.
        :return:
        """
        return self._parser
