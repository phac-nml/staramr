"""
An Exception to be raised during command-line parsing.
"""


class CommandParseException(Exception):

    def __init__(self, msg, parser, print_help=False):
        """
        Constructs a new CommandParseException
        :param msg: The Exception message.
        :param parser: The argparse.ArgumentParser for the particular subparser.
        :param print_help: Whether or not to print a help statement when catching this exception.
        """
        super().__init__(msg)
        self._parser = parser
        self._print_help = print_help

    def get_parser(self):
        """
        Gets the argparse.ArgumentParser for the particular subparser that threw the Exception.
        :return: The argparse.ArgumentParser
        """
        return self._parser

    def print_help(self):
        """
        Whether or not to print a help statement when handling this exception.
        :return: True if help should be printed, False otherwise.
        """
        return self._print_help
