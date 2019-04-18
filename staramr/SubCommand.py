import abc
import logging
import coloredlogs

"""
Abstract class for any sub-commands for the command-line application.
"""


class SubCommand:

    def __init__(self, subparser, script_name):
        """
        Creates a new SubCommand instance.
        :param subparser: The subparser to use.  Generated from argparse.ArgumentParser.add_subparsers().
        :param script_name: The name of the script being run.
        """
        __metaclass__ = abc.ABCMeta
        self._subparser = subparser
        self._script_name = script_name

        arg_parser = self._setup_args(subparser)
        arg_parser.set_defaults(run_command=self.run)

        self._root_arg_parser = arg_parser

    @abc.abstractmethod
    def _setup_args(self, arg_parser):
        pass

    def run(self, args):
        """
        Runs this sub-command with the passed arguments.
        :param args: The dictionary, as returned from argparse.ArgumentParser.parse_args()
        :return: None
        """
        if args.verbose:
            coloredlogs.install(level='DEBUG', fmt='%(asctime)s %(levelname)s %(name)s.%(funcName)s,%(lineno)s: %(message)s')
        else:
            coloredlogs.install(level='INFO', fmt='%(asctime)s %(levelname)s: %(message)s')
