import abc

"""
Abstract class for any sub-commands for the command-line application.
"""


class SubCommand:

    def __init__(self, arg_parser, script_dir):
        """
        Creates a new SubCommand instance.
        :param arg_parser: The argparse.ArgumentParser to use.
        :param script_dir: The directory containing the main application script.
        """
        __metaclass__ = abc.ABCMeta
        self._root_arg_parser = arg_parser
        self._script_dir = script_dir

        self._setup_args(arg_parser)
        arg_parser.set_defaults(run_command=self.run)

    @abc.abstractmethod
    def _setup_args(self, arg_parser):
        pass

    @abc.abstractmethod
    def run(self, args):
        """
        Runs this sub-command with the passed arguments.
        :param args: The dictionary, as returned from argparse.ArgumentParser.parse_args()
        :return: None
        """
        pass
