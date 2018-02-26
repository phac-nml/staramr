import abc


class SubCommand:

    def __init__(self, arg_parser, script_dir):
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
        pass
