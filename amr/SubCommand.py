import abc

class SubCommand:

    def __init__(self, arg_parser):
        __metaclass__ = abc.ABCMeta

        self._setup_args(arg_parser)
        arg_parser.set_defaults(run_command=self.run)

    @abc.abstractmethod
    def _setup_args(self, arg_parser):
        pass

    @abc.abstractmethod
    def run(self, args):
        pass