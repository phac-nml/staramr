"""
An Exception to be raised when an error is encountered while running BLAST/makeblastdb.
"""


class BlastProcessError(Exception):

    def __init__(self, msg, called_process_error_exception):
        """
        Constructs a new BlastProcessError
        :param msg: The Exception message.
        :param called_process_error_exception: The source exception of type subprocess.CalledProcessError
        """
        super().__init__(
            "{}\ncommand={}\nstderr={}".format(msg, called_process_error_exception.cmd,
                                               called_process_error_exception.stderr))
