class CommandParseException(Exception):

    def __init__(self, msg, parser):
        super().__init__(msg)
        self._parser = parser

    def get_parser(self):
        return self._parser
