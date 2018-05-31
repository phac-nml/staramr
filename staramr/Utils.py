from typing import Dict

"""
Utility functions that don't fit anywhere else.
"""


def get_string_with_spacing(data: Dict[str, str]) -> str:
    """
    Gets a string representation of a list of key/value pairs (as OrderedDictionary) with proper spacing between key/values.
    :param data: A Dictionary containing key/value pairs.
    :return: A string representation of the Dictionary.
    """
    max_width = max([len(k) for k in data])
    return '\n'.join(('{} = {}'.format(k.ljust(max_width), v) for k, v in data.items())) + '\n'
