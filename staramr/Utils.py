"""
Utility functions that don't fit anywhere else.
"""


def get_string_with_spacing(data):
    """
    Gets a string representation of a list of key/value pairs (stored as a sub-list) with proper spacing between key/values.
    :param data: A list containing sub-lists of key/value pairs.
    :return: A string representation of the list.
    """
    s = ''
    max_width = max([len(w[0]) for w in data])

    for x, y in data:
        d = dict(x=x.ljust(max_width), y=y)
        s += '{x} = {y}\n'.format(**d)

    return s
