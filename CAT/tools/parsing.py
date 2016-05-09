"""
Library of parsing things related to argparse.
"""

import argparse
import os
import errno
from collections import defaultdict


class HashableNamespace(argparse.Namespace):
    """
    Adds a __hash__ function to argparse's Namespace. Follows best practices for implementation of __hash__.
    """
    def __hash__(self):
        xor_fn = lambda x, y: x ^ hash(y)
        val_iter = self.__dict__.itervalues()
        first = hash(val_iter.next())
        return reduce(xor_fn, val_iter, first) ^ hash(tuple(self.__dict__.values()))


class NamespaceDictAction(argparse.Action):
    """
    http://stackoverflow.com/questions/34930630/grouping-an-unknown-number-of-arguments-with-argparse#34930706
    Asking questions on SO is a good idea!

    This action runs in two modes - dict and defaultdict. In defaultdict mode, a key can be passed more than once.
    In dict mode, the last key will be used. Overall, produces a single Namespace as the value for a argument.
    """
    def __init__(self, *args, **kwargs):
        self.mode = kwargs.pop('mode', 'dict')  # Will use 'dict' as default
        super(NamespaceDictAction, self).__init__(*args, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        try:
            if self.mode == 'dict':
                arg_vals = HashableNamespace(**dict(v.split('=') for v in values))
            elif self.mode == 'defaultdict':
                d = defaultdict(list)
                for v in values:
                    v = v.split('=')
                    d[v[0]].append(v[1])
                df = {x: tuple(y) for x, y in d.iteritems()}
                arg_vals = HashableNamespace(**df)
            else:
                raise NotImplementedError("only dict or defaultdict")
        except TypeError:
            raise RuntimeError('Group {} appears to be incorrectly formatted'.format(values))
        setattr(namespace, self.dest, arg_vals)


class NamespaceAction(argparse.Action):
    """
    http://stackoverflow.com/questions/34930630/grouping-an-unknown-number-of-arguments-with-argparse#34930706
    Asking questions on SO is a good idea!
    This modified action allows me to group together the four key-value pairs passed to --geneSets as a nested
    namespace.
    """
    def __init__(self, *args, **kwargs):
        self.mode = kwargs.pop('mode', 'dict')  # Will use 'dict' as default
        super(NamespaceAction, self).__init__(*args, **kwargs)

    def __call__(self, parser, namespace, values, option_string=None):
        # The default value is often set to None rather than an empty list.
        current_arg_vals = getattr(namespace, self.dest, []) or []
        setattr(namespace, self.dest, current_arg_vals)
        arg_vals = getattr(namespace, self.dest)
        try:
            arg_vals.append(HashableNamespace(**dict(v.split('=') for v in values)))
        except TypeError:
            raise RuntimeError('Group {} appears to be incorrectly formatted'.format(values))


class FileArgumentParser(argparse.ArgumentParser):
    """
    http://codereview.stackexchange.com/questions/28608/checking-if-cli-arguments-are-valid-files-directories-in-python

    This modified parser allows me to validate input files as existing without using the builtin FileType. FileType is
    bad because it wants to actually open a file handle, which is not useful to us here.

    """
    def __is_valid_file(self, parser, arg):
        if not os.path.isfile(arg):
            parser.error('The file {} does not exist!'.format(arg))
        else:
            # File exists so return the filename
            return arg

    def __is_valid_directory(self, parser, arg):
        if not os.path.isdir(arg):
            parser.error('The directory {} does not exist!'.format(arg))
        else:
            # File exists so return the directory
            return arg

    def __mkdir_p(self, parser, arg):
        try:
            os.makedirs(arg)
        except OSError as exc:  # Python >2.5
            if exc.errno == errno.EEXIST and os.path.isdir(arg):
                pass
            else:
                raise RuntimeError('Error: was unable to make {}'.format(arg))
        return arg

    def add_argument_with_check(self, *args, **kwargs):
        # Look for your FILE or DIR settings
        if 'metavar' in kwargs and 'type' not in kwargs:
            if kwargs['metavar'] is 'FILE':
                type=lambda x: self.__is_valid_file(self, x)
                kwargs['type'] = type
            elif kwargs['metavar'] is 'DIR':
                type=lambda x: self.__is_valid_directory(self, x)
                kwargs['type'] = type
        self.add_argument(*args, **kwargs)

    def add_argument_with_mkdir_p(self, *args, **kwargs):
        if 'metavar' in kwargs and 'type' not in kwargs:
            if kwargs['metavar'] is 'FILE':
                type=lambda x: self.__mkdir_p(self, os.path.dirname(x))
                kwargs['type'] = type
            elif kwargs['metavar'] is 'DIR':
                type=lambda x: self.__mkdir_p(self, x)
                kwargs['type'] = type
        self.add_argument(*args, **kwargs)
