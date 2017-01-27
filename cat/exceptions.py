class UserException(Exception):
    """generic exception to use when a user makes a mistake"""
    pass


class ToolMissingException(UserException):
    """exception to use when a tool is missing, usually checked in a task validate() method"""
    pass


class InputMissingException(UserException):
    """exception to use when input data are missing"""
    pass


class InvalidInputException(UserException):
    """exception to use when something about the input is invalid"""
    pass


class MissingFileException(UserException):
    """exception to use when a input file is missing"""

