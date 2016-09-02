def format_ratio(numerator, denominator, num_digits=None):
    """
    Convenience function that converts two numbers to a ratio.
    Handles dividing by zero, as well as transforming values into floats.
    Rounds the number to the number of num_digits, if requested (not None)
    """
    if denominator == 0:
        return float("nan")
    r = float(numerator) / float(denominator)
    if num_digits is not None:
        r = round(r, num_digits)
    return r
