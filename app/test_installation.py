"""
A script that imports all of the modules we will need.
If this script says you don't have something,
you will need to install it using the given command.
"""

try:
    import numpy
except:
    print("Numpy was not found.")
    print("Install it by doing: python -m pip install numpy ")


try:
    import pandas
except:
    print("pandas was not found.")
    print("Install it by doing: python -m pip install pandas ")


try:
    import argparse
except:
    print("argparse was not found.")
    print("Install it by doing: python -m pip install argparse ")