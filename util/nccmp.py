#!/usr/bin/env python
from sys import argv, exit
from getopt import getopt, GetoptError

tol = 0.0   # default tolerance is perfection

def success():
    print "Files are the same within tolerance %.1e" % tol
    exit(0)

def failure():
    print "Files are different."
    exit(1)

usage="""nccmp.py compares NetCDF files by absolute max norms of difference of variables
usage:
  nccmp.py foo.nc bar.nc            compare all variables
  nccmp.py -v A,B foo.nc bar.nc     compare variables A and B
  nccmp.py -x C foo.nc bar.nc       compare all variables except C
  nccmp.py -t 1e-6 foo.nc bar.nc    use tolerance 1e-6 instead of default of 0"""

def usagefailure(message):
    print message
    print
    print usage
    exit(2)

def compare_vars(nc1, nc2, name, tol):
    from numpy import squeeze

    try:
        var1 = squeeze(nc1.variables[name][:])
    except:
        usagefailure("ERROR: VARIABLE '%s' NOT FOUND IN FILE 1" % name)
    try:
        var2 = squeeze(nc2.variables[name][:])
    except:
        usagefailure("ERROR: VARIABLE '%s' NOT FOUND IN FILE 2" % name)

    try:
        delta = abs(var1 - var2).max()
    except:
        usagefailure("ERROR: VARIABLE '%s' OF INCOMPATIBLE SHAPES (?) IN FILES" % name)

    # The actual check:
    if (delta > tol):
        print "name = %s, delta = %e, tol = %e" % (name, delta, tol)
        failure()


def compare(file1, file2, variables, exclude, tol):
    try:
        from netCDF4 import Dataset as NC
    except:
        from netCDF3 import Dataset as NC

    from numpy import unique, r_

    try:
        nc1 = NC(file1, 'r')
    except:
        usagefailure("ERROR: FILE '%s' CANNOT BE OPENED FOR READING" % file1)
    try:
        nc2 = NC(file2, 'r')
    except:
        usagefailure("ERROR: FILE '%s' CANNOT BE OPENED FOR READING" % file2)

    if (exclude == False):
        if len(variables) == 0:
            vars1 = nc1.variables.keys()
            vars2 = nc2.variables.keys()
            variables = unique(r_[vars1, vars2])

        for each in variables:
            compare_vars(nc1, nc2, each, tol)
    else:
        vars1 = nc1.variables.keys()
        vars2 = nc2.variables.keys()
        vars = unique(r_[vars1, vars2])

        for each in vars:
            if (each in variables):
                continue
            compare_vars(nc1, nc2, each, tol)

if __name__ == "__main__":
    from numpy import double
    try:
      opts, args = getopt(argv[1:], "t:v:x", ["help","usage"])
    except GetoptError:
      usagefailure('ERROR: INCORRECT COMMAND LINE ARGUMENTS FOR nccmp.py')
    file1 = ""
    file2 = ""
    variables = []
    exclude = False
    for (opt, arg) in opts:
        if opt == "-t":
            tol = double(arg)
        if opt == "-x":
            exclude = True
        if opt == "-v":
            variables = arg.split(",")
        if opt in ("--help", "--usage"):
            print usage
            exit(0)

    if len(args) != 2:
        usagefailure('ERROR: WRONG NUMBER OF ARGUMENTS FOR nccmp.py')

    compare(args[0],args[1], variables, exclude, tol)

    success()

