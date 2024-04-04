#!/usr/bin/env python3

header = """
..
   DO NOT EDIT: This file was automatically generated by running doc/acknowledge.py

   Edit doc/acknowledge.py, doc/funding.csv, and doc/citing-pism.bib
"""

acknowledgement = """
Acknowledging PISM funding sources
----------------------------------

If you use PISM in a publication then we ask for an acknowledgement of funding and a
citation. However, unless PISM developers are involved in the preparation of the
publication at the usual co-author level, we do not expect co-authorship on PISM-using
papers.

To acknowledge PISM funding please include the statement:
"""

citing = """
Citing
------

To cite PISM please use at least one of Bueler and Brown (2009) or Winkelmann et al.
(2011), below, as appropriate to the application. Also cite PISM as a software as
specified in the file `CITATION.cff <https://github.com/pism/pism/blob/main/CITATION.cff>`_,
or check the 'Cite this repository' link in the sidebar of PISM's `github
repository <https://github.com/pism/pism>`_.

Do not forget to specify the PISM *version* you use. If your results came from source code
modifications to PISM then we request that your publication say so explicitly.

If your study relies heavily on certain PISM sub-models (such as hydrology, calving,
fracture mechanics, thermodynamics) please contact the corresponding author/developer for
information on additional citations.

.. code::
"""

import csv
import time
import sys
import argparse

parser = argparse.ArgumentParser()
parser.description = '''Generate a funding acknowledgment string.'''
parser.add_argument("--manual", action="store_true")
options = parser.parse_args()

year = time.gmtime(time.time())[0]
funding = {}

with open("funding.csv", "r") as f:
    reader = csv.reader(f, skipinitialspace=True, quoting=csv.QUOTE_ALL)

    funding = {}
    for row in reader:
        start_year, end_year, agency, number, _ = row

        try:
            start_year = int(start_year)
            end_year = int(end_year)
        except:
            continue

        # skip grants for which we don't have a number (yet)
        if number.strip() == "":
            continue

        if start_year <= year and year <= end_year:
            try:
                funding[agency].append(number)
            except:
                funding[agency] = [number]


def join(strings):
    assert len(strings) > 0
    if len(strings) == 1:
        return strings[0]
    elif len(strings) == 2:
        return "{} and {}".format(strings[0], strings[1])
    else:
        return join(["{}, {}".format(strings[0], strings[1]),
                     join(strings[2:])])


grants = []
for k, v in funding.items():
    grant = "grant"
    if len(v) > 1:
        grant = "grants"

    grants.append("{agency} {grant} {number}".format(agency=k,
                                                     grant=grant,
                                                     number=join(v)))

if options.manual:
    print(header)
    print("""
Development of PISM is supported by {grants}.""".format(grants=join(grants)))
else:
    print(header)
    print(acknowledgement)
    print("""
    Development of PISM is supported by {grants}.
""".format(grants=join(grants)))
    print(citing)
    with open("citing-pism.bib") as f:
        for line in f:
            # Note: len(line) is one if line == "\n", i.e. it's an empty line.
            sys.stdout.write(line if len(line) == 1 else "   " + line)
