#!/usr/bin/env python

"""
Simple script to test out using pyugrid

Example from Rich Signell.

"""

from __future__ import (absolute_import, division, print_function)

import pyugrid

# FVCOM dataset.
url = 'http://testbedapps-dev.sura.org/thredds/dodsC/in/usf/fvcom/ike/ultralite/vardrag/wave/2d'  # noqa

# Get the datasets:
# Note: this reads the whole thing in to memory at once!
# Maybe we don't want to do that.
print("Loading data: This could take a while...")
ug = pyugrid.UGrid.from_ncfile(url)

# What's in there?
print("There are %i nodes" % ug.nodes.shape[0])
print("There are %i edges" % ug.edges.shape[0])
print("There are %i faces" % ug.faces.shape[0])

print('The start of the "connectivity array":', ug.faces[:10])
