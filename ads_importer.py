#!/usr/bin/env python
# coding: utf-8

import sys, os, re
import requests

# this script will query ADS to collect the citations used in a LaTeX file
# into a BibTeX file

# you can add your own command-line option parsing and help line
# using getopt or optparse; we'll assume a command-line single argument FILE
# (given by sys.argv[1]), read from FILE.aux, and write to FILE.bib

auxfile = sys.argv[1] + '.aux'
bibfile = sys.argv[1] + '.bib'

# FIRST, we'll collect all citation keys from FILE.aux;
# citations will look like \citation{2004PhRvD..69j4017P,2004PhRvD..69j4017P}

cites = set()   # start with an empty set (like list, not ordered, no repetitions)

for line in open(auxfile,'r'):              # Python idiom—loop over every line from file
    m = re.search(r'\\citation\{(.*)\}',line)   # match \citation{...}, collect the ...
                                                # note that we escape \, {, and }
    if m:
        cites.update(m.group(1).split(','))     # if there's a match, split string by commas
                                                # add the resulting keys to set

# check: print "Seek:", cites

# SECOND, we'll check what refs we have already in FILE.bib, to avoid
# repetitive querying of ADS; references will look like
# @TYPE{key,
# ...

haves = []

if os.path.isfile(bibfile):                 # the bibfile exists...
    for line in open(bibfile,'r'):
        m = re.search(r'@.*?\{(.*),',line)  # .*\{ means "any # of any char followed by {";
                                            # .*?\{ means "the shortest string matching
                                            #              any # of any char followed by {"
        if m:
            haves.append(m.group(1))        # remember: append item to list
                                            #           extend list with list
                                            #           add item to set
                                            #           update set with list

# check: print "Have:", haves

# THIRD, we'll query ADS for all the keys that are not in haves,
# and write the resulting bibtex entries to FILE.bib

bibtex = open(bibfile,'a')      # open for appending (very C-like)

for c in cites:
    if c not in haves:
        r = requests.post('http://ukads.nottingham.ac.uk/cgi-bin/nph-bib_query',    # CGI script
                          data = {'bibcode': c, 'data_type': 'BIBTEX'} )        # CGI parameters (note pretty indent)

# we could also have done a (more restrictive) GET HTTP request
# r = requests.get('http://adsabs.harvard.edu/cgi-bin/nph-bib_query?bibcode=%s&data_type=BIBTEX' % c)
# where % yields the Python %-formatting for strings

        if r.ok:                                            # found it!
            m = re.search(r'(@.*\})',r.content.decode('utf-8'),re.DOTALL)   # get everything after the first @
                                                            # until the last }
            bibtex.write(m.group(1) + '\n')                 # write to file with extra newline
            # check: print "Found:", c
        else:
            bibtex.write('@MISC{%s,note="{%s not found in ADS!}"}' % (c,c))
                                                            # record not found,
                                                            # we'll write a useful note in bibtex
            # check: print "Not found:", c

bibtex.close()                              # close file and flush write buffer

# FOURTH, we're done!
