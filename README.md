# pyBlast

pyBlast simplifies working with BLAST from Python by enabling you to directly BLAST two fasta files against each other, without having to format the BLAST database or run BLAST from the command line.

This is a minimal program to BLAST `query.fasta` against `database.fasta` using pyBlast
```
#!/usr/bin/env python3

import sys
from Bio.Blast.Applications import NcbiblastnCommandline
from pyBlast import pyBlast

query = sys.argv[1]
database = sys.argv[2]

cmd = NcbiblastnCommandline(query=query, db=database)

with pyBlast(cmd) as pb:
    for record in pb:
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                print("sequence:", alignment.title)
                print("length:", alignment.length)
                print("e value:", hsp.expect)
                print(hsp.query[0:75] + "...")
                print(hsp.match[0:75] + "...")
                print(hsp.sbjct[0:75] + "...")
```
Behind the scenes, pyBlast will take `database.fasta`, and call `makeblastdb` on the file to create a temporary database the ncbi-blast program can use. It then runs the [NcbiblastnCommandline](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc100) you specified against the database, returning an iterator of [BLAST records](http://biopython.org/DIST/docs/tutorial/Tutorial.html#htoc103) that contain your results. All files will be automatically removed when the `with` block ends, even when the program crashes.

**Important** Because the blast database is created on the fly and removed after each invocation, be careful using pyBlast on huge database files.
