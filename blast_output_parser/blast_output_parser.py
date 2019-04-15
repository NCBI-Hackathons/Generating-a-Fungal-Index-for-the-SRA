#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys


def main(blastoutputfilename, outputfilename):
    """Parse BLAST output
    # > DRR001369.139710 GEKBKKN08JJ5PQ length=777
    # Length=777
    #
    #  Score = 1376 bits (745),  Expect = 0.0
    #  Identities = 745/745 (100%), Gaps = 0/745 (0%)
    #  Strand=Plus/Plus
    """

    assert os.path.exists(blastoutputfilename)
    assert not os.path.exists(outputfilename)

    id1 = 'id1'
    id2 = 'id2'
    length = 'length'
    score = 'score'
    length2 = 'length2'
    e_value = 'e_value'
    identities = 'identities'
    gaps = 'gaps'
    strand = 'strand'

    with open(outputfilename, 'w') as output_fp:

        header_written = False

        def dump():
            nonlocal header_written

            if id1 != 'id1' or not header_written:
                to_dump = "%s %s %s %s %s %s %s %s %s\n" % (
                    id1, id2, length, score, length2, e_value, identities,
                    gaps, strand)
                output_fp.write(to_dump)
                header_written = header_written or id1 == 'id1'

        with open(blastoutputfilename) as blast_fp:
            for line in blast_fp.readlines():

                if line.startswith(">"):
                    dump()
                    _, id1, id2, length = line.replace("length=", "").split()

                elif line.startswith(" Score = "):
                    items = line.split()
                    score = items[2]
                    length2 = items[4][1:-2]
                    e_value = items[7]
                elif line.startswith(" Identities = "):
                    items = line.split()
                    identities = items[2]
                    gaps = items[6]
                elif line.startswith(" Strand="):
                    strand = line[8:-1]
            dump()


if __name__ == "__main__":
    assert len(sys.argv) == 3
    main(sys.argv[1], sys.argv[2])
