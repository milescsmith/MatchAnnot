#!/usr/bin/env python

# Read a .sam file, output several fasta files containing reads which
# do not overlap in the SAM alignment.

# The SAM file is assumed to be sorted by chr and position.

# AUTHOR: Tom Skelly (thomas.skelly@fnlcr.nih.gov)

import optparse
import os
import string
import sys

from matchannot import matchannot_logger as logger

from . import CigarString as cs

VERSION = "20150501.01"

FLAG_NOT_ALIGNED = 0x04  # SAM file flags
FLAG_REVERSE = 0x10
FLAG_SECONDARY = 0x100

FASTA_WRAP = 60
BOWTIE_BUILD = (
    "/is2/projects/pacbio/static/software/packages/bowtie2-2.2.3/bowtie2-build"
)


def main():

    logger.debug(f"version {VERSION} starting")

    opt, args = getParms()

    if len(args) > 0:
        logger.debug(f"reading SAM file {args[0]}")
        handle = open(args[0], "r")
    else:
        logger.debug("reading SAM data from stdin")
        handle = sys.stdin

    totReads = 0  # assorted counters
    totAlign = 0
    totUnsup = 0
    lastPos = dict()  # current position in SAM file by chr (for sort check)
    fileSeq = 0
    fastaList = list()  # list of currently open output fasta files

    compltab = string.maketrans("ACGT", "TGCA")  # for reverse-complementing reads

    fastaNonAligned = tiledFasta(f"{opt.prefix}.nonaligned.fasta")

    for line in handle:  # main loop: read the SAM file

        if line.startswith("@"):
            continue

        totReads += 1

        lineFields = line.split("\t")  # just split it once, not 6 times in list comp
        readName, flags, chr, start, cigarString, bases = [
            lineFields[i] for i in (0, 1, 2, 3, 5, 9)
        ]
        flags = int(flags)
        start = int(start)

        if readName.find("|f1p0|") != -1:
            totUnsup += 1
            continue  # skip unsupported reads

        if flags & FLAG_NOT_ALIGNED:
            fastaNonAligned.addRead(readName, "*", 0, 0, bases)
            continue

        totAlign += 1

        strand = "-" if (flags & FLAG_REVERSE) else "+"
        cigar = cs.CigarString(cigarString, start)
        end = start + cigar.genomicLength() - 1
        # -1 to report last base, rather than last+1

        if start < lastPos.get(chr, 0):
            raise RuntimeError("SAM file is not sorted by position")
        lastPos[chr] = start

        if strand == "-":
            bases = bases[::-1].translate(compltab)

        found = False
        for fasta in fastaList:
            if fasta.fits(chr, start):
                fasta.addRead(readName, chr, start, end, bases)
                found = True
                break

        if not found:
            fileSeq += 1
            fasta = tiledFasta(f"{opt.prefix}.{int(fileSeq):03}.fasta")
            fastaList.append(fasta)
            fasta.addRead(readName, chr, start, end, bases)

    handle.close()

    fastaNonAligned.close()
    for fasta in fastaList:
        fasta.close()
        fasta.makeRef()

    logger.debug(
        "found %d reads, of which %d aligned, %d were f1p0",
        totReads,
        totAlign,
        totUnsup,
    )
    logger.debug("finished")


def getParms():  # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage="%prog [options] <SAM_file> ... ")

    parser.add_option("--prefix", help="start of fasta file name (required)")

    parser.set_defaults(
        prefix=None,
    )

    opt, args = parser.parse_args()

    return opt, args


class tiledFasta(object):
    def __init__(self, name):

        self.name = name
        self.handle = open(name, "w")
        self.lastPos = dict()

        logger.debug(f"opened {name}")

    def fits(self, chr, start):

        return start > self.lastPos.get(chr, 0)

    def addRead(self, readName, chr, start, end, bases):

        self.lastPos[chr] = end
        self.handle.write(f">{readName}  {chr}  {int(start)}  {int(end)}\n")
        for ix in range(0, len(bases), FASTA_WRAP):
            self.handle.write(f"{bases[ix:ix + FASTA_WRAP]}\n")

    def close(self):

        self.handle.close()
        logger.debug(f"closed {self.name}")

    def makeRef(self):
        """Invoke bowtie-build on the fasta file."""

        if not self.handle.closed:
            self.close()

        command = f"{BOWTIE_BUILD} {self.name} {self.name} > {self.name}.out 2>&1"
        logger.debug(command)

        buildOut = os.popen(
            command
        )  # this should return nothing, since we've redirected the output
        rc = buildOut.close()
        if rc is not None:
            raise RuntimeError(f"bowtie2-build failed: {int(rc)}")


if __name__ == "__main__":
    main()
