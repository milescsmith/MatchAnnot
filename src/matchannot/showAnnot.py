#!/usr/bin/env python

# Read annotation file, print selected stuff in human-readable format.

# AUTHOR: Tom Skelly (thomas.skelly@fnlcr.nih.gov)

import optparse
import pickle as pickle

from matchannot import matchannot_logger as logger

from . import Annotations as anno

VERSION = "20150417.01"


def main():

    logger.debug(f"version {VERSION} starting")

    opt, args = getParms()

    if opt.gtfpickle is not None:
        handle = open(opt.gtfpickle, "r")
        pk = pickle.Unpickler(handle)
        annotList = pk.load()
        handle.close()
    else:
        annotList = anno.AnnotationList(opt.gtf)

    geneList = annotList.getGene(opt.gene)
    if geneList is None:
        print(f"gene {opt.gene} not found in annotations")
    elif len(geneList) != 1:
        print(
            "there are %d occurrences of gene %s in annotations"
            % (len(geneList), opt.gene)
        )
    else:
        geneEnt = geneList[0]

        print("gene:    ", end=" ")
        printEnt(geneEnt)

        for transEnt in geneEnt.getChildren():
            print("\ntr:      ", end=" ")
            printTran(transEnt)

            for exonEnt in transEnt.getChildren():
                print("exon:    ", end=" ")
                printEnt(exonEnt)

    logger.debug("finished")

    return


def printEnt(ent):

    print(
        "%-15s  %9d  %9d  %6d" % (ent.name, ent.start, ent.end, ent.end - ent.start + 1)
    )

    return


def printTran(ent):

    print(
        "%-15s  %9d  %9d  %6d"
        % (ent.name, ent.start, ent.end, ent.end - ent.start + 1),
        end=" ",
    )
    if hasattr(ent, "startcodon"):
        print(f" start: {int(ent.startcodon):9}", end=" ")
    if hasattr(ent, "stopcodon"):
        print(f" stop:  {int(ent.stopcodon):9}", end=" ")
    print()

    return


def getParms():  # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage="%prog [options] <fasta_file> ... ")

    parser.add_option("--gtf", help="annotations in gtf format")
    parser.add_option("--gtfpickle", help="annotations in pickled gtf format")
    parser.add_option("--gene", help="gene to print")

    parser.set_defaults(
        gtf=None,
        gtfpickle=None,
        gene=None,
    )

    opt, args = parser.parse_args()

    return opt, args


if __name__ == "__main__":
    main()
