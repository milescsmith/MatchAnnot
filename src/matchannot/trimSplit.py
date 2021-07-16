#!/usr/bin/env python

# Given a fasta file of reads-of-insert, qsub N trimPrimers jobs, each
# of which will process a subset of the file. Combine the outpts when
# all are done.

import datetime
import optparse
import os
import re  # regular expressions

from matchannot import matchannot_logger as logger

VERSION = "20150403.01"

DEF_NJOBS = 16
DEF_TMPDIR = "tmp"

STARTWAIT = 30  # seconds to delay start of qsub'd jobs


def main():

    logger.debug(f"version {VERSION} starting")

    opt, args = getParms()

    makeTempDir(opt.tmpdir)

    nSeqs = countSeqs(opt.input)
    logger.debug(f"{opt.input} contains {int(nSeqs)} sequences")

    seqsPerJob = (nSeqs + opt.njobs - 1) / opt.njobs
    logger.debug(
        f"each of {int(opt.njobs)} jobs will process {int(seqsPerJob)} sequences"
    )

    chunkList = makeFastaChunks(opt, nSeqs, seqsPerJob)

    for chunk in chunkList:
        chunk.makeScript()
        chunk.submitScript()

    submitFinalJobs(opt, chunkList)

    logger.debug("finished")

    return


def makeTempDir(dir):

    if os.path.isdir(dir):
        logger.warning(f"WARNING: temp directory {dir} already exists")
    else:
        os.makedirs(dir)

    return


def countSeqs(filename):
    """Run grep -c ">" on a fasta file to count the sequences it contains."""

    command = f'grep -c ">" {filename}'
    popen_file = os.popen(command)
    response = popen_file.read().strip()
    rc = popen_file.close()
    if rc is not None:
        logger.error(f"command failed, rc={int(rc)}")
        raise RuntimeError
    if not response.isdigit():
        logger.error("grep -c returned:" % response)
        raise RuntimeError

    return int(response)


def makeFastaChunks(opt, nSeqs, seqsPerJob):

    curPos = 1
    thisRec = 0
    chunkList = list()

    fastaIn = open(opt.input, "r")
    line = fastaIn.readline()
    if not line.startswith(">"):
        raise RuntimeError(f"first fasta line is not a header: {line}")

    while line:  # outer loop reads the whole input file

        endPos = min(curPos + seqsPerJob - 1, nSeqs)
        fastaChunk = Chunk(opt, curPos, endPos)
        logger.debug(f"writing chunk {fastaChunk.inputChunkName}")
        with open(fastaChunk.inputChunkName, "w") as chunkOut:

            while line:  # inner loop writes one chunk

                if line.startswith(">"):

                    thisRec += 1
                    if thisRec > seqsPerJob:
                        chunkList.append(fastaChunk)
                        thisRec = 0
                        break  # break inner loop and exit with block, closing output file

                    curPos += 1

                chunkOut.write(line)
                line = fastaIn.readline()

    if thisRec > 0:
        chunkList.append(fastaChunk)

    fastaIn.close()

    return chunkList


def submitFinalJobs(opt, chunkList):

    chunkFiles = [f"{chk.trimmedChunkName} \\\n" for chk in chunkList]

    sh = list()
    sh.append("#!/bin/bash\n\n")
    sh.append("set -o errexit\n")
    sh.append("set -o nounset\n\n")

    sh.append("cat \\\n")
    sh.extend(chunkFiles)
    sh.append(f" > {opt.output}\n")

    if opt.report is not None:
        reportFiles = [f"{chk.reportChunkName} \\\n" for chk in chunkList]
        sh.append("\ncat \\\n")
        sh.extend(reportFiles)
        sh.append(f" > {opt.report}\n")

    finalScriptName = f"{opt.tmpdir}/trim_final.sh"
    handle = open(finalScriptName, "w")
    handle.writelines(sh)
    handle.close()

    deps = ":".join([chk.jobno for chk in chunkList])

    cmd = list()
    cmd.append("qsub")
    cmd.append("-N trim_final")  # job name
    cmd.append("-o trim_final.out")  # output file
    cmd.append("-j oe")  # combine stdout and stderr
    cmd.append("-l nodes=1:ppn=1,walltime=4:00:00")  # resources required
    cmd.append("-d . ")  # working directory (strangely, ./ is not the default)
    cmd.append("-r n")  # do NOT attempt to restart on failure
    cmd.append("-V")  # export all environment variables to job
    cmd.append("-W umask=0002")  # make logs rw-rw-r--
    cmd.append("-m n")  # don't send any mail
    cmd.append(f"-W depend=afterok:{deps}")
    cmd.append(finalScriptName)  # script to run

    command = " ".join(cmd)
    logger.debug(f"running {command}")

    popen_file = os.popen(command)
    response = popen_file.read().strip()
    rc = popen_file.close()
    if rc is not None:
        logger.error(f"command failed, rc={int(rc)}")
        raise RuntimeError

    logger.debug(f"jobno is {response}")

    return response


def getParms():  # use default input sys.argv[1:]

    parser = optparse.OptionParser(usage="%prog [options]", version=VERSION)

    parser.add_option("--input", help="fasta file to be trimmed (required)")
    parser.add_option("--primers", help="fasta file of primers (required)")
    parser.add_option("--output", help="output fasta file (required)")
    parser.add_option("--report", help="output alignments report text file (optional)")
    parser.add_option(
        "--njobs", help="number of jobs to submit (def: %default)", type="int"
    )
    parser.add_option("--tmpdir", help="temporary directory (def: %default)")

    parser.set_defaults(
        njobs=DEF_NJOBS,
        tmpdir=DEF_TMPDIR,
    )

    opt, args = parser.parse_args()

    return opt, args


class Chunk(object):
    """Manage a chunk, including keeping all the file names consistent and in one place."""

    JOBNO_PATTERN = re.compile(r"^(\d+)")  # the number part of job number

    def __init__(self, opt, start, end):

        self.opt = opt
        self.start = start
        self.end = end
        self.jobno = None  # filled in when submitted
        self.jobName = f"trim.{int(start):06}"
        self.scriptName = f"{opt.tmpdir}/trim.{int(start):06}_{int(end):06}.sh"
        self.scriptOutput = f"{opt.tmpdir}/trim.{int(start):06}_{int(end):06}.out"
        self.inputChunkName = "%s/input_chunk.%06d_%06d.fasta" % (
            opt.tmpdir,
            start,
            end,
        )
        self.trimmedChunkName = "%s/trimmed_chunk.%06d_%06d.fasta" % (
            opt.tmpdir,
            start,
            end,
        )
        if opt.report is None:
            self.reportChunkName = None
        else:
            self.reportChunkName = "%s/report_chunk.%06d_%06d.txt" % (
                opt.tmpdir,
                start,
                end,
            )

    def makeScript(self):

        sh = list()
        sh.append("#!/bin/bash\n\n")
        sh.append("set -o errexit\n")
        sh.append("set -o nounset\n\n")
        sh.append("~/work/MatchAnnot/trimPrimers.py \\\n")
        sh.append(f"    --input {self.inputChunkName} \\\n")
        sh.append(f"    --primers {self.opt.primers} \\\n")

        if self.opt.report is not None:
            sh.append(f"    --report {self.reportChunkName} \\\n")

        sh.append(f"    > {self.trimmedChunkName}\n")

        handle = open(self.scriptName, "w")
        handle.writelines(sh)
        handle.close()

        return

    # TODO: we don't need to keep this as it only works for specific clusters
    def submitScript(self):

        # Dependent job submission will fail if parent has already
        # completed. So delay all job startups by a short amount of time.

        startAt = datetime.datetime.now() + datetime.timedelta(0, STARTWAIT)
        startAtStr = startAt.strftime("%Y%m%d%H%M.%S")

        cmd = list()
        cmd.append("qsub")
        cmd.append(f"-N {self.jobName}")  # job name
        cmd.append(f"-o {self.scriptOutput}")  # output file
        cmd.append("-j oe")  # combine stdout and stderr
        cmd.append("-l nodes=1:ppn=1,walltime=4:00:00")  # resources required
        cmd.append(f"-a {startAtStr}")  # delay start, see above
        cmd.append("-d . ")  # working directory (strangely, ./ is not the default)
        cmd.append("-r n")  # do NOT attempt to restart on failure
        cmd.append("-V")  # export all environment variables to job
        cmd.append("-W umask=0002")  # make logs rw-rw-r--
        cmd.append("-m n")  # don't send any mail
        cmd.append(self.scriptName)  # script to run

        command = " ".join(cmd)
        logger.debug(f"running {command}")

        popen_file = os.popen(command)
        response = popen_file.read().strip()
        rc = popen_file.close()
        if rc is not None:
            logger.error(f"command failed, rc={int(rc)}")
            raise RuntimeError

        match = re.match(Chunk.JOBNO_PATTERN, response)
        if match is None:
            logger.error(f"invalid job sequence number: {jobSeqStr}")
            raise RuntimeError

        response = match.group(1)
        logger.debug(f"jobno is {response}")
        self.jobno = response
        return response


if __name__ == "__main__":
    main()
