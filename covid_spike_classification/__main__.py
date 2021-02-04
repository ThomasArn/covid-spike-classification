#!/usr/bin/env python3

import argparse
import datetime
import os
import shutil
import sys
import tempfile

from .config import CSCConfig

from .core import (
    REGIONS,
    DELETIONS,
    basecall,
    map_reads,
    align_fasta,
    check_variants,
    parse_vcf,
    get_deletions
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-D", "--datadir", default=os.path.join(os.getcwd(), "data"),
                        help="Directory containing the ab1 files to call variants on (default: %(default)s).")
    parser.add_argument("--genome", default=False,
                        help="input is a genome file as fasta")
    parser.add_argument("-r", "--reference", default=os.path.join(os.getcwd(), "ref", "NC_045512.fasta"),
                        help="Reference FASTA file to use (default: %(default)s).")
    parser.add_argument("-o", "--outdir",
                        default=datetime.datetime.now().strftime("%Y-%m-%d"),
                        help="File to write result CSV and fastq files to (default: %(default)s).")
    parser.add_argument("-q", "--quiet", action="store_true", default=False,
                        help="Suppress noisy output from the tools run")
    parser.add_argument("-d", "--debug", action="store_true", default=False,
                        help="Debug mode: Keep bam file around when the parsing crashes")
    parser.add_argument("--show-unexpected", action="store_true", default=False,
                        help="Show unexpected mutations instead of reporting 'no known mutation'")
    parser.add_argument("-z", "--zip-results", action="store_true", default=False,
                        help="Create a zipfile from the output directory instead of the output directory.")
    args = parser.parse_args()

    config = CSCConfig.from_args(args)
    
    os.makedirs(args.outdir, exist_ok=True)
    if args.genome is False:
        with tempfile.TemporaryDirectory() as tmpdir:
            basecall(tmpdir, config)
            map_reads(tmpdir, config)
            check_variants(tmpdir, config)
            if config.zip_results:
                shutil.make_archive(config.outdir, "zip", root_dir=config.outdir)
                shutil.rmtree(config.outdir, ignore_errors=True)

    else:
        with tempfile.TemporaryDirectory() as tmpdir:
            align_fasta(tmpdir, config)
            parse_vcf(tmpdir, config)
            if config.zip_results:
                shutil.make_archive(config.outdir, "zip", root_dir=config.outdir)
                shutil.rmtree(config.outdir, ignore_errors=True)


if __name__ == "__main__":
    main()
