#!/usr/bin/env python


"""Provide a command line tool to validate and transform tabular samplesheets."""


import argparse
import csv
import logging
import sys
import re
from collections import Counter
from pathlib import Path


logger = logging.getLogger()


class RowChecker:
    """
    Define a service that can validate and transform each given row.

    Attributes:
        modified (list): A list of dicts, where each dict corresponds to a previously
            validated and transformed row. The order of rows is maintained.

    """

    VALID_FORMATS = (
        ".cram",
        ".bai",
        ".crai",
        ".bed",
        ".ped"
    )

    def __init__(
        self,
        sample_col="sample",
        first_col="cram",
        second_col="crai",
        third_col="bed",
        fourth_col="ped",
        **kwargs,
    ):
        """
        Initialize the row checker with the expected column names.

        Args:
            sample_col (str): The name of the column that contains the sample name
                (default "sample").
            first_col (str): The name of the column that contains the full path to the CRAM file (default "cram")
            second_col (str): The name of the column that contains the full path to the CRAI file (default "crai").
            third_col (str): The name of the column that contains the full path to the BED file (default "bed").
            fourth_col (str): The name of the column that contains the full path to the PED file (default "ped").

        """
        super().__init__(**kwargs)
        self._sample_col = sample_col
        self._first_col = first_col
        self._second_col = second_col
        self._third_col = third_col
        self._fourth_col = fourth_col
        self._seen = set()
        self.modified = []

    def validate_and_transform(self, row):
        """
        Perform all validations on the given row and insert the read pairing status.

        Args:
            row (dict): A mapping from column headers (keys) to elements of that row
                (values).

        """
        self._validate_sample(row)
        self._validate_first(row)
        self._validate_third(row)
        self._validate_fourth(row)
        self._seen.add((row[self._sample_col], row[self._first_col]))
        self.modified.append(row)

    def _validate_sample(self, row):
        """Assert that the sample name exists and convert spaces to underscores."""
        assert len(row[self._sample_col]) > 0, "Sample input is required."
        # Sanitize samples slightly.
        row[self._sample_col] = row[self._sample_col].replace(" ", "_")

    def _validate_first(self, row):
        """Assert that the CRAM entry is non-empty and has the right format."""
        assert len(row[self._first_col]) > 0, "A CRAM file is required."
        self._validate_format(row[self._first_col],[".cram"])

    def _validate_third(self, row):
        """Assert that the BED entry has the right format if it exists."""
        assert len(row[self._third_col]) > 0, "A BED file is required"
        self._validate_format(row[self._third_col],[".bed", ".bed.gz"])

    def _validate_fourth(self, row):
        """Assert that the PED entry has the right format if it exists."""
        assert len(row[self._fourth_col]) > 0, "A PED file is required"
        self._validate_format(row[self._fourth_col],[".ped"])

    def _validate_format(self, filename, extensions):
        """Assert that a given filename has one of the expected extensions."""
        assert any(filename.endswith(extension) for extension in extensions), (
            f'The {str(filename).split(".")[-1].upper()} file has an unrecognized extension: {filename}\n'
            f"It should be one of: {', '.join(extensions)}"
        )

def sniff_format(handle):
    """
    Detect the tabular format.

    Args:
        handle (text file): A handle to a `text file`_ object. The read position is
        expected to be at the beginning (index 0).

    Returns:
        csv.Dialect: The detected tabular format.

    .. _text file:
        https://docs.python.org/3/glossary.html#term-text-file

    """
    peek = handle.read(4096)
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(peek)
    handle.seek(0)
    return dialect


def check_samplesheet(file_in, file_out):
    """
    Validate the general shape of the table, expected columns, and each row.

    Args:
        file_in (pathlib.Path): The given tabular samplesheet. The format can be either
            CSV, TSV, or any other format automatically recognized by ``csv.Sniffer``.
        file_out (pathlib.Path): Where the validated and transformed samplesheet should
            be created; always in CSV format.

    Example:
        This function checks that the samplesheet follows the following structure:

            sample,cram,crai,bed,ped
            SAMPLE_1,SAMPLE_1.cram,SAMPLE_1.crai,SAMPLE_1.bed,FILE.ped
            SAMPLE_2,SAMPLE_2.cram,SAMPLE_2.crai,SAMPLE_2.bed,FILE.ped
            SAMPLE_3,SAMPLE_3.cram,SAMPLE_3.crai,SAMPLE_3.bed,FILE2.ped

    """
    required_columns = {"sample", "cram", "crai", "bed", "ped"}
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_in.open(newline="") as in_handle:
        reader = csv.DictReader(in_handle, dialect=sniff_format(in_handle))
        # Validate the existence of the expected header columns.
        if not required_columns.issubset(reader.fieldnames):
            logger.critical(f"The sample sheet **must** contain the column headers: {', '.join(required_columns)}.")
            sys.exit(1)
        # Validate each row.
        checker = RowChecker()
        for i, row in enumerate(reader):
            try:
                checker.validate_and_transform(row)
            except AssertionError as error:
                logger.critical(f"{str(error)} On line {i + 2}.")
                sys.exit(1)
    header = list(reader.fieldnames)
    # See https://docs.python.org/3.9/library/csv.html#id3 to read up on `newline=""`.
    with file_out.open(mode="w", newline="") as out_handle:
        writer = csv.DictWriter(out_handle, header, delimiter=",")
        writer.writeheader()
        for row in checker.modified:
            writer.writerow(row)


def parse_args(argv=None):
    """Define and immediately parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Validate and transform a tabular samplesheet.",
        epilog="Example: python check_samplesheet.py samplesheet.csv samplesheet.valid.csv",
    )
    parser.add_argument(
        "file_in",
        metavar="FILE_IN",
        type=Path,
        help="Tabular input samplesheet in CSV or TSV format.",
    )
    parser.add_argument(
        "file_out",
        metavar="FILE_OUT",
        type=Path,
        help="Transformed output samplesheet in CSV format.",
    )
    parser.add_argument(
        "-l",
        "--log-level",
        help="The desired log level (default WARNING).",
        choices=("CRITICAL", "ERROR", "WARNING", "INFO", "DEBUG"),
        default="WARNING",
    )
    return parser.parse_args(argv)


def main(argv=None):
    """Coordinate argument parsing and program execution."""
    args = parse_args(argv)
    logging.basicConfig(level=args.log_level, format="[%(levelname)s] %(message)s")
    if not args.file_in.is_file():
        logger.error(f"The given input file {args.file_in} was not found!")
        sys.exit(2)
    args.file_out.parent.mkdir(parents=True, exist_ok=True)
    check_samplesheet(args.file_in, args.file_out)


if __name__ == "__main__":
    sys.exit(main())
