'''
My awesome script.
'''
import argparse
import logging

import pandas as pd

from schrodinger.application.desmond.constants import WT_IDENTIFIER
from schrodinger.utils import log

# Local, depends on SCHRODINGER_PYTHONPATH
from pptc.util import mutation


# Configure logging
logger = log.get_output_logger(__name__)
logger.level = logging.WARNING


def rebase_variant(variant_str: str,
                   base_variant_str: str) -> str:
    '''
    Rebase mutations in a DataFrame.

    This functions similarly to the `git rebase` command.

    It updates the specified 1-letter `variant` mutation string such that its
    mutation(s) are made in the context of the `rebase_variant`.  This
    essentially applies the given `variant` mutations on top of the
    `rebase_variant` reference point, then reports the final mutation string
    relative to the true wild-type variant.

    For example, if `base_variant_str` is given as `A-A123T`, this means the string
    "WT" actually represents a protein with a mutation of A:A123T relative to
    the "true" wild-type protien.  Similarly, the string `A-S128Q` represents
    a double mutant, and the string `A-T123A` represents a WT-reversion
    mutation.  In this case, the function would convert the following input
    variants to produce the corresponding results:

        variant         result
        -------         ------
        WT              A-A123T
        A-S128Q         A-A123T,A-S128Q
        A-T123A         WT

    :param variant: variant string (1-letter code) to be updated
    :param base_variant: variant string (1-letter code) on top of which the
        variant mutations should be rebased
    :return: rebased variant string (1-letter code)
    '''
    variant = mutation.variant_from_str(variant_str)
    base_variant = mutation.variant_from_str(base_variant_str)
    rebased_variant = base_variant + variant
    return rebased_variant.to_str(chain_sep="-")


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('in_csv', type=str,
                        help='input CSV file')
    parser.add_argument('-v', type=str, dest='variant_colname',
                        help='variant column name (1-letter residue code)')
    parser.add_argument('-m', type=str, dest='mutation_str',
                        help='1-letter mutation string for mutations present '
                        'in the WT node of the input CSV')
    parser.add_argument('-debug', action='store_true',
                        help='print debugging info')
    return parser.parse_args(argv)


def main(argv=None):
    '''Main workflow, optionally callable like subprocess with list of args.'''
    args = parse_args(argv)

    # Adjust logging
    if args.debug:
        logger.setLevel(logging.DEBUG)

    # Read the DataFrame and extract the list of variants
    logger.info(f'Reading {args.in_csv} to DataFrame')
    df = pd.read_csv(args.in_csv)
    variants = df[args.variant_colname].tolist()

    # Rebase the variants and update the DataFrame
    logger.info(f'Rebasing variants using {args.mutation_str}')
    rebased_variants = [
        rebase_variant(variant, args.mutation_str)
        for variant in variants
    ]
    df[args.variant_colname] = rebased_variants

    # Write the updated DataFrame to a new CSV file
    out_csv = args.in_csv.replace('.csv', '_rebased.csv')
    df.to_csv(out_csv, index=False)
    logger.warning(f'- Wrote {out_csv}')


if __name__ == '__main__':
    main()
