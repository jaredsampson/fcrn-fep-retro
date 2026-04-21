'''
Prepare experimental data
'''
import argparse
import glob
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import gmean, gstd

from schrodinger.application.scisol.fep_gui import unit_convert
from schrodinger.application.scisol.ph_dependence import phdep_utils
from schrodinger.utils import log


# Configure logging
logger = log.get_output_logger(__name__)
logger.level = logging.WARNING

# List of valid affinity units for experimental binding data
VALID_UNITS = [units.name for units in unit_convert.AffinityUnits]

# Input column labels
VARIANT_NAME = 'variant_name'
FC_MUTATIONS = 'fc_mutations'
FCRN_MUTATIONS = 'fcrn_mutations'
SOURCE = 'source'
SOURCE_VARIANT_NAME = 'source_variant_name'
ANTIBODY = 'antibody'
TEMP_C = 'temp_C'
PH = 'pH'
KD_NM = 'Kd_nM'
ERR_KD_NM = 'err_Kd_nM'

# Output column labels
INEQ = 'ineq'
DG = 'exp_dg'
ERR_DG = 'exp_dg_err'
WT_DG = 'wt_exp_dg'
ERR_WT_DG = 'wt_exp_dg_err'
DDG = 'exp_ddg'
ERR_DDG = 'exp_ddg_err'
WT_DG_OFFSET = 'wt_exp_dg_offset'
OFFSET_DG = 'offset_exp_dg'
NMEAS = 'n_meas'
VARIANT = 'variant'

# Output file names
FULL_OUT_CSV = 'full_retrospective_dataset.csv'
VARIANTS_OUT_CSV = 'variants_dataset.csv'

# Other constants
WT = 'WT'
FEP_CHAIN_SEP = '-'
FEP_MUTATION_SEP = ','
CURRENT_WORK = 'Current work'


##################################################
# Utility functions
##################################################


def is_wt(mut_str: str) -> bool:
    '''
    Check if a mutation string is empty or contains only the WT string.
    '''
    return str(mut_str).strip().upper() in [WT, '']


def rms(values):
    '''
    Calculate the root mean square of a list of values.

    :param values: List of values
    :returns: root mean square of the values
    '''
    return np.sqrt(np.mean(np.square(values)))


##################################################
# Variants dataframe functions
##################################################


def read_csv_files(source_path: Path):
    '''
    Read all csv files in the indicated directory and return a single dataframe

    :param source_path: Path to the directory containing the input CSV files
    '''
    glob_str = str(source_path / '*.csv')
    source_csv_files = glob.glob(glob_str)
    source_dfs = []
    for source_csv_file in source_csv_files:
        logger.info(f'Reading {source_csv_file}')
        source_df = pd.read_csv(source_csv_file)
        source_dfs.append(source_df)
        logger.info(f'- added {len(source_df.columns)} columns x '
                    f'{len(source_df.columns)} rows')
    return pd.concat(source_dfs, ignore_index=True)


def cleanup_raw_data(df):
    '''
    Clean up the column names in the dataframe.
    '''
    # Remove NaN strings
    df[VARIANT_NAME] = df[VARIANT_NAME].fillna('')
    df[SOURCE_VARIANT_NAME] = df[SOURCE_VARIANT_NAME].fillna('')


def convert_to_numeric(df):
    '''
    Convert input columns to numeric as needed.
    '''
    # Ensure all temp and pH values are numeric
    df[TEMP_C] = pd.to_numeric(df[TEMP_C], errors='coerce')
    df[PH] = pd.to_numeric(df[PH], errors='coerce')

    # Specialy handling of inequality Kd values; use pH dependence
    # infrastructure to extract inequalities as SidedError integers then
    # convert back to inequality strings.
    df[[KD_NM, INEQ]] = pd.DataFrame(
        df[KD_NM].apply(
            lambda x: phdep_utils.parse_sided_error_value_str(str(x))).tolist(),
        index=df.index)
    df[INEQ] = df[INEQ].apply(lambda x: phdep_utils.SIDED_ERROR_TO_INEQUALITY_MAP[x])

    # Ensure all Kd error values are numeric
    df[ERR_KD_NM] = pd.to_numeric(df[ERR_KD_NM], errors='coerce')

    # Replace Kd error values of 0 with NaN
    df[ERR_KD_NM] = df[ERR_KD_NM].replace({'0':np.nan, 0:np.nan})


def convert_kd_err_to_dg(kd_val, kd_err, dg_val, units):
    '''
    Convert Kd error values to DG error values.
    '''
    logger.debug("Converting Kd error to DG error")
    logger.debug(f'kd_val={kd_val}, kd_err={kd_err}, dg_val={dg_val}, units={units}')
    if any([pd.isna(kd_val), pd.isna(kd_err), pd.isna(dg_val)]):
        return pd.NA

    # Calculate dG range defined by Kd error
    lower_bound = kd_val - kd_err
    upper_bound = kd_val + kd_err

    dg_upper_bound = unit_convert.convert(upper_bound, units)
    upper_diff = abs(dg_val - dg_upper_bound)

    # Lower bound must be positive to be physically meaningful
    if lower_bound <= 0:
        logger.debug(f'Lower bound is nonphysical, using upper bound only')
        logger.debug(f'Kd error = {upper_diff}')
        return upper_diff

    dg_lower_bound = unit_convert.convert(lower_bound, units)
    lower_diff = abs(dg_val - dg_lower_bound)
    max_diff = max(lower_diff, upper_diff)
    logger.debug(f'Kd error = {max_diff}')
    return max_diff


def convert_kd_column_to_dg(df, in_colname, units, out_colname):
    '''
    Convert a column of Kd values to DG values in kcal/mol.

    Units must be one of the values from fep_gui.unit_convert.AffinityUnits.

    :param colname: Name of the column to convert
    :param units: Units of the input values
    '''
    df[out_colname] = df[in_colname].apply(
        lambda x: unit_convert.convert(x, units))


def convert_kd_err_column_to_dg_err(df, in_colname, in_err_colname, units,
                                    out_colname):
    '''
    Convert a column of Kd values to DG values in kcal/mol.

    Units must be one of the values from fep_gui.unit_convert.AffinityUnits.

    '''
    df[out_colname] = df.apply(
        lambda row: convert_kd_err_to_dg(row[in_colname], row[in_err_colname],
                                         row[DG], units),
        axis=1)


def merge_variants_for_current_work(df):
    '''
    Merge all measurements of each variant from "Current work" source into a
    single row.
    '''
    previous_work_df = df[df[SOURCE] != CURRENT_WORK]
    current_work_df = df[df[SOURCE] == CURRENT_WORK]
    meta_cols = [VARIANT_NAME, FC_MUTATIONS, FCRN_MUTATIONS, SOURCE,
                 SOURCE_VARIANT_NAME, ANTIBODY, TEMP_C, PH]
    gb = current_work_df.groupby(meta_cols)
    agg_dict = {
        KD_NM: gmean,
        ERR_KD_NM: gstd,
        INEQ: 'first',
        DG: 'mean',
        ERR_DG: rms,
    }
    merged_variants_df = gb.agg(agg_dict).reset_index()

    out_df = pd.concat([previous_work_df, merged_variants_df],
                       ignore_index=True)
    return out_df


def get_wt_dg_df(df: pd.DataFrame) -> pd.DataFrame:
    '''
    Get the wild type dG and dG error value for each source.
    '''
    # Get the rows where both Fc and FcRn are WT.
    fc_is_wt = df[FC_MUTATIONS].apply(is_wt)
    fcrn_is_wt = df[FCRN_MUTATIONS].apply(is_wt)
    wt_df = df[fc_is_wt & fcrn_is_wt]
    wt_dg_df = pd.DataFrame()
    # Calculate the mean
    wt_dg_df[WT_DG] = wt_df.groupby([SOURCE, ANTIBODY])[DG].mean()
    # Placeholder: use max error value for now
    wt_dg_df[ERR_WT_DG] = (wt_df.dropna(subset=[ERR_DG])
                                .groupby([SOURCE, ANTIBODY])[ERR_DG]
                                .agg(max))
    return wt_dg_df


def add_wt_dg_columns(df):
    '''
    Add wild type dG and dG error values to the dataframe and return the merged
    dataframe.
    '''
    wt_dg_df = get_wt_dg_df(df)
    merged_df = df.join(wt_dg_df, on=[SOURCE, ANTIBODY])
    return merged_df


def add_wt_dg_offset_col(df):
    '''
    Calculate the offset to be added to the individual source wt dG values to
    make them equal to the mean wt dG value for the entire dataset.

    i.e. source_wt_dG + wt_offset = global_mean_wt_dG
    '''
    # Get the mean WT exp dG value for each source.
    source_wt_df = df.groupby([SOURCE, ANTIBODY]).agg({WT_DG: 'mean'})
    # Calculate the mean over all sources.
    global_mean_wt_dg = source_wt_df[WT_DG].mean()
    # Add the offset column in-place.
    df[WT_DG_OFFSET] = global_mean_wt_dg - df[WT_DG]


def add_offset_dg_columns(df):
    '''
    Add offset dG columns to the dataframe.
    '''
    df[OFFSET_DG] = df[DG] + df[WT_DG_OFFSET]


def add_ddg_columns(df):
    '''
    Add ddG and ddG error values to the dataframe and return the merged
    dataframe.
    '''
    df[DDG] = df[DG] - df[WT_DG]
    df[ERR_DDG] = (df[ERR_DG] ** 2 + df[ERR_WT_DG] ** 2) ** 0.5


def cleanup_full_dataset(df: pd.DataFrame) -> pd.DataFrame:
    '''
    Clean up the full dataset, standardizing data formats and removing
    unwanted rows.  Returns a copy of the dataframe.
    '''
    # Drop rows with NaN dG values
    df = df.dropna(subset=[DG])

    # Drop rows with insertion or deletion mutations
    df = df[~df[FC_MUTATIONS].str.contains('ins|del', case=False, na=False)]

    # Keep only rows that are not inequality values
    df = df[df[INEQ] == '']

    # Standardize NaN values
    df.replace(pd.NA, np.nan, inplace=True)

    return df


def prepare_full_expt_df(source_path, units):
    '''
    Prepare experimental data
    '''
    df = read_csv_files(source_path)
    cleanup_raw_data(df)
    convert_to_numeric(df)
    convert_kd_column_to_dg(df, KD_NM, units, DG)
    convert_kd_err_column_to_dg_err(df, KD_NM, ERR_KD_NM, units, ERR_DG)
    # add nmeas column
    df[NMEAS] = pd.NA

    # Calculate ddG values from mean WT dG, averaged per-source
    df = merge_variants_for_current_work(df)
    df = add_wt_dg_columns(df)  # returns a copy after join
    add_ddg_columns(df)

    # Shift each source's dG values so the WT dG values are equal to the
    # global mean WT dG value
    add_wt_dg_offset_col(df)
    add_offset_dg_columns(df)

    # Write the full dataset to a CSV file
    df = cleanup_full_dataset(df)
    return df


##################################################
# Variants dataframe functions
##################################################


def get_fep_variant_str(mut_str: str, chain_id: str,
                        mut_sep: str = '/') -> str:
    '''
    Get the FEP variant string from a mutation string and chain ID.

    If the mutation string is empty or contains only the WT string, return an
    empty string.

    :param mut_str: Mutation string
    :param chain_id: Chain ID to add to each mutation
    :param mut_sep: Separator between mutations in the mutation string
    :returns: FEP variant string
    '''
    if is_wt(mut_str):
        return ''

    mutations = mut_str.strip().split(mut_sep)
    chain_mutations = [f'{chain_id}{FEP_CHAIN_SEP}{mut}' for mut in mutations]
    return FEP_MUTATION_SEP.join(chain_mutations)


def get_fep_variant_str_for_row(row, fc_chain_id, fcrn_chain_id, mut_sep='/'):
    '''
    Get the variant name for a row in the dataframe.
    '''
    fc_mut_str = row[FC_MUTATIONS]
    fcrn_mut_str = row[FCRN_MUTATIONS]
    fc_variant = get_fep_variant_str(fc_mut_str, fc_chain_id, mut_sep)
    fcrn_variant = get_fep_variant_str(fcrn_mut_str, fcrn_chain_id, mut_sep)
    if '' not in [fc_variant, fcrn_variant]:
        logger.warning(f'WARNING: Unsupported row with mutations in both Fc and FcRn: '
                       f'Fc={fc_mut_str}, FcRn={fcrn_mut_str}')
    variant_str = fc_variant + fcrn_variant  # fine as long as one is empty
    if variant_str == '':
        return WT
    return variant_str


def make_variant_column(df, fc_chain, fcrn_chain, mut_sep='/'):
    '''
    Add a column to the dataframe that combines the Fc/FcRn mutations and the
    chain IDs.
    '''
    df[VARIANT] = df.apply(
        lambda row: get_fep_variant_str_for_row(row, fc_chain, fcrn_chain,
                                                mut_sep=mut_sep),
        axis=1)


def prepare_variants_df(df, fc_chain, fcrn_chain):
    '''
    Prepare a dataframe of mean dG values for each variant.
    '''
    make_variant_column(df, fc_chain, fcrn_chain)
    agg_dict = {
        OFFSET_DG: 'mean',
        ERR_DG: rms,
        DDG: 'mean',
        ERR_DDG: rms,
    }
    gb = df.groupby(VARIANT)
    counts_df = gb.size().to_frame(name=NMEAS)
    agg_df = gb.agg(agg_dict)
    variants_df = counts_df.join(agg_df).reset_index()
    return variants_df


def parse_args(argv):
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument('-u', '--units',
                        dest='units',
                        type=str,
                        default='KI_nM',
                        choices=VALID_UNITS,
                        help='units for experimental binding affinity data in '
                        'in the user-provided CSV files; '
                        'the same units must be used for both pH values')
    parser.add_argument('--source-path', type=Path, default='source',
                        help='directory containing the source CSV files')
    parser.add_argument('--fc-chain-id', type=str, default='H',
                        help='chain ID for the Fc mutations')
    parser.add_argument('--fcrn-chain-id', type=str, default='A',
                        help='chain ID for the FcRn mutations')
    parser.add_argument('-v', dest='verbose', action='store_true',
                        help='increase verbosity')
    parser.add_argument('-debug', action='store_true',
                        help='print debugging info')
    return parser.parse_args(argv)


def main(argv=None):
    '''Main workflow, optionally callable like subprocess with list of args.'''
    args = parse_args(argv)

    # Adjust logging
    if args.debug:
        logger.setLevel(logging.DEBUG)
    elif args.verbose:
        logger.setLevel(logging.INFO)

    units = unit_convert.AffinityUnits[args.units]

    # Prepare the full experimental dataframe and write to a CSV file
    df = prepare_full_expt_df(args.source_path, units)
    df.to_csv(FULL_OUT_CSV, index=False, float_format='%0.4f')
    logger.debug(df.to_string(float_format='%.4f'))
    logger.warning(f'Wrote {len(df)} rows to {FULL_OUT_CSV}.')

    # Prepare the variants dataframe and write to a CSV file
    variants_df = prepare_variants_df(df, fc_chain=args.fc_chain_id,
                                      fcrn_chain=args.fcrn_chain_id)
    variants_df.to_csv(VARIANTS_OUT_CSV, float_format='%0.4f', index=False)
    logger.info(f'Wrote {len(variants_df)} rows to variants.csv.')

    logger.info('Done.')


if __name__ == '__main__':
    main()
