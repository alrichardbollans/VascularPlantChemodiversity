import os
import pickle

import pandas as pd
from phytochemMiner.extending_model_outputs import is_valid_inchikey, is_probably_valid_organic_smiles
from tqdm import tqdm

from configs import repo_path, WCVP_VERSION

# Register `pandas.progress_apply` and `pandas.Series.map_apply` with `tqdm`
# (can use `tqdm.gui.tqdm`, `tqdm.notebook.tqdm`, optional kwargs, etc.)
tqdm.pandas()
from phytochempy.compound_properties import simplify_inchi_key, fill_match_ids, standardise_SMILES
from phytochempy.wikidata_searches import get_wikidata_id_for_taxon, generate_wikidata_search_query, submit_query
from wcvpy.wcvp_download import get_all_taxa, wcvp_columns, wcvp_accepted_columns
from wcvpy.wcvp_name_matching import get_accepted_wcvp_info_from_ipni_ids_in_column, output_record_col_names, \
    get_accepted_info_from_names_in_column

# Vascular plant ID. Note may not contain angiosperms? I suspect as https://www.wikidata.org/wiki/Q14832431 doesn't have it as a parent.
# OK the wikidata taxonomy is a mess actually.
tracheophyte = 'Q27133'

compound_occurence_path = os.path.join(repo_path, 'temp_outputs', 'wikidata')
_temp_path = os.path.join(compound_occurence_path, 'temp')

# Might be best to do per family
family_pkl_file = os.path.join(repo_path, 'families.pkl')

wikidata_plantae_compounds_csv = os.path.join(compound_occurence_path, 'wikidata_compounds.csv')
# A version which keeps duplicates if found across different DOIs
wikidata_plantae_reference_data_csv = os.path.join(compound_occurence_path, 'wikidata_compound_doi.csv')


def get_all_families(rerun=False):
    """
    Fetches and organizes all families and their associated Wikidata identifiers.

    This function retrieves a list of all unique families from the WCVP taxonomy
    database. For each family, it generates corresponding Wikidata identifiers,
    either by re-computing the data or loading it from a pre-existing pickle
    file, depending on the `rerun` flag. The results are saved back to a pickle
    file for future use. If a family already exists in the dictionary, it is
    skipped.

    :param rerun: A flag that determines whether to regenerate all family-Wikidata
        relationships. If True, the function computes all data from scratch.
        Otherwise, it attempts to load saved data from a pickle file.
    :type rerun: bool
    :return: A dictionary containing family names as keys and their corresponding
        Wikidata identifier information as values.
    :rtype: dict
    """
    if rerun:
        out_dict = {}
    else:
        out_dict = pickle.load(open(family_pkl_file, 'rb'))

    wcvp = get_all_taxa(version=WCVP_VERSION)
    families = wcvp['family'].unique().tolist()

    for family in families:
        if family not in out_dict:
            fam_outputs = get_wikidata_id_for_taxon(family)
            if len(fam_outputs) > 0:
                out_dict[family] = fam_outputs

            with  open(family_pkl_file, 'wb') as pfile:
                pickle.dump(out_dict, pfile)
        else:
            print(f'{family} already in dict')


def tidy_wikidata_output(wikidata_results: pd.DataFrame):
    important_cols = ['organism_name', 'ipniID', 'structureLabel', 'structure_inchikey', 'structure_smiles',
                      'structure_cas', 'chembl_id', 'refDOI']
    wikidata_results = wikidata_results[important_cols]
    # To improve wikidata only include records with a doi reference
    wikidata_results = wikidata_results.dropna(subset=['structureLabel', 'organism_name', 'refDOI'], how='any')

    # rename to match other data sources
    wikidata_results = wikidata_results.rename(
        columns={'structureLabel': 'example_compound_name', 'structure_cas': 'CAS ID',
                 'structure_inchikey': 'InChIKey',
                 'structure_smiles': 'SMILES',
                 'ipniID': 'wikidata_ipniID'})

    tidy_final_output(wikidata_results, wikidata_plantae_compounds_csv, ipniid_col='wikidata_ipniID')
    tidy_final_output(wikidata_results, wikidata_plantae_reference_data_csv, ipniid_col='wikidata_ipniID',
                      for_paper_analysis=True)


def get_compounds_for_families():
    """
    Processes families and retrieves associated compounds by executing a query for
    each family ID. Results are saved as CSV files if they do not already exist
    in the specified temporary path.

    :raises FileNotFoundError: If the 'family_pkl_file' path does not exist or
        the file cannot be read.
    :raises KeyError: If the expected keys are not found in the 'fam_dict'.
    """
    fam_dict = pickle.load(open(family_pkl_file, 'rb'))
    for family in fam_dict:
        for id in fam_dict[family]:
            temp_output_csv = os.path.join(_temp_path, f'{family}_{id}.csv')
            if not os.path.exists(temp_output_csv):
                my_query = generate_wikidata_search_query(id, 1000000)
                submit_query(my_query, temp_output_csv, 1000000)
            else:
                print(f'{family} already in path')


def tidy_final_output(wikidata_results: pd.DataFrame, output_csv: str, ipniid_col=None,
                      for_paper_analysis: bool = False):
    '''

    Resolve names in wikidata query output and tidy.

    :return:
    '''
    for smiles in wikidata_results['SMILES'].dropna().unique():
        try:
            assert is_probably_valid_organic_smiles(smiles)
        except AssertionError:
            print(smiles)
    for ink in wikidata_results['InChIKey'].dropna().unique():
        assert is_valid_inchikey(ink)

    wikidata_results['Standard_SMILES'] = wikidata_results['SMILES'].progress_apply(standardise_SMILES)
    for c_id in ['InChIKey', 'Standard_SMILES', 'CAS ID']:
        wikidata_results = fill_match_ids(wikidata_results, c_id)
    wikidata_results['InChIKey_simp'] = wikidata_results['InChIKey'].apply(simplify_inchi_key)

    wikidata_results = wikidata_results.dropna(subset=['InChIKey_simp', 'Standard_SMILES', 'CAS ID'], how='all')

    if for_paper_analysis:
        wikidata_results = wikidata_results.drop_duplicates(subset=['organism_name', 'example_compound_name', 'refDOI'],
                                                            keep='first')
    else:
        wikidata_results = wikidata_results.drop_duplicates(
            subset=['organism_name', 'example_compound_name'],
            keep='first')
        if 'refDOI' in wikidata_results.columns:
            wikidata_results.drop(columns=['refDOI'])

    if ipniid_col is not None:
        all_taxa = get_all_taxa(version=WCVP_VERSION)
        ipni_matches = get_accepted_wcvp_info_from_ipni_ids_in_column(wikidata_results,
                                                                      'wikidata_ipniID',
                                                                      all_taxa)[
            wikidata_results.columns.tolist() + output_record_col_names]
        unmatched = ipni_matches[ipni_matches[wcvp_columns['wcvp_id']].isna()][wikidata_results.columns]

        ipni_matched = ipni_matches[~ipni_matches[wcvp_columns['wcvp_id']].isna()]
        ipni_matched['matched_by'] = 'ipni_id'
        name_matched = get_accepted_info_from_names_in_column(unmatched, 'organism_name', wcvp_version=WCVP_VERSION,
                                                              all_taxa=all_taxa)

        acc_df = pd.concat([ipni_matched, name_matched])
    else:
        acc_df = get_accepted_info_from_names_in_column(wikidata_results, 'organism_name', wcvp_version=WCVP_VERSION)

    outcols = ['example_compound_name', 'InChIKey', 'InChIKey_simp', 'SMILES', 'Standard_SMILES', 'CAS ID',
               'organism_name', wcvp_accepted_columns['name'],
               wcvp_accepted_columns['name_w_author'],
               wcvp_accepted_columns['species'], 'accepted_family']
    if 'refDOI' in acc_df.columns:
        outcols.append('refDOI')
    acc_df = acc_df[outcols]

    acc_df = acc_df.dropna(subset=[wcvp_accepted_columns['name']], how='any')
    acc_df = acc_df.sort_values(by=[wcvp_accepted_columns['name'], 'example_compound_name'])
    acc_df.to_csv(output_csv)
    acc_df.describe(include='all').to_csv(output_csv.replace('.csv', '_summary.csv'))


def tidy_outputs():
    all_family_compounds = pd.DataFrame()
    import glob
    path = _temp_path + "/*.csv"
    for fname in glob.glob(path):
        # print(fname)
        df = pd.read_csv(fname)
        all_family_compounds = pd.concat([all_family_compounds, df])

    tidy_wikidata_output(all_family_compounds)


if __name__ == '__main__':
    get_all_families()
    get_compounds_for_families()
    tidy_outputs()
