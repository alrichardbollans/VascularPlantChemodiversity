import os
import pickle

import numpy as np
import pandas as pd
from phytochemMiner.extending_model_outputs import resolve_name_to_smiles, resolve_name_to_inchi
from phytochempy.compound_properties import add_CAS_ID_translations_to_df
from phytochempy.knapsack_searches import get_knapsack_compounds_in_family
from tqdm import tqdm
import sys

from configs import repo_path
from get_wikidata import family_pkl_file, tidy_final_output

sys.path.append('..')
# Register `pandas.progress_apply` and `pandas.Series.map_apply` with `tqdm`
# (can use `tqdm.gui.tqdm`, `tqdm.notebook.tqdm`, optional kwargs, etc.)
tqdm.pandas()

knapsack_data_path = os.path.join(repo_path, 'temp_outputs', 'knapsack_data')
_temp_path = os.path.join(knapsack_data_path, 'temp')
knapsack_plantae_compounds_csv = os.path.join(knapsack_data_path, 'knapsack_compounds.csv')


def get_knapsack_data_for_each_family():
    fam_dict = pickle.load(open(family_pkl_file, 'rb'))
    families_of_interest = list(set(fam_dict.keys()))
    for i in tqdm(range(len(families_of_interest)), desc=f"Searching families…"):
        family = families_of_interest[i]
        out_csv = os.path.join(_temp_path, f'{family}.csv')
        if not os.path.exists(out_csv):
            get_knapsack_compounds_in_family(family, out_csv)
        else:
            print(f'{family} already exists')


def tidy_knapsack_output(knapsack_results: pd.DataFrame, output_csv: str):
    important_cols = ['CAS ID', 'example_compound_name', 'Organism']
    knapsack_results = knapsack_results[important_cols]
    knapsack_results = knapsack_results.dropna(subset=['example_compound_name', 'Organism'], how='any')
    knapsack_results = knapsack_results.drop_duplicates(subset=important_cols,
                                                        keep='first')

    knapsack_results = knapsack_results.rename(
        columns={'Organism': 'organism_name'})
    knapsack_results = add_CAS_ID_translations_to_df(knapsack_results, 'CAS ID',
                                                     os.path.join(knapsack_data_path, 'cirpycache'))
    problems = knapsack_results[knapsack_results['InChIKey'] == '']
    assert len(problems) == 0
    unresolved_inchikey = knapsack_results[knapsack_results['InChIKey'].isna()][
        ['example_compound_name']].drop_duplicates(subset=['example_compound_name'])
    unresolved_inchikey['InChIKey_from_name'] = unresolved_inchikey['example_compound_name'].progress_apply(
        resolve_name_to_inchi)
    knapsack_results = pd.merge(knapsack_results, unresolved_inchikey, on='example_compound_name', how='left')
    knapsack_results['InChIKey'] = np.where(knapsack_results['InChIKey'].isna(),
                                            knapsack_results['InChIKey_from_name'], knapsack_results['InChIKey'])
    knapsack_results = knapsack_results.drop(columns=['InChIKey_from_name'])

    unresolved_smiles = knapsack_results[knapsack_results['SMILES'].isna()][['example_compound_name']].drop_duplicates(
        subset=['example_compound_name'])
    unresolved_smiles['SMILES_from_name'] = unresolved_smiles['example_compound_name'].progress_apply(
        resolve_name_to_smiles)
    knapsack_results = pd.merge(knapsack_results, unresolved_smiles, on='example_compound_name', how='left')
    knapsack_results['SMILES'] = np.where(knapsack_results['SMILES'].isna(),
                                          knapsack_results['SMILES_from_name'], knapsack_results['SMILES'])
    knapsack_results = knapsack_results.drop(columns=['SMILES_from_name'])

    tidy_final_output(knapsack_results, output_csv)


def compile_family_data():
    big_df = pd.DataFrame()
    for file in os.listdir(_temp_path):
        out_csv = os.path.join(_temp_path, file)
        df = pd.read_csv(out_csv)
        big_df = pd.concat([big_df, df])
    tidy_knapsack_output(big_df, knapsack_plantae_compounds_csv)


if __name__ == '__main__':
    get_knapsack_data_for_each_family()
    compile_family_data()
