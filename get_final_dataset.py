import os

import pandas as pd
from phytochempy.chemical_diversity_metrics import get_pathway_based_diversity_measures, calculate_FAD_measures, \
    compile_rarified_calculations
from phytochempy.compound_properties import get_npclassifier_classes_from_df, get_npclassifier_pathway_columns_in_df, \
    read_manual_npclassifier_input, \
    NP_PATHWAYS
from phytochempy.data_compilation_utilities import tidy_final_dataset

from configs import repo_path

COMPOUND_ID_COL = 'Standard_SMILES'
COMPOUND_GROUP_COL = 'accepted_species'

compound_occurrence_csv = os.path.join(repo_path, 'outputs', 'compound_occurrences.csv')
species_chemodiversity_csv = os.path.join(repo_path, 'outputs', 'species_chemodiversity.csv')
species_chemodiversity_transformed_csv = os.path.join(repo_path, 'outputs', 'species_chemodiversity_transformed.csv')
max_number_of_pathways = len(NP_PATHWAYS)
FAD_INDICES = ['FAD', 'MFAD', 'APWD']
PATHWAY_INDICES = ['H', 'Hbc', 'G', 'J']
CHEMODIV_METRICS = PATHWAY_INDICES + FAD_INDICES


def get_final_occcurence_data():
    phytochem_downloads_path = os.path.join(repo_path, 'temp_outputs')
    knapsack_plantae_compounds_csv = os.path.join(phytochem_downloads_path, 'knapsack_data', 'knapsack_compounds.csv')
    wikidata_plantae_compounds_csv = os.path.join(phytochem_downloads_path, 'wikidata', 'wikidata_compounds.csv')
    tidy_wiki_data = pd.read_csv(wikidata_plantae_compounds_csv, index_col=0)
    tidy_knapsack_data = pd.read_csv(knapsack_plantae_compounds_csv, index_col=0)
    all_compounds_in_taxa = pd.concat([tidy_wiki_data, tidy_knapsack_data], axis=0, ignore_index=True)
    _temp_outputs_path = os.path.join(repo_path, 'temp_outputs')

    # Upload this using https://gnps.ucsd.edu/ProteoSAFe/index.jsp?params=%7B%22workflow%22:%22NPCLASSIFIER%22%7D
    # need a 'smiles' column in input file: https://github.com/CCMS-UCSD/GNPS_Workflows/blob/master/npclassifier/tools/npclassifier/bin/npclassify.py
    for_manual_upload_to_gnps = all_compounds_in_taxa[['Standard_SMILES']]
    for_manual_upload_to_gnps = for_manual_upload_to_gnps.dropna(subset=['Standard_SMILES']).drop_duplicates(
        subset=['Standard_SMILES']).rename(
        columns={'Standard_SMILES': 'smiles'})
    for_manual_upload_to_gnps.to_csv(os.path.join(_temp_outputs_path, 'compounds_for_manual_upload_to_gnps.csv'),
                                     index=False)

    # Check job using https://gnps.ucsd.edu/ProteoSAFe/jobs.jsp#%7B%22table_sort_history%22%3A%22create_time_millis_dsc%22%7D
    # Save manual file to _temp_outputs_path beginning with npclassifierinfo_manual_
    # Note I had to fix the reading function changing line 179, 180 to
    # np_classif_results['col_split'] = np_classif_results[col].str.split(',')
    # max_splits = int(np_classif_results['col_split'].str.len().max())
    assert os.path.isfile(os.path.join(_temp_outputs_path, 'npclassifierinfo_manual_gnps_results.tsv'))

    ## Add extra information related to the compound properties

    with_npclass_classes = get_npclassifier_classes_from_df(all_compounds_in_taxa, 'Standard_SMILES',
                                                            _temp_outputs_path)
    with_npclass_classes.to_csv(os.path.join(_temp_outputs_path, 'npclassifier.csv'))
    pway_cols = get_npclassifier_pathway_columns_in_df(with_npclass_classes)

    with_npclass_classes['pairs'] = with_npclass_classes[COMPOUND_GROUP_COL] + with_npclass_classes[COMPOUND_ID_COL]
    ### Then tidy and output final dataset (remove duplicates)
    tidy_final_dataset(with_npclass_classes, compound_occurrence_csv, COMPOUND_ID_COL, 'accepted_species')

    summary = pd.read_csv(compound_occurrence_csv, index_col=0)
    summary.describe(include='all').to_csv(
        os.path.join('outputs', 'compound_occurrences_summary.csv'))


def transform_compiled_data(compiled_data: pd.DataFrame):
    from sklearn.preprocessing import PowerTransformer
    transformer = PowerTransformer(method='yeo-johnson')
    compiled_data = compiled_data.set_index(COMPOUND_GROUP_COL, drop=True)
    transformed_data = transformer.fit_transform(
        compiled_data[CHEMODIV_METRICS + ['GroupSize_FAD', 'GroupSize_Pathways']])
    # Convert the transformed data back into a DataFrame
    df_transformed = pd.DataFrame(transformed_data,
                                  columns=CHEMODIV_METRICS + ['GroupSize_FAD', 'GroupSize_Pathways'])
    df_transformed[COMPOUND_GROUP_COL] = compiled_data.index
    df_transformed = df_transformed[
        [COMPOUND_GROUP_COL] + CHEMODIV_METRICS + ['GroupSize_FAD', 'GroupSize_Pathways']]
    df_transformed.to_csv(species_chemodiversity_transformed_csv)


def add_diversity_data():
    working_data = pd.read_csv(compound_occurrence_csv, index_col=0)

    # Remove duplicates and remove species with fewer than 7 compounds
    working_data = working_data.drop_duplicates(subset=[COMPOUND_GROUP_COL, COMPOUND_ID_COL], keep='first')
    species_counts = working_data.groupby(COMPOUND_GROUP_COL)[COMPOUND_ID_COL].transform('count')
    working_data = working_data[species_counts >= max_number_of_pathways]

    abundance_diversity = get_pathway_based_diversity_measures(working_data, COMPOUND_GROUP_COL, COMPOUND_ID_COL)

    FAD_measures = calculate_FAD_measures(working_data, compound_grouping=COMPOUND_GROUP_COL)

    ## For this many groups, rarifying is extremely computationally expensive.
    ## As seen in richard-bollans et al 2025, it is likely that rarified FAD, MFAD and APWD resemble unrarified APWD.
    ## And that rarified H, Hbc and G resemble unrarified versions.
    # fad_rare, pathway_rare = compile_rarified_calculations(working_data, COMPOUND_GROUP_COL,
    #                                                        max_number_of_pathways, COMPOUND_ID_COL,
    #                                                        1000)

    compiled_data = pd.merge(abundance_diversity, FAD_measures, how='left', on=COMPOUND_GROUP_COL,
                             validate='one_to_one')

    for g in working_data[COMPOUND_GROUP_COL].tolist():
        if g not in compiled_data[COMPOUND_GROUP_COL].values:
            print(f'Warning: {g} not in compiled data for some reason')

    compiled_data = compiled_data.dropna(subset=CHEMODIV_METRICS, how='all')
    compiled_data.to_csv(species_chemodiversity_csv)
    compiled_data.describe(include='all').to_csv(
        os.path.join('outputs', f'species_chemodiversity_summary.csv'))

    transform_compiled_data(compiled_data)


if __name__ == '__main__':
    get_final_occcurence_data()
    add_diversity_data()
