from .loading_utils import *


uniqueness_df = load_csv('../spectranalysis/22Apr2023_computed_uniqueness_scores.csv')
uniqueness_df['hashed_smiles'] = [str(hash_smiles(s)) for s in uniqueness_df['smiles']]

pred_spec_bars = load_json('../biospectral/data/spectra/all_accessible_metabolites/7Apr2023_all_accesible_bars.json')
for k in list(pred_spec_bars.keys()):
  pred_spec_bars[k.replace('_b3lyp','').split('/')[-1]] = pred_spec_bars[k]

pred_spec_df = load_spec_df('../biospectral/data/spectra/all_accessible_metabolites/7Apr2023_all_done_spectra.csv', transpose=False)
pred_spec_df.columns = [x.replace('_b3lyp','').split('/')[-1] for x in pred_spec_df.columns] 
pred_spec_df = pred_spec_df.loc[:, ~pred_spec_df.columns.duplicated()]

GRAPH_PATH = '../noboi/noboi/construct_network_of_biosyntheses/08Mar2023_whole_metabolic_network_labeled.pkl'
met_graph = load_nx_graph(GRAPH_PATH)


ec_building_blocks = load_building_blocks('../noboi/noboi/parse_organism_dbs/e_coli/e_coli_metabolites_from_pathways.csv', sep='\t') 
ec_building_blocks.add('[Mg+2]')
