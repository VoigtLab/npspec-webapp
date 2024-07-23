from .loading_utils import *

name_converter = load_json('data/clean_name_to_display_name.json')

uniqueness_df = load_csv('data/computed_uniqueness_scores.csv', sep='\t')
uniqueness_df = get_family_labels(uniqueness_df)
uniqueness_df['hashed_smiles'] = [str(hash_smiles(s)) for s in uniqueness_df['smiles']]
uniqueness_df['display_name'] = uniqueness_df['name'].apply(lambda x: name_converter[x] if x in name_converter.keys() else x)

pred_spec_bars = load_json('data/all_done_bars.json')
for k in list(pred_spec_bars.keys()):
  pred_spec_bars[k.replace('_b3lyp','').split('/')[-1]] = pred_spec_bars[k]

pred_spec_df = load_spec_df('data/all_done_spectra.csv', transpose=False)
pred_spec_df.columns = [x.replace('_b3lyp','').split('/')[-1] for x in pred_spec_df.columns] 
pred_spec_df = pred_spec_df.loc[:, ~pred_spec_df.columns.duplicated()]

peaks_info = []
for c in uniqueness_df['name']:
    if c in pred_spec_df.columns:
        spec = pred_spec_df[c].values
        peak_info = list(get_peaks(spec/np.max(spec), pred_spec_df.index,  prominence = .1, width = 5))
        peaks_info.append(peak_info)
    else:
        peaks_info.append([])
uniqueness_df['peak_info'] = peaks_info

GRAPH_PATH = 'data/whole_metabolic_network_labeled.pkl'
met_graph = load_nx_graph(GRAPH_PATH)

ec_building_blocks = load_building_blocks('data/e_coli_metabolites_from_pathways.csv', sep='\t') 
ec_building_blocks.add('[Mg+2]')
