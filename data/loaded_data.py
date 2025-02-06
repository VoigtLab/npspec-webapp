from .loading_utils import *

name_converter = load_json('data/clean_name_to_display_name.json')
reverse_name_converter = {v:k for k,v in name_converter.items()}

uniqueness_df = load_csv('data/computed_uniqueness_scores.csv', sep='\t')
uniqueness_df = get_family_labels(uniqueness_df)
uniqueness_df['hashed_smiles'] = [str(hash_smiles(s)) for s in uniqueness_df['smiles']]
uniqueness_df['display_name'] = uniqueness_df['name'].apply(lambda x: name_converter[x] if x in name_converter.keys() else x)

pred_spec_bars = load_json('data/all_done_bars.json')
for k in list(pred_spec_bars.keys()):
  pred_spec_bars[k.replace('_b3lyp','').split('/')[-1]] = pred_spec_bars[k]

GRAPH_PATH = 'data/whole_metabolic_network_labeled.pkl'
met_graph = load_nx_graph(GRAPH_PATH)

ec_building_blocks = load_building_blocks('data/e_coli_metabolites_from_pathways.csv', sep='\t') 
ec_building_blocks.add('[Mg+2]')

reaction_df = pd.read_csv('data/all_reactions.csv', sep='\t')

for c in ['level_0', 'Unnamed: 0.1', 'Unnamed: 0']:
    if c in reaction_df.columns:
        reaction_df = reaction_df.drop(columns=[c])
