import os

import pandas as pd

tclin_data_file = 'output/drug_central_nodes.tsv' # Drug central nodes, with Tclin/Tchem/Tbio info
interaction_data = 'output/query_interactors_2nd.result'
output_file = 'output/query_interactors_2nd_druggable.tsv'
output_uniq_file = 'output/query_interactors_2nd_druggable_uniq.tsv'

tclin_data = pd.read_csv(tclin_data_file, sep="\t")
interaction_data = pd.read_csv(interaction_data, sep="\t")

# filter out stuff that isn't
tclin_data = tclin_data.loc[tclin_data['category'] == 'biolink:Protein']

interaction_data['?humanp2'] = interaction_data['?humanp2'].str.replace(r'<http://identifiers.org/uniprot/', r'UniProtKB:')
interaction_data['?humanp2'] = interaction_data['?humanp2'].str.replace(r'>', r'')

interaction_data['id'] = interaction_data['?humanp2']

result = pd.merge(interaction_data,
                 tclin_data[['id', 'TDL']],
                 on='id',
                 how='left')
result = result[result['TDL'].notna()]
result.to_csv(output_file, sep='\t', index=False)
os.system("cut -f7,8 {} | grep -v 'TDL' | sort | uniq | sort -k 1 > {}".format(output_file, output_uniq_file))
pass
