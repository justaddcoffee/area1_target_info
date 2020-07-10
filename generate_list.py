import os
import tarfile

import wget
from kg_covid_19.query import run_query
from pandas import read_csv
import urllib.parse
import urllib.request
import pandas as pd
import io
import requests
from tqdm import tqdm
from retrieve_uniprot import get_uniprot_data, add_uniprot_info, add_uniprot_txt_info, \
    add_tclin_tchem_info

print("Area 1 Target info")
data_dir = 'data'

transformed_dir = 'transformed'
intact_path = '/'.join([transformed_dir, 'intact'])
sars_genes_path = '/'.join([transformed_dir, 'sars_cov_2_gene_annot'])
drug_central = '/'.join([transformed_dir, 'drug_central'])

for this_dir in [data_dir,
                 os.path.join(data_dir, intact_path),
                 os.path.join(data_dir, sars_genes_path)]:
    os.makedirs(this_dir, exist_ok=True)

kg_nodes = os.path.join(data_dir, 'merged-kg_nodes.tsv')
kg_edges = os.path.join(data_dir, 'merged-kg_edges.tsv')
if not os.path.exists(kg_nodes) or not os.path.exists(kg_edges):
    kg_tar = os.path.join(data_dir, 'kg-covid-19.tar.gz')
    if not os.path.exists(kg_tar):
        wget.download('http://kg-hub.berkeleybop.io/kg-covid-19.tar.gz', data_dir)
    with tarfile.open(kg_tar) as tar:
        tar.extractall(data_dir)

for path in [intact_path, sars_genes_path, drug_central]:
    node_file = '/'.join([data_dir, path, 'nodes.tsv'])
    edge_file = '/'.join([data_dir, path, 'edges.tsv'])
    node_url = '/'.join(['http://kg-hub.berkeleybop.io', path, 'nodes.tsv'])
    edge_url = '/'.join(['http://kg-hub.berkeleybop.io', path, 'edges.tsv'])
    if not os.path.exists(node_file):
        wget.download(node_url, node_file)
    if not os.path.exists(edge_file):
        wget.download(edge_url, edge_file)

run_query('TargetCandidates', input_dir=data_dir, output_dir=data_dir)

target_info = read_csv('data/target_candidates.tsv', sep="\t")

pnnl_data = read_csv('data/pnnl/PNNLTargetList_2020_06_08 - TargetList_2020_06_05.tsv', sep="\t")
pnnl_data.rename(columns={'Uniprot': 'protein ID'}, inplace=True)

pnnl_data['protein ID'] = \
    pnnl_data['protein ID'].apply(lambda x: "{}{}".format('UniProtKB:', x))

target_info.drop(['confidence score','comments'], axis=1, inplace=True)

target_info = pd.concat([pnnl_data, target_info], axis=0, join='outer', ignore_index=False, keys=None,
          levels=None, names=None, verify_integrity=False, copy=True)

host_only = False
if host_only:
    indexNames = target_info[target_info['species'] == 'V'].index
    target_info.drop(indexNames, inplace=True)

# add druggability info
add_tclin_tchem_info(viral_data=target_info, tclin_col='tclin', tchem_col='tchem',
                     new_col_name='druggability',
                     tclin_data_file='data/transformed/drug_central/nodes.tsv')

# retrieve uniprot data for viral proteins
dl_uniprot_gff = False
uniprot_data_dir = os.path.join(data_dir, 'uniprot')
if dl_uniprot_gff:
    get_uniprot_data(uniprot_data_dir, target_info)

# query for glycosylation info
add_uniprot_info(uniprot_data_dir, target_info, 2, 'Glycosylation',
             'is_glycosylated')
add_uniprot_info(uniprot_data_dir, target_info, 2, 'Transmembrane',
             'has_transmembrane_domain')

# add PDB info, if any
add_uniprot_txt_info(uniprot_data_dir, target_info, 'DR', 'PDB;', 'pdb_data')
# add GO annotations
add_uniprot_txt_info(uniprot_data_dir, target_info, 'DR', 'GO\:', 'go_annotations')

target_info.drop_duplicates(subset=['protein ID', ], inplace=True)
target_info.sort_values(by=['species', 'protein ID'], inplace=True)

target_info.to_csv(path_or_buf=os.path.join(data_dir, 'target_candidates_plus.tsv'),
                   sep="\t",
                   index=False)


