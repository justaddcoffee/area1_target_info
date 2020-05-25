import os
import tarfile
import wget

data_dir = 'data'
intact_path = '/'.join(['transformed', 'intact'])
sars_genes_path = '/'.join(['transformed', 'sars_cov_2_gene_annot'])

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

for path in [intact_path, sars_genes_path]:
    node_file = '/'.join([data_dir, path, 'nodes.tsv'])
    edge_file = '/'.join([data_dir, path, 'edges.tsv'])
    node_url = '/'.join(['http://kg-hub.berkeleybop.io', path, 'nodes.tsv'])
    edge_url = '/'.join(['http://kg-hub.berkeleybop.io', path, 'edges.tsv'])
    if not os.path.exists(node_file):
        wget.download(node_url, node_file)
    if not os.path.exists(edge_file):
        wget.download(edge_url, edge_file)
