import os
import tarfile
import wget

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
        def is_within_directory(directory, target):
            
            abs_directory = os.path.abspath(directory)
            abs_target = os.path.abspath(target)
        
            prefix = os.path.commonprefix([abs_directory, abs_target])
            
            return prefix == abs_directory
        
        def safe_extract(tar, path=".", members=None, *, numeric_owner=False):
        
            for member in tar.getmembers():
                member_path = os.path.join(path, member.name)
                if not is_within_directory(path, member_path):
                    raise Exception("Attempted Path Traversal in Tar File")
        
            tar.extractall(path, members, numeric_owner=numeric_owner) 
            
        
        safe_extract(tar, data_dir)

for path in [intact_path, sars_genes_path, drug_central]:
    node_file = '/'.join([data_dir, path, 'nodes.tsv'])
    edge_file = '/'.join([data_dir, path, 'edges.tsv'])
    node_url = '/'.join(['http://kg-hub.berkeleybop.io', path, 'nodes.tsv'])
    edge_url = '/'.join(['http://kg-hub.berkeleybop.io', path, 'edges.tsv'])
    if not os.path.exists(node_file):
        wget.download(node_url, node_file)
    if not os.path.exists(edge_file):
        wget.download(edge_url, edge_file)
