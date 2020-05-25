import os
import re
import requests
from pandas.errors import EmptyDataError
from tqdm import tqdm
import pandas as pd


def get_uniprot_data(data_dir: str, viral_data: pd.DataFrame) -> None:
    uniprot_base_url = 'https://www.uniprot.org/uniprot/'
    for this_id in tqdm(list(viral_data['protein ID'])):
        if re.match("UniProtKB", this_id):
            (prefix, uniprot_id) = this_id.split(":")
            url = uniprot_base_url + uniprot_id + ".gff"  # GFFs have glycosylation info
            s = requests.get(url).content
            outfile = os.path.join(data_dir, uniprot_id + ".gff")
            with open(outfile, 'wb') as f:
                f.write(s)


def add_uniprot_info(data_dir: str, viral_data: pd.DataFrame, col_idx: int,
                     gff_string: str, new_col_name: str) -> None:
    new_col_list: list = []
    for this_id in list(viral_data['protein ID']):
        if re.match("UniProtKB", this_id):
            (prefix, uniprot_id) = this_id.split(":")
            infile = os.path.join(data_dir, uniprot_id + ".gff")
            try:
                dat = pd.read_csv(infile, sep="\t", comment='#', header=None)
            except EmptyDataError:
                new_col_list.append('Unknown')
                continue

            glyc_rows = dat[list(dat.iloc[:, col_idx].str.contains(gff_string))]
            if glyc_rows.shape[0] > 0:
                new_col_list.append(True)
            else:
                new_col_list.append(False)
        else:
            new_col_list.append(False)
    viral_data[new_col_name] = new_col_list
    return None
