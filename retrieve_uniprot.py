import csv
import os
import re
import pandas
import requests
from pandas.errors import EmptyDataError
from tqdm import tqdm
import pandas as pd


def get_uniprot_data(data_dir: str, viral_data: pd.DataFrame) -> None:
    uniprot_base_url = 'https://www.uniprot.org/uniprot/'

    # gffs
    for this_id in tqdm(list(viral_data['protein ID'])):
        if re.match("UniProtKB", this_id):
            (prefix, uniprot_id) = this_id.split(":")
            url = uniprot_base_url + uniprot_id + ".gff"  # GFFs have glycosylation info
            s = requests.get(url).content
            outfile = os.path.join(data_dir, uniprot_id + ".gff")
            with open(outfile, 'wb') as f:
                f.write(s)

    # txt
    for this_id in tqdm(list(viral_data['protein ID'])):
        if re.match("UniProtKB", this_id):
            (prefix, uniprot_id) = this_id.split(":")
            url = uniprot_base_url + uniprot_id + ".txt"  # GFFs have glycosylation info
            s = requests.get(url).content
            outfile = os.path.join(data_dir, uniprot_id + ".txt")
            with open(outfile, 'wb') as f:
                f.write(s)

    # rdf
    for this_id in tqdm(list(viral_data['protein ID'])):
        if re.match("UniProtKB", this_id):
            (prefix, uniprot_id) = this_id.split(":")
            url = uniprot_base_url + uniprot_id + ".rdf"  # GFFs have glycosylation info
            s = requests.get(url).content
            outfile = os.path.join(data_dir, uniprot_id + ".rdf")
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

            matching_rows = dat[list(dat.iloc[:, col_idx].str.contains(gff_string))]
            if matching_rows.shape[0] > 0:
                new_col_list.append(True)
            else:
                new_col_list.append(False)
        else:
            new_col_list.append(False)
    viral_data[new_col_name] = new_col_list
    return None


def add_tclin_tchem_info(viral_data: pd.DataFrame,
                         tclin_col: str,
                         tchem_col: str,
                         new_col_name: str,
                         tclin_data_file) -> None:

    new_col_list: list = []

    d = {}
    with open(tclin_data_file) as fd:
        rd = csv.reader(fd, delimiter="\t", quotechar='"')
        header = next(rd)
        for row in rd:
            d[row[0]] = dict(zip(header, row))

    for this_id in list(viral_data['protein ID']):
        druggability = ''

        if re.match("UniProtKB", this_id) and this_id in d:
            entries = []
            if d[this_id]['tclin'] == 'True':
                entries.append('tclin')
            if d[this_id]['tchem'] == 'True':
                entries.append('tchem')
            druggability += "|".join(entries)

        new_col_list.append(druggability)

    viral_data[new_col_name] = new_col_list
    return None


def add_uniprot_txt_info(data_dir: str, viral_data: pd.DataFrame,
                         col_1_str: str, col_2_str: str, new_col_name: str,
                         file_ext: str = ".txt") -> None:
    new_col_list: list = []
    for this_id in list(viral_data['protein ID']):
        if re.match("UniProtKB", this_id):
            (prefix, uniprot_id) = this_id.split(":")
            infile = os.path.join(data_dir, uniprot_id + file_ext)

            # find line number where SQ starts for txt files, because pandas won't parse
            # them correctly in their entirety
            if infile.endswith(".txt"):
                sq_start = 0
                with open(infile) as fp:
                    for line in fp:
                        if line.startswith('SQ'):
                            break
                        sq_start += 1
                try:
                    dat = pd.read_csv(infile, sep="   ", comment='FT                   ',
                                      header=None, nrows=sq_start)
                except pandas.errors.ParserError:
                    pass
                except EmptyDataError:
                    new_col_list.append('')
                    continue

            else:
                try:
                    dat = pd.read_csv(infile, sep="   ", comment='#', header=None)
                except EmptyDataError:
                    new_col_list.append('')
                    continue

            matching_rows = dat[list(dat.iloc[:, 0].str.contains(col_1_str) & dat.iloc[:, 1].str.contains(col_2_str))]
            if matching_rows.shape[0] > 0:
                new_col_list.append(matching_rows.iloc[:,1].to_string(index=False))
            else:
                new_col_list.append('')
        else:
            new_col_list.append('')
    viral_data[new_col_name] = new_col_list
    return None
