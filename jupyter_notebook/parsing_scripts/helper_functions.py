import csv
import pandas
from Bio import SeqIO
from pysradb import SRAdb, download_sradb_file

db = SRAdb('SRAmetadb.sqlite')

def get_biosample_title(query_accession, db):
    """Retrieve biosample record title for the given query accession number (ID beginning with SRR, DRR, or ERR.)"""
    run_accession=query_accession.split(".")[0]
    df = db.query('select experiment_title from sra where run_accession="{run_accession_id}";'.format(run_accession_id=run_accession))
    return df

def get_biosample_attribute(query_accession, db):
    """Retrieve biosample record attributes for the given query accession number (ID beginning with SRR, DRR, or ERR.)"""
    run_accession=query_accession.split(".")[0]
    df = db.query('select sample.sample_attribute from sample INNER JOIN sra ON sra.sample_accession=sample.sample_accession WHERE sra.run_accession="{run_accession_id}";'.format(run_accession_id=run_accession))
    return df

def add_biosample_title(df, query_accession, db):
    df['biosample_title'] = df.apply(lambda x: get_biosample_title(x, db).get("experiment_title")[0], axis=1)
    return df

def add_biosample_attribute(df, query_accession, db):
    df['biosample_attr'] = df.apply(lambda x: get_biosample_title(x, db).get("sample_attribute")[0], axis=1)
    return df

def get_query_taxid(query_id, query_records, tax_dict):
    taxid = tax_dict[query_id]
    name = [rec.description for rec in query_records if query_id == rec.id][0]
    return taxid, name[name.index(' ') + 1:name.index(' ITS')]


def get_taxon_dict(df, query_fasta_file):
    with open(query_fasta_file, 'r') as query_fasta:
        names = list(SeqIO.parse(query_fasta, 'fasta'))
    with open('../nr.nucl_gb.acc2taxid.txt', 'r') as csvfile:
        rows = csv.reader(csvfile, delimiter='\t')
        taxon_dict = {row[1]: row[2] for row in rows}
    accs = list(set(df['query']))
    return {acc: get_query_taxid(acc, names, taxon_dict) for acc in accs}


def add_taxon_info(df, query_fasta_file):
    taxid_dict = get_taxon_dict(df, query_fasta_file)
    df['taxid'] = df.apply(lambda x: taxid_dict[x['query']][0], axis=1)
    df['scientific_name'] = df.apply(lambda x: taxid_dict[x['query']][1], axis=1)
    return df
