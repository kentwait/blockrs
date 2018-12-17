# -*- coding: utf-8 -*-
"""Functions for reading files.
"""
from collections import OrderedDict
from hashlib import md5
import os
import pandas as pd

def read_fasta_file_to_dict(path, description_parser=None):
    """Reads a FASTA formatted text file into an ordered dictionary.

    Each entry is represented as a dictionary composed of
    the ID string <id>, description <description>, and
    the sequence <sequence>.

    Parameters
    ----------
    path : str
    description_parser : function
        Function that parses the description section of the ID line
        and adds additional keys to the dictionary. The expected input is
        the description as a string, and outputs a dictionary of key-value
        pairs parsed from the description. If None, the dictionary
        will only contain 'id', 'description', and 'sequence' key-values.

    Returns
    -------
    OrderedDict
        Dictionary where keys are the ID string and values are dictionaries
        composed of the ID string <id>, description <description>, and
        the sequence <sequence>.

    """
    name = ''
    description = ''
    seq = []
    seq_d = OrderedDict()
    with open(path, 'r') as f:
        for line in f.readlines():
            line = line.rstrip()
            if line.startswith('>'):
                if seq:
                    seq_string = ''.join(seq)
                    if name in seq_d.keys():
                        raise Exception('{} key already present'.format(name))
                    seq_d[name] = {'id': name,
                                   'description': description,
                                   'sequence': seq_string}
                    if description_parser:
                        seq_d[name] = dict(list(seq_d[name].items()) +
                                           list(description_parser(description).items()))
                    seq = []
                try:
                    name, description = line[1:].split(' ', 1)
                except ValueError:
                    name = line[1:]
                    description = ''
            else:
                seq += [line]

        if seq:
            seq_string = ''.join(seq)
            if name in seq_d.keys():
                raise Exception('{} key already present'.format(name))
            seq_d[name] = {'id': name,
                           'description': description,
                           'sequence': seq_string}
            if description_parser:
                seq_d[name] = dict(list(seq_d[name].items()) +
                                   list(description_parser(description).items()))
    return seq_d


def read_fasta_dir_to_dict(dirpath, suffix='.aln',
                           expected_count=None, filename_parser=None,
                           description_parser=None):
    """Reads a directory of FASTA files and stores data as
    nested dictionaries.

    The outer dictionary is used to identify contents from each
    FASTA file. Each FASTA file is then parsed as an ordered dictionary
    of sequences composed of the ID string <id>, description <description>,
    and the sequence <sequence>.

    Parameters
    ----------
    dirpath : str
    suffix : str
    expected_count : int
    filename_parser : function
        Function that transforms the filename into something for the key.
        The expected input is the filename (string) and outputs a
        string that will be used as the key. If None, the filename becomes
        the key.
    description_parser : function
        Function that parses the description section of the ID line
        and adds additional keys to the dictionary. The expected input is
        the description as a string, and outputs a dictionary of key-value
        pairs parsed from the description. If None, the dictionary
        will only contain 'id', 'description', and 'sequence' key-values.

    Returns
    -------
    dict
        Dictionary of OrderedDict where key is the filename or produced
        by the filename_parser function.

    """
    # Check if dirpath exists
    if not os.path.exists(dirpath):
        raise Exception('{} does not exist'.format(dirpath))
    else:
        if not os.path.isdir(dirpath):
            raise Exception('{} is not a directory'.format(dirpath))

    cnt = 0
    sequence_d = {}
    for fname in os.listdir(dirpath):
        if not fname.endswith(suffix):
            continue

        key = filename_parser(fname) if filename_parser else fname
        if key in sequence_d.keys():
            raise KeyError('{} already exists'.format(key))
        sequence_d[key] = read_fasta_file_to_dict(os.path.join(dirpath, fname),
                                                  description_parser=description_parser)

        cnt += 1

    if expected_count is not None:
        if expected_count != cnt:
            raise Exception('expected {} items, instead '
                            'processed {} file'.format(expected_count, cnt) +
                            's' if cnt != 1 else ''
                           )
    return sequence_d


def read_blast_table(path, sep='\t',
                     col_labels=['qaccver', 'saccver', 'pident', 'length',
                                 'mismatch', 'gapopen', 'qstart', 'qend',
                                 'sstart', 'send', 'evalue', 'bitscore', 'qlen',
                                 'slen', 'qcovs', 'sstrand']):
    """Reads BLAST results as a pandas DataFrame.

    Parameters
    ----------
    path : str
    col_labels : list of str
        Column labels of blast results.

    Returns
    -------
    DataFrame

    """
    df = pd.read_csv(path, sep=sep, header=None, index_col=None)
    df.columns = col_labels
    return df


def read_geneinfo_seq_file_to_dict(path, id_parser=lambda x: x.split('$')[-1]):
    """Reads a geneinfo sequence file into an ordered dictionary of sequences.

    Parameters
    ----------
    path : str
    id_parser : function
        Function that takes in the entire ID line and outputs the
        desired gene name/ID.
    
    Returns
    -------
    OrderedDict
        Each entry is identified by its gene name

    """
    gene_name = ''
    seq_name = ''
    seq = ''
    in_cds = False
    in_int = False
    total_cnt = None
    cnt = 0
    sequence_d = OrderedDict()
    with open(path, 'r') as f:
        for i, line in enumerate(f.readlines()):
            line = line.rstrip()
            # Get item count
            if i == 0:
                try:
                    total_cnt = int(line)
                except ValueError as e:
                    print('Expected item count, '
                          'instead encountered error {}'.format(e))

            if line.startswith('>'):
                if seq:
                    sequence_d[gene_name][seq_name] = {
                        'sequence': seq,
                        'baselen': len(seq),
                        'seqtype': 'cds' if seq_name == 'cds' else 'int',
                        'md5': md5(seq.encode('utf-8')).hexdigest(),
                    }
                    # Reset
                    in_cds = False
                    in_int = False
                    seq_name = ''
                    seq = ''
                _, gene_name = id_parser(line[1:]) if id_parser else line[1:]
                if gene_name in sequence_d.keys():
                    raise KeyError('{} key already present'.format(gene_name))
                sequence_d[gene_name] = OrderedDict({
                    'id': gene_name,
                    'description': 'sid={}, seq_type=cds'.format(line[1:]),
                })
                cnt += 1
            elif line.startswith('cod'):
                in_cds = True
                in_int = False
                seq_name = 'cds'
                sequence_d[gene_name][seq_name] = ''
            elif line.startswith('int'):
                in_cds = False
                in_int = True
                seq_name, _ = line.split(':')
                if seq_name in sequence_d[gene_name]:
                    raise KeyError('{} intron key already present in {}'.format(
                        seq_name, gene_name)
                    )
                sequence_d[gene_name][seq_name] = ''
            else:
                if in_cds:
                    seq = line.rstrip('$')
                    in_cds = False
                elif in_int:
                    seq = line.rstrip('$')
                    in_int = False
        if seq:
            sequence_d[gene_name][seq_name] = {
                'sequence': seq,
                'length': len(seq),
                'md5': md5(seq.encode('utf-8')).hexdigest(),
            }
        if total_cnt is not None and total_cnt != cnt:
            raise Exception('Expected {} entries, instead got {}'.format(total_cnt, cnt))
    return sequence_d    
