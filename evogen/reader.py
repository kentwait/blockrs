# -*- coding: utf-8 -*-
"""Helper functions for reading files.
"""
from collections import OrderedDict
import os

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
