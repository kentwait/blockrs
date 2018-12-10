# -*- coding: utf-8 -*-
"""Helper functions for reading files.
"""
from collections import OrderedDict

def read_fasta_file_to_dict(path):
    """Reads a FASTA formatted text file into an ordered dictionary.

    Each entry is represented as a dictionary composed of
    the ID string <id>, description <description>, and
    the sequence <sequence>.

    Parameters
    ----------
    path : str

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
    return seq_d
