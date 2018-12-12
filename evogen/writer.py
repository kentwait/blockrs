# -*- coding: utf-8 -*-
"""Functions for writing files.
"""
import os
import pickle

def write_dict_to_fasta_file(d, path, line_width=None, description_parser=None):
    """Writes a dictionary of dictionaries into a FASTA-formatted text file.

    This function expects a dictionary of dictionaries where the inner
    dictionary has the following keys: `id`, `description`, and `sequence`.
    An exception is raised if one of these keys are not found.

    Parameters
    ----------
    d : dict
    path : str
    line_width : int
        Maximum number of characters per line of sequence.
    description_parser : function
        Function that uses the current sequence dictionary to create a
        string description. If None, the existing description value will be used.

    Returns
    -------
    int
        Number of sequences written to file.

    """
    # Check if path exists
    dirpath = os.path.dirname(path)
    if not os.path.exists(dirpath):
        raise Exception('{} does not exist'.format(dirpath))

    with open(path, 'w') as f:
        cnt = 0
        for k, seq_d in d.items():
            # Check if keys exist
            keys = seq_d.keys()
            for required_k in ['id', 'description', 'sequence']:
                if required_k not in keys:
                    raise KeyError('{} key not found in <{}> sequence'.format(required_k, k))

            if description_parser:
                seq_d['description'] = description_parser(seq_d)

            if seq_d['description']:
                print('>{} {}'.format(seq_d['id'], seq_d['description']),
                      file=f)
            else:
                print('>{}'.format(seq_d['id']), file=f)

            if line_width:
                s = [seq_d['sequence'][i:i+line_width]
                     for i in range(0, len(seq_d['sequence']), line_width)]
                print('\n'.join(s), file=f)
            else:
                print(seq_d['sequence'], file=f)

            cnt += 1
    return cnt


def write_dict_to_fasta_dir(d, dirpath, suffix='.aln',
                            filename_parser=None, description_parser=None,
                            line_width=None, verbose=False):
    """Writes dictionary of sequence dictionaries into FASTA files to be
    saved into a specified directory.

    Parameters
    ----------
    d : dict of dict
    dirpath : str
    suffix : str
        Appended to the end of the key to create the filename.
        When a filename_parser is specified, this parameter is ignored.
    filename_parser : function
        Function that uses the current sequence dictionary to create a filename.
        The expected input is a dictionary and outputs a string that will be
        used as the filename. If None, the key and the suffix will be the filename.
    description_parser : function
        Function that uses the current sequence dictionary to create a
        string description. If None, the existing description value will be used.
    line_width : int
        Maximum number of characters per line of sequence.
    verbose : bool

    Returns
    -------
    int
        Number of files written to the filesystem.

    """
    # Check if dirpath exists
    if not os.path.exists(dirpath):
        os.makedirs(dirpath)

    cnt = 0
    for k, seq_d in d.items():
        fname = str(k) + suffix
        if filename_parser:
            fname = filename_parser(k)

        path = os.path.join(dirpath, fname)
        c = write_dict_to_fasta_file(seq_d, path, line_width=line_width,
                                     description_parser=description_parser)

        if verbose:
            print(path, c)

        cnt += 1

    return cnt


def pickle_fasta_dir_dict(description, fasta_d, path, usage=None, **kwargs):
    """Pickles a nested dicts derived from a directory of FASTA files
    together with some metadata describing the data.

    This function encapsules the fasta_d within an outer dictionary
    that is separated into 'data', and 'meta' sections. The fasta_d dictionary
    is stored under 'data', while 'meta' contains information about the data
    such as a descriptive text under 'description', and optionally,
    details about data organization and schema in 'usage'.

    Other sections can be added to 'meta' by passing unlisted keyword arguments.

    By default, this encapsulation assumes two top-level keywords 'meta',
    and 'data'. Within 'meta', the 'description' key-value pair always exists,
    and 'usage' may sometimes exist. Other keys are user-specified and
    are not always guaranteed to exist and its value has no standard form.

    Parameters
    ----------
    description : str
        Describes the dataset
    fasta_d : dict
    usage : str, optional
        Details how the data is organized within the dictionary (schema)
    **kwargs
        Other key-value pairs that should be appended to the meta section

    """
    meta_list = [('description', description)]
    if usage:
        meta_list.append(('usage', usage))
    meta_list += list(kwargs.items())

    packaged_d = {'meta': dict(meta_list), 'data': fasta_d}
    with open(path, 'wb') as f:
        pickle.dump(packaged_d, f)
