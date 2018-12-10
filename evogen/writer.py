import json

def write_fasta_dict_to_json():
    pass

def write_fasta_dict_to_fasta(fasta_dict, path):
    """Writes a dictionary as a FASTA-formatted text file.

    Parameters
    ----------
    fasta_dict : dict
    path : str

    Returns
    -------
    int
        Number of sequences written to file.
    """
    cnt = 0
    with open(path, 'w') as f:
        for d in fasta_dict.values():
            # Error handling
            if 'id' not in d:
                raise KeyError('"id" key not found')
            if 'description' not in d:
                raise KeyError('"description" key not found')
            if 'sequence' not in d:
                raise KeyError('"sequence" key not found')

            if d['description']:
                print('>{} {}'.format(d['id'], d['description']), file=f)
            else:
                print('>{}'.format(d['id']), file=f)
            print(d['sequence'], file=f)
            cnt += 1
    return cnt
    