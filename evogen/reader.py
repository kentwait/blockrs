# -*- coding: utf-8 -*-
"""Functions for reading files.
"""
from collections import OrderedDict
from hashlib import md5
from copy import deepcopy
import os
import re
import pandas as pd
import numpy as np
from evogen.block import combine_exon_blocks

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


def read_genpos_file_to_dicts(path, convert_to_zero_based_index=True):
    # Flags
    geneinfo_finished = False
    transcript_metadata_flag = False
    tr_exon_num_finished = False
    cds_metadata_flag = False
    cds_exon_num_finished = False
    cds_exon_coords_finished = False
    intron_metadata_flag = False
    intron_num_finished = False
    intron_coords_finished = False
    # temp dicts
    tr_meta = {}
    cds_data = {}
    int_data = {}
    # temp vars
    geneinfo_id = None

    # Read to ordered dicts
    tr_meta_d = OrderedDict()
    cds_d = OrderedDict()
    intron_d = OrderedDict()
    with open(path, 'r') as f:
        for line in f.readlines():
            line = line.rstrip()
            # Skip comment lines
            if line.startswith('/*'):
                continue
            # ID line
            if line.startswith('>'):
                try:
                    geneinfo_id = int(line[1:])
                    tr_meta = {
                        'geneinfo_id': geneinfo_id
                    }
                except ValueError as e:
                    raise e
                geneinfo_finished = True
            # Chromosome metadata
            elif geneinfo_finished:
                chr_pos, chr_scaffold, orientation = line.split('_')
                tr_meta['chromosome'] = chr_pos
                # dont convert scaffold number to 0-index
                tr_meta['scaffold'] = int(chr_scaffold)
                tr_meta['forward'] = True if orientation == '+' else False
                geneinfo_finished = False
            # Transcript metadata
            elif line.startswith('tr_dat'):
                transcript_metadata_flag = True
            elif transcript_metadata_flag and \
                line.startswith('tr_start_site_in_genome'):
                try:
                    start_loc = re.search(r'\:\s*(\d+)', line).group(1)
                    tr_meta['genome_start'] = int(start_loc)
                except AttributeError:
                    raise Exception('Failed to match '
                                    'tr_start_site_in_genome field: {}'.format(line))
            elif transcript_metadata_flag and \
                line.startswith('tr_end_site_in_genome'):
                try:
                    stop_loc = re.search(r'\:\s*(\d+)', line).group(1)
                    tr_meta['genome_stop'] = int(stop_loc)
                except AttributeError:
                    raise Exception('Failed to match '
                                    'tr_end_site_in_genome field: {}'.format(line))
            elif transcript_metadata_flag and \
                line.startswith('tr_length'):
                try:
                    baselen = re.search(r'\:\s*(\d+)', line).group(1)
                    tr_meta['baselen'] = int(baselen)
                except AttributeError:
                    raise Exception('Failed to match for '
                                    'tr_length field: {}'.format(line))
            elif transcript_metadata_flag and \
                line.startswith('tr_exon_num'):
                try:
                    exon_num = re.search(r'\:\s*(\d+)', line).group(1)
                    tr_meta['exon_count'] = int(exon_num)
                except AttributeError:
                    raise Exception('Failed to match for '
                                    'tr_exon_num field: {}'.format(line))
                tr_exon_num_finished = True
            elif transcript_metadata_flag and tr_exon_num_finished:
                try:
                    str_blocks = re.findall(r'\d+\.\.\d+', line)
                    slice_blocks = \
                        [slice(*tuple(map(int, s.split('..'))))
                         for s in str_blocks]
                    tr_meta['exon_blocks'] = slice_blocks
                except AttributeError:
                    raise Exception('Failed to match for'
                                    'exon coordinates: {}'.format(line))
                tr_exon_num_finished = False
                # Add to ordered dict
                tr_meta_d[geneinfo_id] = deepcopy(tr_meta)
                tr_meta = {}
                transcript_metadata_flag = False
                
            # CDS metadata
            elif line.startswith('cds_dat'):
                cds_metadata_flag = True
                transcript_metadata_flag = False  # redundant
                cds_data = {
                    'geneinfo_id': geneinfo_id,
                }
            elif cds_metadata_flag and \
                line.startswith('cds_len'):
                try:
                    baselen = re.search(r'\:\s*(\d+)', line).group(1)
                    cds_data['cds_len'] = int(baselen)
                except AttributeError:
                    raise Exception('Failed to match for'
                                    'cds_len field: {}'.format(line))
            elif cds_metadata_flag and \
                line.startswith('start_loc_cds'):
                try:
                    start_loc = re.search(r'\:\s*(\d+)', line).group(1)
                    cds_data['genome_start'] = int(start_loc)
                except AttributeError:
                    raise Exception('Failed to match for'
                                    'start_loc_cds: {}'.format(line))
            elif cds_metadata_flag and \
                line.startswith('end_loc_cds'):
                try:
                    stop_loc = re.search(r'\:\s*(\d+)', line).group(1)
                    cds_data['genome_stop'] = int(stop_loc)
                except AttributeError:
                    raise Exception('Failed to match for'
                                    'end_loc_cds: {}'.format(line))
            elif cds_metadata_flag and \
                line.startswith('cds_exon_num'):
                try:
                    exon_num = re.search(r'\:\s*(\d+)', line).group(1)
                    cds_data['exon_count'] = int(exon_num)
                except AttributeError:
                    raise Exception('Failed to match for'
                                    'cds_exon_num: {}'.format(line))
                cds_exon_num_finished = True
            elif cds_metadata_flag and \
                cds_exon_num_finished:
                try:
                    str_blocks = re.findall(r'\d+\.\.\d+', line)
                    slice_blocks = \
                        [slice(*tuple(map(int, s.split('..')))) 
                         for s in str_blocks]
                    cds_data['exon_blocks'] = slice_blocks
                except AttributeError:
                    raise Exception('Failed to match for'
                                    'cds exon coordinates: {}'.format(line))
                cds_exon_num_finished = False
                cds_exon_coords_finished = True
            elif cds_metadata_flag and \
                cds_exon_coords_finished:
                cds_data['sequence'] = line
                cds_exon_coords_finished = False
                # Add to ordered dict
                if geneinfo_id in cds_d.keys():
                    raise KeyError('{} already exists in cds dictionary'.format(geneinfo_id))
                cds_d[geneinfo_id] = deepcopy(cds_data)
                cds_data = {}
                cds_metadata_flag = False
            # Intron metadata
            elif line.startswith('intron_dat'):
                intron_metadata_flag = True
                cds_metadata_flag = False  # redundant
            elif intron_metadata_flag and line.startswith('intron #'):
                try:
                    int_num = re.search(r'\:\s*(\d+)', line).group(1)
                    tr_meta_d[geneinfo_id]['intron_count'] = int(int_num)  # add to tr_meta_d
                except AttributeError:
                    raise Exception('Failed to match for'
                                    'intron_num field: {}'.format(line))
                intron_num_finished = True
            elif intron_metadata_flag and intron_num_finished:
                if line.startswith('-'):
                    geneinfo_id = ''
                    transcript_metadata_flag = False
                    cds_metadata_flag = False
                    intron_metadata_flag = False
                    continue
                if not intron_coords_finished:
                    try:
                        match = list(map(int, re.findall(r'\d+', line)))
                        int_data = {
                            'geneinfo_id': geneinfo_id,
                            'int_id': match[0],
                            # 'coding_type_id': match[1],  # 1: utr, 2: cds
                            # 'coding_int_id': match[2],   # 4xx: utr, 5xx: cds
                            'cds_int_id': None if match[1] == 1 else (match[2] - 500),
                            'transcribed': False if match[1] == 1 else True,
                            'baselen': match[3],
                            # Convert to 0-based
                            'from_tss_start': match[4],
                            'from_tss_stop': match[5],
                        }
                        
                    except Exception as e:
                        print(e, line, sep='\n')
                    intron_coords_finished = True
                else:
                    int_data['sequence'] = line
                    # Add to transcript_metadata
                    key = (geneinfo_id, int_data['int_id'])
                    if key in intron_d.keys():
                        raise KeyError('{} already exists in intron dictionary'.format(
                            str(key)))
                    intron_d[key] = deepcopy(int_data)
                    int_data = {}
                    intron_coords_finished = False
    if convert_to_zero_based_index:
        for gid, d in tr_meta_d.items():
            if tr_meta_d[gid]['genome_start'] < tr_meta_d[gid]['genome_stop']:
                tr_meta_d[gid]['genome_start'] -= 1
            else:
                tr_meta_d[gid]['genome_stop'] -= 1
            tr_meta_d[gid]['exon_blocks'] = [slice(s.start - 1, s.stop)
                                            for s in d['exon_blocks']]
        for gid, d in cds_d.items():
            if tr_meta_d[gid]['genome_start'] < tr_meta_d[gid]['genome_stop']:
                cds_d[gid]['genome_start'] -= 1
            else:
                cds_d[gid]['genome_stop'] -= 1
            cds_d[gid]['exon_blocks'] = [slice(s.start - 1, s.stop)
                                        for s in d['exon_blocks']]
        for int_id in intron_d.keys():
            intron_d[int_id]['from_tss_start'] -= 1

    return tr_meta_d, cds_d, intron_d

def read_genpos_file_to_dataframes(path, convert_to_zero_based_index=True):
    tr_meta_d, cds_d, intron_d = \
        read_genpos_file_to_dicts(path,convert_to_zero_based_index)

    # create transcript metadata df
    # transcript metadata
    tr_cols = ['geneinfo_id', 'chromosome', 'scaffold',
               'forward', 'genome_start', 'genome_stop',
               'baselen', 'exon_count', 'intron_count']
    tr_df = pd.DataFrame.from_dict(tr_meta_d,
                                   orient='index',
                                   columns=tr_cols)
    # add 'exon_blocks' column
    # tr_df['exon_blocks'] = [';'.join(['{}:{}'.format(s.start, s.stop)
    #                                   for s in d['exon_blocks']])
    #                         for k, d in tr_meta_d.items()]
    # Set index to geneinfo_id
    tr_df = tr_df.set_index('geneinfo_id')

    # exon blocks df
    # can reconstruct cds_df from this dataframe
    tr_exon_blocks_list = []
    for k, d in cds_d.items():
        all_exon_blocks, exon_ids = \
            combine_exon_blocks(tr_meta_d[k]['exon_blocks'], d['exon_blocks'])
        assert len(all_exon_blocks) == len(exon_ids)
        from_cds_start = 0
        cds_cnt = 1  # 1-based index
        for i, s in enumerate(all_exon_blocks):
            baselen = s.stop - s.start
            if tr_meta_d[k]['forward']:
                gstart = tr_meta_d[k]['genome_start'] + s.start
                gstop = tr_meta_d[k]['genome_start'] + s.stop
            else:
                gstart = tr_meta_d[k]['genome_start'] - s.start
                gstop = tr_meta_d[k]['genome_start'] - s.stop
            if s.step == 1:
                tr_exon_blocks_list.append({
                    'geneinfo_id': k,
                    'exon_id': exon_ids[i],
                    'transcribed': True,
                    # 'cds_exon_id': cds_cnt,
                    'baselen': baselen,
                    'from_tss_start': s.start,
                    'from_tss_stop': s.stop,
                    'from_cds_start': from_cds_start,
                    'from_cds_stop': from_cds_start + baselen,
                    'genome_start': gstart,
                    'genome_stop': gstop,
                    'sequence': d['sequence'][from_cds_start:
                                              from_cds_start + baselen]
                })
                from_cds_start += baselen
                cds_cnt += 1
            else:
                tr_exon_blocks_list.append({
                    'geneinfo_id': k,
                    'exon_id': exon_ids[i],
                    'transcribed': False,
                    # 'cds_exon_id': None,
                    'baselen': baselen,
                    'from_tss_start': s.start,
                    'from_tss_stop': s.stop,
                    'from_cds_start': None,
                    'from_cds_stop': None,
                    'genome_start': gstart,
                    'genome_stop': gstop,
                    'sequence': 'N'*baselen
                })
    exon_df = pd.DataFrame(tr_exon_blocks_list)
    # Set index to geneinfo_id + exon_id (0-based) + transcribed
    exon_df = exon_df.set_index(
        ['geneinfo_id', 'exon_id', 'transcribed'])
    exon_df = exon_df[['baselen',
                       'genome_start', 'genome_stop',
                       'from_tss_start', 'from_tss_stop',
                       'from_cds_start', 'from_cds_stop',
                       'sequence']]

    # intron blocks df
    for k, d in intron_d.items():
        gid, _ = k
        if tr_meta_d[gid]['forward']:
            gstart = tr_meta_d[gid]['genome_start'] + d['from_tss_start']
            gstop = tr_meta_d[gid]['genome_start'] + d['from_tss_stop']
        else:
            gstart = tr_meta_d[gid]['genome_start'] - d['from_tss_start']
            gstop = tr_meta_d[gid]['genome_start'] - d['from_tss_stop']
        intron_d[k]['int_id'] -= 1  # convert to 0-based indexing
        intron_d[k]['genome_start'] = gstart
        intron_d[k]['genome_stop'] = gstop
    intron_df = pd.DataFrame(list(intron_d.values()))
    # Set index to geneinfo_id + int_id (1-based)
    intron_df = intron_df.set_index(['geneinfo_id', 'int_id'])
    # Add cds postions
    # Get cds relative start position
    cds_start = exon_df.query('transcribed == True')['from_tss_start'] \
        .groupby('geneinfo_id').min()
    # Create a temp join df
    tmp_df = intron_df.join(cds_start, on='geneinfo_id', rsuffix='_cds')
    tmp_df['from_cds_start'] = np.nan
    tmp_df['from_cds_stop'] = np.nan
    tmp_df.loc[tmp_df['transcribed'] == True, 'from_cds_start'] = \
        tmp_df['from_tss_start'] - tmp_df['from_tss_start_cds']
    tmp_df.loc[tmp_df['transcribed'] == True, 'from_cds_stop'] = \
        tmp_df['from_tss_stop'] - tmp_df['from_tss_start_cds']
    # Put columns into intron_df
    intron_df['from_cds_start'] = tmp_df['from_cds_start'].copy()
    intron_df['from_cds_stop'] = tmp_df['from_cds_stop'].copy()
    intron_df = intron_df[['transcribed', 'baselen',
                           'genome_start', 'genome_stop',
                           'from_tss_start', 'from_tss_stop',
                           'from_cds_start', 'from_cds_stop',
                           'sequence']]
    return tr_df, exon_df, intron_df
