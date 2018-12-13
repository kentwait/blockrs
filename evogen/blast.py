# -*- coding: utf-8 -*-
"""Functions to filter and analyze BLAST results.
"""
import pandas as pd
import numpy as np


def summarize_group(group):
    """Flattens a group of rows from a blast DataFrame.

    Paramters
    ---------
    group
        Grouped rows in a DataFrame.
        Assumes that columns "qaccver", "saccver", "pident",
        "length", "mismatch", "gapopen", "bitscore", "qlen", "slen"
        are present.

    Returns
    -------
    pandas.Series

    """
    if np.max(group['pident']) == 100 and \
       (np.max(group['length']) == np.max(group['slen'])):
        return pd.Series({
            'qaccver': np.max(group['qaccver']),
            'saccver': np.max(group['saccver']),
            'pident':  np.max(group['pident']),
            'length': np.max(group['length']),
            'mismatch': np.min(group['mismatch']),
            'gapopen': np.min(group['gapopen']),
            'bitscore': np.max(group['bitscore']),
            'qlen': np.max(group['qlen']),
            'slen': np.max(group['slen']),
        })
    return pd.Series({
        'qaccver': np.max(group['qaccver']),
        'saccver': np.max(group['saccver']),
        'pident':  np.sum(group['length'] - group['mismatch'] -
                          group['gapopen']) /
                   np.max(group['slen']) * 100,
        'length': group['length'].sum(),
        'mismatch': group['mismatch'].sum(),
        'gapopen': group['gapopen'].sum(),
        'bitscore': group['bitscore'].sum(),
        'qlen': np.max(group['qlen']),
        'slen': np.max(group['slen']),
    })


def summarize_blast_df(df, summarize_func=summarize_group, 
                       eval_threshold=1e-10, pident_threshold=100.0,
                       summary_col_labels=['qaccver', 'saccver', 'pident',
                                           'length', 'mismatch', 'gapopen',
                                           'bitscore', 'qlen', 'slen'],
                       eval_col_name='evalue',
                       pident_col_name='pident'):
    """Collates the matches from the same match pair into one and
    summarizes selected columns of the blast result.

    Parameters
    ----------
    df : pandas.DataFrame
        DataFrame containing blast matches.
        Assumes that query and subject ID columns are named
        "qaccver" and "saccver" respectively.
    summarize_func : function, optional
        Function to summarize grouped match results.
        By default uses evogen.blast.summarize_group function.
    eval_threshold : float, optional
        Only results less than this value will be including in
        the summary. Other results will be ignored.
        Value must be greater than 0.
    pident_threshold : float, optional
        Only summarized results greater than or equal to this value
        will be returned.
        Value must be greater than or equal to zero and less than or
        equal to 100.
    summary_col_labels : list of str, optional
        Columns to be retained in the summarized DataFrame.
    eval_col_name : str, optional
        Name of the "e-value" column.
    pident_col_name : str, optional
        Name of the "percent identity" column.

    Returns
    -------
    pandas.DataFrame

    """
    # Check values
    if eval_threshold <= 0:
        raise ValueError('eval_threshold ({}) must be greater than '
                         'zero'.format(eval_threshold))
    if pident_threshold > 100:
        raise ValueError('pident_threshold ({}) must be less than or '
                         'equal to 100.0'.format(pident_threshold))
    elif pident_threshold < 0:
        raise ValueError('pident_threshold ({}) must be greater than '
                         'zero'.format(pident_threshold))

    # Filter raw results by e-value
    filtered_df = df[df[eval_col_name] < eval_threshold]
    if len(filtered_df) < 1:
        return filtered_df

    # Group and summarize
    grp_df = filtered_df.groupby(['qaccver', 'saccver'])[summary_col_labels]
    summ_df = grp_df.apply(summarize_func)[summary_col_labels[2:]]

    # Filter summarized results by pident
    filtered_summ_df = summ_df[summ_df[pident_col_name] >= pident_threshold]

    return filtered_summ_df


def reciprocal_blast_match_df(forward_blast_df, reverse_blast_df, join='inner'):
    """Returns a concatenated DataFrame containg the forward and reverse
    blast matches joined by an inner join.

    Parameters
    ----------
    forward_blast_df : pandas.DataFrame
    reverse_blast_df : pandas.DataFrame
    join : str, optional
        Type of join to use. By default, an inner join is performed.

    Returns
    -------
    pandas.DataFrame
        MultiIndex follows the query and subject IDs in forward_blast_df.

    """
    forward_blast_df.columns = ['f_' + label
                                for label in forward_blast_df.columns]
    reverse_blast_df.columns = ['r_' + label
                                for label in reverse_blast_df.columns]

    # Concat df
    concat_df = pd.concat([forward_blast_df, reverse_blast_df],
                          axis=1, join=join)
    concat_df.columns = pd.MultiIndex.from_product(
        [['forward', 'reverse'],
         ['pident', 'length', 'mismatch', 'gapopen', 'bitscore', 'qlen', 'slen']
        ],
        names=['match', 'property']
    )
    return concat_df


def reciprocal_blast_match_pairs(forward_blast_df, reverse_blast_df):
    """Returns the pair of ids that reciprocally match based on
    forward and reverse blast results.

    Parameters
    ----------
    forward_blast_df : pandas.DataFrame
    reverse_blast_df : pandas.DataFrame

    Returns
    -------
    set of tuple
        Sorted by value of the qeury ID in the forward blast matching.

    """
    concat_df = reciprocal_blast_match_df(forward_blast_df, reverse_blast_df,
                                          join='inner')
    return sorted(concat_df.index.get_values())
