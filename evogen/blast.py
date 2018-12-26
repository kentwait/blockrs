# -*- coding: utf-8 -*-
"""Functions to filter and analyze BLAST results.
"""
import pandas as pd
import numpy as np
import subprocess as proc
from evogen.reader import blast_table_to_df


MAKEBLASTDB_CMD_TEMPLATE = 'makeblastdb -in {ffn} -input_type fasta -dbtype ' \
                           'nucl -out {out}'
BLASTN_CMD_TEMPLATE = 'blastn -task {task} -query {query} -db {db} ' \
                      '-out {report} -outfmt "{outfmt}" ' \
                      '-max_target_seqs 1 -num_threads {threads}'
BLASTN_TBL_TEMPLATE = '6 std qlen slen qcovs sstrand'


def make_blast_db(path, output_path):
    """Runs makeblastdb command within Python.
    Assumes that makeblastdb is in $PATH and can be called anywhere.

    Parameters
    ----------
    path : str
    output_path : str

    Returns
    -------
    int
        Return code. Returns 0 if makeblastdb successfully finished.

    Raises
    ------
    proc.CalledProcessError
        If the return code from running makeblastdb is non-zero, raises a
        CalledProcessError.

    """
    cmd = MAKEBLASTDB_CMD_TEMPLATE.format(ffn=path, out=output_path)
    try:
        p = proc.run(cmd, check=True, shell=True)
    except proc.CalledProcessError as e:
        raise e
    return p.returncode


def blastn_match(fasta_path, db_path, report_path,
                 task='blastn', outfmt='6 std qlen slen qcovs sstrand',
                 threads=1):
    """Runs BLAST search within Python using the blastn command.
    Assumes that blastn is in $PATH and can be called anywhere.

    Parameters
    ----------
    fasta_path : str
        Path to the fasta file used as query.
    db_path : str
        Path to the database files created by makeblastdb.
    report_path : str
        Path to where to store blast results.
    task : str, optional
        blastn task or mode to sue. By default, uses "blastn".
        Other tasks include "blastn-short", "megablast", and "dc-megablast".
    outfmt : str, optional
        Output format string. By default, uses the table format (6) with
        additional columns "qlen", "slen", "qcovs", "sstrand".
    threads : int, optional
        Number of threads blastn will use. It is recommended to change this
        to the number of cores in the computer.

    Returns
    -------
    pandas.DataFrame
        Blast

    Raises
    ------
    proc.CalledProcessError
        If the return code from running blastn is non-zero, raises a
        CalledProcessError.

    """
    cmd = BLASTN_CMD_TEMPLATE.format(
        task=task,
        query=fasta_path,
        db=db_path,
        report=report_path,
        outfmt=outfmt,
        threads=threads
    )
    try:
        proc.run(cmd, check=True, shell=True)
    except proc.CalledProcessError as e:
        raise e
    return blast_table_to_df(report_path)


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


def reciprocal_blast_match_df(forward_blast_df, reverse_blast_df, join='inner',
                              col_labels=['pident', 'length', 'mismatch',
                                          'gapopen', 'bitscore', 'qlen',
                                          'slen']):
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
    reverse_blast_df = reverse_blast_df.reorder_levels(['saccver', 'qaccver'])

    # Concat df
    concat_df = pd.concat([forward_blast_df, reverse_blast_df],
                        axis=1, join=join, names=['qaccver', 'saccver'])
    concat_df.columns = pd.MultiIndex.from_product(
        [['forward', 'reverse'], col_labels],
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
