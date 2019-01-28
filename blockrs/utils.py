# -*- coding: utf-8 -*-
"""Helper functions and utilities
"""

def parse_shortform_block_annotation(description):
    """Parses shortform version of the block annotation from a string.

    Parameters
    ----------
    description : str

    Returns
    -------
    dict
        Returns a dictionary with keys 'ref_version', 'chromosome',
        'chromosome_scaffold', 'genome_pos', 'orientation', 'blocks'.

    Example
    -------
    >>> description = 'Dm528_2L_3_+_8366009_1_1:1773'
    >>> parse_shortform_block_annotation(description)
    {'ref_version': 'Dm528',
     'chromosome': '2L',
     'chromosome_scaffold': 3,
     'genome_pos': 8366008,
     'orientation': '+',
     'blocks': [slice(0, 1773, None)]}

    """
    if not description:
        return {}
    # TODO: Test encoded string formatting using regexp

    ref_version, chrom, scaffold, ori, gpos, _, coords = description.split('_')
    coords = [slice(int(start) - 1, int(stop))
              for start, stop in [c.split(':') for c in coords.split(';')]]

    return {
        'ref_version': ref_version,
        'chromosome': chrom,
        'chromosome_scaffold': int(scaffold),
        'genome_pos': int(gpos) - 1,
        'orientation': ori,
        'blocks': coords,
    }

def encode_shortform_block_annotation(d):
    """Creates the shortform version of the block annotation
    using information from the sequence dictionary.

    Parameters
    ----------
    d : dict
        Dictionary representation individuals entries in a FASTA file.

    Returns
    -------
    str

    """
    # Sample: 'Dm528_2L_3_+_8366009_1_1:1773'
    template = '{ref_version}_{chr_loc}_{chr_scaffold}_{ori}_{gpos_oneb}_' \
               '{num_blocks}_{blocks}'
    try:
        return template.format(
            ref_version=d['ref_version'],
            chr_loc=d['chromosome'],
            chr_scaffold=d['chromosome_scaffold'],
            ori=d['orientation'],
            gpos_oneb=d['genome_pos'] + 1,
            num_blocks=len(d['blocks']),
            blocks=';'.join(
                ['{}:{}'.format(s[0].start + 1, s[0].stop) for s in d['blocks']]
            ),
        )
    except KeyError:
        # If any one of the keys is not found, returns an empty string
        # Because "marker" will not have any extra keys in its per sequence dict
        return ''


def summarize_ancestral_prob_df(df):
    """Combines multiple entries that have identical ancestral species MRCA
    and population MRCA states into one using pandas groupby.

    Parameters
    ----------
    df : pandas.DataFrame

    Returns
    -------
    pandas.DataFrame

    """
    df = df.groupby(['pattern', 'allele_count_a', 'allele_count_b',
                     'anc_species_state', 'anc_pop_state',
                     'anc_species_pop']) \
           .apply(lambda x: x['joint_prob'].sum()) \
           .reset_index() \
           .set_index(['pattern', 'allele_count_a', 'allele_count_b'])
    df.columns = ['anc_species_state', 'anc_pop_state',
                  'anc_species_pop', 'prob']
    return df
