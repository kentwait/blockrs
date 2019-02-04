# -*- coding: utf-8 -*-
"""Functions to create and modify block data.
"""
import numpy as np
from blockrs import libblock

def array_to_blocks(range_list, debug=False):
    """Converts an explicit list of positions into a list of blocks.

    Parameters
    ----------
    range_list : list of int

    Returns
    -------
    list of slice
        Each block is represented as a slice (0-indexed)

    Example
    -------
    >>> lst = [1, 3, 5, 6, 7]
    >>> array_to_blocks(lst)
    [slice(1, 2, None), slice(3, 4, None), slice(5, 8, None)]
    >>> lst = [0, 1, -1, -2, 5, 6, 7]
    >>> array_to_blocks(lst)
    [slice(0, 2, None), slice(-1, -3, None), slice(5, 8, None)]

    """
    return libblock.array_to_blocks(range_list)


def blocks_to_array(block_list, debug=False):
    """Converts a list of blocks (slices) into an explicit
    listing of positions.

    Parameters
    ----------
    block_list : list of slice

    Returns
    -------
    numpy.array

    """
    pos_array = libblock.blocks_to_array(block_list)
    return np.array(pos_array, dtype=np.int32)


def to_zero_index_slice(s):  # pylint: disable=invalid-name
    """Convert 1-based to 0-based position.

    Parameters
    ----------
    s : slice

    Returns
    -------
    slice

    Example
    -------
    Assumes that values are both positive or both negative.
    Coordinate values are positive:
    >>> s1 = slice(1, 5)
    >>> to_zero_index_slice(s1)
    slice(0, 5, None)
    >>> s2 = slice(10, 50)
    >>> to_zero_index_slice(s2)
    slice(9, 50, None)

    Coordinate values are negative:
    >>> s3 = slice(-1, -5)
    >>> to_zero_index_slice(s3)
    slice(-1, -6, None)
    >>> s4 = slice(-5, -10)
    >>> to_zero_index_slice(s4)
    slice(-5, -11, None)

    """
    if s.start < 0:
        return slice(s.start, s.stop-1, s.step)
    return slice(s.start-1, s.stop, s.step)


def to_one_index_slice(s):  # pylint: disable=invalid-name
    """Convert 0-based to 1-based position.

    Parameters
    ----------
    s : slice

    Returns
    -------
    slice

    Example
    -------
    Assumes that values are both positive or both negative.
    Coordinate values are positive:
    >>> s1 = slice(0, 5)
    >>> to_one_index_slice(s1)
    slice(1, 5, None)
    >>> s2 = slice(10, 50)
    >>> to_one_index_slice(s2)
    slice(11, 50, None)

    Coordinate values are negative:
    >>> s3 = slice(-1, -5)
    >>> to_one_index_slice(s3)
    slice(-1, -4, None)
    >>> s4 = slice(-5, -10)
    >>> to_one_index_slice(s4)
    slice(-5, -9, None)

    """
    if s.start < 0:
        return slice(s.start, s.stop+1, s.step)
    return slice(s.start+1, s.stop, s.step)


def pairwise_to_blocks(ref_seq, other_seq, gap_char='-', debug=False):
    """Creates a list of blocks as slices relativeto a reference sequence.

    Parameters
    ----------
    ref_seq : str
    other_seq : str
    gapchar : str, optional

    Returns
    -------
    list

    Notes
    -----
    There are three possible scenarios when encoding the position of a sequence
    relative to a reference.

    Scenario 1:
    The reference sequence in ungapped while the other sequence has
    gaps.

                012345678
        ref   : ATGGAGCTG
        other : ATG---CAG
                012   678

        output : [slice(0,3), slice(6,9)]

    In this case, counting simply skips the positions where the
    gaps occur and resumes when a non-gap character is encountered.

    Since the output contains two slice blocks and the end
    position of the first does not match the start position of the second,
    this indicates that a gap exists in the sequence relative to the reference.
    The size of the gap can be calculated by the difference between the start
    of the latter slice and the end of the former.

    Scenario 2:
    The reference sequence and the other sequence contains gaps in the
    same positions.

                012   345
        ref   : ATG---CTG
        other : ATG---CAG
                012   345

        output : [slice(0,3), slice(3, 6)]

    This is also a simple case to deal with as gaps can be ignored when
    counting. This pairwise alignment is similar to the following when
    the gaps are removed:

                012345
        ref   : ATGCTG
        other : ATGCAG
                012345

        output : [slice(0, 6)]

    While the treatment is the same, the output is different between the two.
    To indicate that a gap existed, the slice block is split into two. Unlike
    the first scenario, the end of the first block has the same value as the
    start of the second block.

    Scenario 3:
    The reference sequence is gapped but the other sequence is not.

                012   345
        ref   : ATG---CTG
        other : ATGGAGCAG
                012???345

        output : [slice(0,3), slice(-4, -1), slice(3,6)]

    Since counting is relative to the reference, the positions where non-gap
    characters exist in the other sequence have an "undefined" position.
    To solve this problem, a GapBlock is used to indicate "nothingness".

    Note the similarity to scenario two where the end of the first slice block
    and the start of the second slice block have the same value. The presence
    of the "negative" slice block of length 3 between the two positive
    slice blocks indicates that there are three positions in between
    that have an undefined position in the reference.

    """
    assert len(ref_seq) == len(other_seq), 'sequence lengths are not the same'
    return libblock.pairwise_to_blocks(ref_seq, other_seq, gap_char, debug)


def remove_sites(seq, block_list, removed_pos_list, gap_char='-', debug=False):
    """Removes parts of the sequence using a list of positions and
    updated associated block list.

    Parameters
    ----------
    seq : str
    block_list : list of slice
        Each block is represented as a slice (0-indexed).
        The combined range of all blocks must be equal to the
        length of seq.
    removed_pos_list : list of int
        Positions to be removed from the sequence.
        0 refers to the first character in seq and
        len(seq) - 1 refers to the last character.
    zero_indexed : bool
        If True, block list positions use Python 0-based indexing, such that
        slice(0,3) -> [0,1,2], or slice(-1, -4) -> [-1,-2,-3]
        Otherwise, uses 1-based indexing, such that slice(1,3) -> [1,2,3],
        or slice(-1, -3) -> [-1,-2,-3]

    Returns
    -------
    str, list of slice
        Returns a tuple of the sequence and the block list.

    """
    return libblock.remove_sites(seq, block_list, removed_pos_list, gap_char)


# Combine exon block lists to get UTR exon blocks
def combine_exon_blocks(tr_exon_blocks, cds_exon_blocks):
    # Assumes 0-indexed slices

    # Create transcript exon d
    raw_exon_blocks = sorted(
        [slice(s.start, s.stop, -1) for s in tr_exon_blocks] +
        [slice(s.start, s.stop, 1) for s in cds_exon_blocks]
    )
    all_exon_blocks = []
    exon_id_list = []
    i = 0
    id_cnt = 0
    while i < len(raw_exon_blocks):
        try:
            # Case 1: (0, 10, cds=False) and (5, 10, cds=True)
            # Split to (0, 5, cds=False), (5, 10, cds=True)
            if (raw_exon_blocks[i].start < raw_exon_blocks[i+1].start) and \
               (raw_exon_blocks[i].stop == raw_exon_blocks[i+1].stop):
                # Split
                first_block = slice(raw_exon_blocks[i].start,
                                    raw_exon_blocks[i+1].start,
                                    raw_exon_blocks[i].step)
                all_exon_blocks.append(first_block)
                all_exon_blocks.append(raw_exon_blocks[i+1])
                exon_id_list += [id_cnt, id_cnt]
                i += 2
            # Case 2: (0, 5, cds=True), (0, 10, cds=False)
            # Split to (0, 5, cds=True), (5, 10, cds=False)
            elif (raw_exon_blocks[i].start == raw_exon_blocks[i+1].start) and \
                 (raw_exon_blocks[i].stop < raw_exon_blocks[i+1].stop):
                # Split
                second_block = slice(raw_exon_blocks[i].stop,
                                     raw_exon_blocks[i+1].stop,
                                     raw_exon_blocks[i+1].step)
                all_exon_blocks.append(raw_exon_blocks[i])
                all_exon_blocks.append(second_block)
                exon_id_list += [id_cnt, id_cnt]
                i += 2
            # Case 3: (0, 5, cds=True), (0, 5, cds=False)
            # Adopt the cds block
            elif (raw_exon_blocks[i].start == raw_exon_blocks[i+1].start) and \
                 (raw_exon_blocks[i].stop == raw_exon_blocks[i+1].stop):
                # Overwrite
                cds_i = i if raw_exon_blocks[i].step == 1 else i+1
                cds_block = slice(raw_exon_blocks[cds_i].start,
                                  raw_exon_blocks[cds_i].stop,
                                  raw_exon_blocks[cds_i].step)
                all_exon_blocks.append(cds_block)
                exon_id_list += [id_cnt]
                i += 2
            # Default case: No overlap with next block
            # Add current block
            else:
                all_exon_blocks.append(raw_exon_blocks[i])
                exon_id_list.append(id_cnt)
                i += 1
            id_cnt += 1
        except IndexError:
            all_exon_blocks.append(raw_exon_blocks[i])
            exon_id_list.append(id_cnt)
            i += 1
            id_cnt += 1
    return all_exon_blocks, exon_id_list

def append_blocks_to_aln(geneinfo_id, fasta_d, ref_func, ignore_func, template):
    ref_seq = ''
    for seq_id, seq in fasta_d.items():
        if ref_func(seq_id):
            ref_seq = seq.sequence
        elif ignore_func(seq_id):
            continue
        
        block_list = pairwise_to_blocks(ref_seq, seq.sequence, False)
        block_str = ';'.join([str(s) for s in block_list])
        
        seq.description = template.format(geneinfo_id=geneinfo_id,
                                          num_blocks=len(block_list),
                                          block_str=block_str) + \
                          ((' ' + seq.description) if seq.description else '')

