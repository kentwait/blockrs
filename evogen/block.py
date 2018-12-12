# -*- coding: utf-8 -*-
"""Functions to create and modify block data.
"""
import numpy as np

def range_to_blocks(range_list):
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
    [slice(1, 2, None), slice(3, 4, None), slice(5, 8, None)]
    >>> lst = [0, 1, -1, -2, 5, 6, 7]
    [slice(0, 2, None), slice(-1, -3, None), slice(5, 8, None)]

    """
    block_list = []
    start = range_list[0]
    prev = range_list[0]
    for current in range_list[1:]:
        if prev >= 0 and current >= 0:
            if prev + 1 != current:
                block_list.append(slice(start, prev+1))
                start = current
        elif prev >= 0 and current < 0:
            block_list.append(slice(start, prev+1))
            start = current
        elif prev < 0 and current < 0:
            if prev - 1 != current:
                block_list.append(slice(start, prev-1))
                start = current
        elif prev < 0 and current >= 0:
            block_list.append(slice(start, prev-1))
            start = current
        prev = current
    block_list.append(slice(start, prev+1 if prev >= 0 else prev-1))
    return block_list


def to_zero_index_slice(s):
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


def to_one_index_slice(s):
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


def remove_sites(seq, block_list, removed_pos_list, zero_indexed=True):
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
    # Unroll block_list into an positional array
    # [slice(0:5), slice(8:10)] -> [0, 1, 2, 3, 4, 8, 9]
    range_value = lambda x, y: (x, y, 1) if x >= 0 and y >= 0 else (x, y, -1)
    if not zero_indexed:
        block_list = list(map(to_zero_index_slice, block_list))
    block_pos_array = np.array(
        [i for block_slice in block_list
         for i in range(*range_value(block_slice.start, block_slice.stop))],
        dtype=np.int32)

    # Associate seq position to positions to be removed
    rel_pos_array = np.arange(0, len(seq), dtype=np.int32)
    if len(block_pos_array) != len(rel_pos_array):
        raise Exception('total length of blocks ({}) is not equal to the '
                        'length of the sequence ({})'.format(
                            len(block_pos_array), len(rel_pos_array)))

    # Remove positions
    rel_pos_array = np.delete(rel_pos_array, removed_pos_list)

    # Edit sequence
    seq_array = np.array(list(seq))[rel_pos_array]

    # Edit block
    block_pos_array = block_pos_array[rel_pos_array]
    new_block_list = range_to_blocks(block_pos_array)
    if not zero_indexed:
        new_block_list = list(map(to_one_index_slice, new_block_list))

    return ''.join(seq_array), new_block_list
