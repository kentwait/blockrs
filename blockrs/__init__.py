# -*- coding: utf-8 -*-
from libblockrs import block as libblock
from libblockrs.block import Block, CatBlock
from libblockrs.block import to_block_str, remove_sites_from_blocks
from blockrs.block import array_to_blocks, blocks_to_array, pairwise_to_blocks
from blockrs.block import remove_sites

__author__ = 'Kent Kawashima'
__version__ = '0.6.1'
__all__ = ['libblock', 'Block', 'CatBlock',
           'to_block_str', 'remove_sites_from_blocks',
           'array_to_blocks', 'blocks_to_array',
           'pairwise_to_blocks', 'remove_sites']
