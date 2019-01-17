use pyo3::prelude::*;

use std::str;
use std::io::Write;

#[pyclass]
/// Block(start, stop)
/// 
/// Block represents an interval of positions from a given start position
/// until the stop position, exclusive.
/// 
/// For example, Block(0, 3) represents the range [0, 1, 2].
/// If the start value is >= 0, the stop value must be greater than the start.
/// 
/// For negative values, Block(-1, -4) represents the range [-1, -2, -3].
/// If the start value is < 0, the stop value must also be < 0 and less than the start.
pub struct Block {

    #[prop(get, set)]
    start: i32,

    #[prop(get, set)]
    stop: i32,
}

#[pymethods]
impl Block {
    #[new]
    /// Creates a new Block object from start and stop values.
    fn __new__(obj: &PyRawObject, start: i32, stop: i32) -> PyResult<()> {
        obj.init(|_| {
            Block { start, stop }
        })
    }
}

impl Block {
    fn new(start: i32, stop: i32) -> Block {
        Block { start, stop }
    }
}

#[pyfunction]
/// array_to_blocks(range_list)
/// 
/// Converts an explicit list of positions into a list of blocks.
/// Returns a list of Block objects.
fn array_to_blocks(range_list: Vec<i32>) -> Vec<Block> {
    let mut block_list: Vec<Block> = Vec::new();
    let mut start: i32 = range_list[0];
    let mut prev: i32 = range_list[0];

    for current in range_list.iter().skip(1) {
        let current: i32 = *current;
        if (prev >= 0) & (current >= 0) {
            if prev + 1 != current {
                let block = Block::new(start, prev + 1);
                block_list.push(block);
                start = current;
            }
        } else if (prev >= 0) & (current < 0) {
            let block = Block::new(start, prev + 1);
            block_list.push(block);
            start = current;
        } else if (prev < 0) & (current < 0) {
            if prev - 1 != current {
                let block = Block::new(start, prev - 1);
                block_list.push(block);
                start = current;
            }
        } else if (prev < 0) & (current >= 0) {
            let block = Block::new(start, prev - 1);
            block_list.push(block);
            start = current;
        }
        prev = current;
    }

    if prev >= 0 {
        let block = Block::new(start, prev + 1);
        block_list.push(block);
    } else {
        let block = Block::new(start, prev - 1);
        block_list.push(block);
    }

    block_list
}

fn option_array_to_blocks(range_list: Vec<Option<i32>>) -> Vec<Block> {
    let mut block_list: Vec<Block> = Vec::new();
    let mut prev_item: Option<i32> = range_list[0];

    let mut start: i32 = -1e9 as i32;

    for (a, b) in range_list.iter().zip(range_list.iter().skip(1)) {
        // Check if prev has a value or not
        if let Some(u) = a {
            let prev = *u;
            if start == -1e9 as i32 {
                start = prev;
            }
            if let Some(v) = b {
                // prev and current have values
                let current = *v;
                if (prev >= 0) && (current >= 0) {
                    if prev + 1 != current {
                        let block = Block::new(start, prev + 1);
                        block_list.push(block);
                        start = current;
                    }
                } else if (prev >= 0) && (current < 0) {
                    let block = Block::new(start, prev + 1);
                    block_list.push(block);
                    start = current;
                } else if (prev < 0) && (current < 0) {
                    if prev - 1 != current {
                        let block = Block::new(start, prev - 1);
                        block_list.push(block);
                        start = current;
                    }
                } else if (prev < 0) && (current >= 0) {
                    let block = Block::new(start, prev - 1);
                    block_list.push(block);
                    start = current;
                }
            } else {
                // prev has value, current is None
                if prev >= 0 {
                    let block = Block::new(start, prev + 1);
                    block_list.push(block);
                    start = -1e9 as i32;
                } else {
                    let block = Block::new(start, prev - 1);
                    block_list.push(block);
                    start = -1e9 as i32;
                }
            }
        } else {
            if let Some(v) = b {
                // prev is None, current has a value
                start = *v;
            }
        }
        prev_item = *b;
    }

    if let Some(u) = prev_item {
        // prev has a value
        let prev = u;
        if prev >= 0 {
            let block = Block::new(start, prev + 1);
            block_list.push(block);
        } else {
            let block = Block::new(start, prev - 1);
            block_list.push(block);
        }
    }
    

    block_list
}

#[pyfunction]
/// blocks_to_array(block_list)
/// 
/// Converts a list of Block objects into an explicit listing of positions.
/// Returns a list of integers.
fn blocks_to_array(block_list: Vec<&Block>) -> Vec<i32> {
    let mut pos_array: Vec<i32> = Vec::new();
    for block in block_list {
        let Block { start, stop } = *block;
        if start <= stop {
            for i in start..stop {
                pos_array.push(i);
            }
        } else {
            let start = start + 1;
            let stop = stop + 1;
            for i in (stop..start).rev() {
                pos_array.push(i);
            }
        }
    }
    pos_array
}

#[pyfunction]
/// pairwise_to_blocks(ref_seq, other_seq, debug)
/// 
/// Applies the positions of ungapped sites in the reference sequence unto the target.
/// Returns a list of Block objects
fn pairwise_to_blocks(ref_seq: &str, other_seq: &str, debug: bool) -> Vec<Block> {
    // Check if sequence lengths are the same
    // TODO: Change into an assert
    if ref_seq.len() != other_seq.len() {
        println!("sequence lengths are not the same");
    }

    // Declare variables
    // let ref_seq: &str = ref_seq.as_str();
    // let other_seq: &str = other_seq.as_str();
    let mut block_list: Vec<Block> = Vec::new();
    let mut seq_cnt: i32 = 0;
    let mut gap_cnt: i32 = 0;
    let mut start: i32 = 0;
    let mut prev_ref: char = ref_seq.chars().next().unwrap();
    let mut prev_seq: char = other_seq.chars().next().unwrap();
    let gap_char: char = '-';

    // Handle first column
    if prev_ref != gap_char && prev_seq != gap_char {
        if debug == true {
            print!("{} ", seq_cnt);
        }
        seq_cnt += 1;
    } else if prev_ref != gap_char && prev_seq == gap_char {
        if debug == true {
            print!("{} ", seq_cnt);
        }
        seq_cnt += 1;
    } else if prev_ref == gap_char && prev_seq == gap_char {
        if debug == true {
            print!("{} ", seq_cnt);
        }
    } else if prev_ref == gap_char && prev_seq != gap_char {
        if debug == true {
            print!("{} ", seq_cnt);
        }
        gap_cnt += 1;
    }
    for (a, b) in ref_seq.chars().skip(1).zip(other_seq.chars().skip(1)) {
        if a != gap_char && b != gap_char {
        // Current site are both filled.
            if gap_cnt == 0 && {prev_ref != gap_char && prev_seq != gap_char} {
                if debug == true {
                    print!("| {} ", seq_cnt);
                }
            } else if gap_cnt == 0 && {prev_ref != gap_char && prev_seq == gap_char} {
                if debug == true {
                    print!("_ {} ", seq_cnt);
                }
                start = seq_cnt;
            } else if gap_cnt == 0 && {prev_ref == gap_char && prev_seq == gap_char} {
                if debug == true {
                    print!("= {} ", seq_cnt);
                }
                start = seq_cnt;
            } else if gap_cnt > 0 && {prev_ref == gap_char && prev_seq != gap_char} {
                if debug == true {
                    print!("- {} ", seq_cnt);
                }
                let block = Block::new(-1, -1 - gap_cnt);
                block_list.push(block);
                gap_cnt = 0;
                start = seq_cnt;
            }
            seq_cnt += 1;

        } else if a != gap_char && b == gap_char {
        // Current site is filled in the reference and
        // gapped in the other sequence.
            if gap_cnt == 0 && {prev_ref != gap_char && prev_seq != gap_char} {
                if debug == true {
                    print!("| {} ", seq_cnt);
                }
                let block = Block::new(start, seq_cnt);
                block_list.push(block);
            } else if gap_cnt == 0 && {prev_ref != gap_char && prev_seq == gap_char} {
                if debug == true {
                    print!("_ {} ", seq_cnt);
                }
            } else if gap_cnt == 0 && {prev_ref == gap_char && prev_seq == gap_char} {
                if debug == true {
                    print!("= {} ", seq_cnt);
                }
            } else if gap_cnt > 0 && {prev_ref == gap_char && prev_seq != gap_char} {
                if debug == true {
                    print!("- {} ", seq_cnt);
                }
                let block = Block::new(-1, -1 - gap_cnt);
                block_list.push(block);
                gap_cnt = 0;
            }
            seq_cnt += 1;

        } else if a == gap_char && b == gap_char {
        // Current site is gapped for both.
            if gap_cnt == 0 && {prev_ref != gap_char && prev_seq != gap_char} {
                if debug == true {
                    print!("| {} ", seq_cnt);
                }
                let block = Block::new(start, seq_cnt);
                block_list.push(block);
            } else if gap_cnt == 0 && {prev_ref != gap_char && prev_seq == gap_char} {
                if debug == true {
                    print!("_ {} ", seq_cnt);
                }
            } else if gap_cnt == 0 && {prev_ref == gap_char && prev_seq == gap_char} {
                if debug == true {
                    print!("= {} ", seq_cnt);
                }
            } else if gap_cnt > 0 && {prev_ref == gap_char && prev_seq != gap_char} {
                if debug == true {
                    print!("- {} ", seq_cnt);
                }
                let block = Block::new(-1, -1 - gap_cnt);
                block_list.push(block);
                gap_cnt = 0;
            }
        } else if a == gap_char && b != gap_char {
        // Current site is gapped in the reference and
        // filled in the other sequence.
            if gap_cnt == 0 && {prev_ref != gap_char && prev_seq != gap_char} {
                if debug == true {
                    print!("| {} ", seq_cnt);
                }
                let block = Block::new(start, seq_cnt);
                block_list.push(block);
                start = seq_cnt;
            } else if gap_cnt == 0 && {prev_ref != gap_char && prev_seq == gap_char} {
                if debug == true {
                    print!("_ {} ", seq_cnt);
                }
            } else if gap_cnt == 0 && {prev_ref == gap_char && prev_seq == gap_char} {
                if debug == true {
                    print!("= {} ", seq_cnt);
                }
            } else if gap_cnt > 0 && {prev_ref == gap_char && prev_seq != gap_char} {
                if debug == true {
                    print!("- {} ", seq_cnt);
                }
            }
            gap_cnt += 1;
        }
    
        // Assign current characters
        prev_ref = a;
        prev_seq = b;
    }

    if gap_cnt > 0 {
        if debug == true {
            print!("- {}", seq_cnt);
        }
        let block = Block::new(-1, -1 - gap_cnt);
        block_list.push(block);
    } else {
        if gap_cnt == 0 && {prev_ref != gap_char && prev_seq != gap_char} {
            if debug == true {
                print!("| {} ", seq_cnt);
            }
            let block = Block::new(start, seq_cnt);
            block_list.push(block);
        } else if gap_cnt == 0 && {prev_ref != gap_char && prev_seq == gap_char} {
            if debug == true {
                print!("_ {} ", seq_cnt);
            }
        } else if gap_cnt == 0 && {prev_ref == gap_char && prev_seq == gap_char} {
            if debug == true {
                print!("= {} ", seq_cnt);
            }
        }
    }
    // Add a newline and flush
    if debug == true {
        print !("\n");
        std::io::stdout().flush().expect("error flushing stdout line buffer")
    }

    // TODO: Convert block_list vector to array
    block_list
}

#[pyfunction]
/// remove_sites(seq, block_list, remove_pos_list, gap_char)
/// 
/// Removes parts of the sequence using a list of positions and updated associated block list.
fn remove_sites(seq: &str, block_list: Vec<&Block>, mut remove_pos_list:  Vec<usize>, gap_char: &str) -> (String, Vec<Block>) {
    let gap_char = gap_char.chars().next().unwrap();
    // Unrolled blocks into positional array
    let block_array: Vec<i32> = blocks_to_array(block_list);

    // Create
    // absolute position array (block)
    let mut abs_pos_array: Vec<Option<i32>> = Vec::new();

    // When a seq is filled (not gap), put the corresponding block array value
    let mut x = 0;
    let mut new_seq_array: Vec<char> = seq.chars().collect();
    for c in seq.chars() {
        if c != gap_char {
            abs_pos_array.push(Some(block_array[x]));
            x += 1
        } else {
            abs_pos_array.push(None);
        }
    }

    // Remove positions
    remove_pos_list.sort_unstable();
    remove_pos_list.reverse();
    for i in remove_pos_list {
        abs_pos_array.remove(i);
        new_seq_array.remove(i);
    }

    // Reconstruct string
    let new_seq: String = new_seq_array.iter().collect();
    
    // Reconstruct blocks
    let new_block_list = option_array_to_blocks(abs_pos_array);

    (new_seq, new_block_list)
}

#[pymodinit]
fn block(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_function!(array_to_blocks)).unwrap();
    m.add_function(wrap_function!(blocks_to_array)).unwrap();
    m.add_function(wrap_function!(pairwise_to_blocks)).unwrap();
    m.add_function(wrap_function!(remove_sites)).unwrap();

    // Add Block class
    m.add_class::<Block>()?;

    Ok(())
}