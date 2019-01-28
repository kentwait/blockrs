use pyo3::prelude::*;
use pyo3::PyObjectProtocol;

use std::str;
use std::io::Write;

#[pyclass]
#[derive(Copy, Clone)]
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
    pub start: i32,

    #[prop(get, set)]
    pub stop: i32,
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

// Block::new method for Rust
impl Block {
    pub fn new(start: i32, stop: i32) -> Block {
        Block { start, stop }
    }
}

// Customizes __repr__ and __str__ of PyObjectProtocol trait
#[pyproto]
impl PyObjectProtocol for Block {
    fn __repr__(&self) -> PyResult<String> {
        Ok(format!("Block({start}, {stop})", start=self.start, stop=self.stop))
    }

    fn __str__(&self) -> PyResult<String> {
        Ok(format!("{start}:{stop}", start=self.start, stop=self.stop))
    }
}

#[pyfunction]
/// array_to_blocks(range_list)
/// 
/// Converts an explicit list of positions into a list of blocks.
/// Returns a list of Block objects.
pub fn array_to_blocks(range_list: Vec<i32>) -> Vec<Block> {
    // Declare variables
    let mut block_list: Vec<Block> = Vec::new();
    let mut start: i32 = range_list[0];
    let mut prev: i32 = range_list[0];  // prev is the first value

    // Iterate over a list of positions skipping the first one
    // current is &i32
    for current in range_list.iter().skip(1) {
        // Shadows current, assigns concrete value of current
        let current: i32 = *current;
        // Checks current value and current relative to prev
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
        // prev gets the value of current
        prev = current;
    }

    // End of the array
    // push last remaining block
    if prev >= 0 {
        let block = Block::new(start, prev + 1);
        block_list.push(block);
    } else {
        let block = Block::new(start, prev - 1);
        block_list.push(block);
    }

    // Returns block_list
    block_list
}

/// Converts an explicit list of Option into a list of blocks.
/// Returns a list of Block objects.
fn option_array_to_blocks(range_list: Vec<Option<i32>>) -> Vec<Block> {
    // Declare variables
    let mut block_list: Vec<Block> = Vec::new();
    let mut prev_item: Option<i32> = range_list[0];  // option instead of i32
    let mut start: i32 = -1e9 as i32;

    // Iterate over pairs of prev (a) and currrent (b) positions
    // such that if range_list is [0, 1, 2, 3, 5, 6, 7] the the zipped list is
    // [(0, 1), (1, 2), (2, 3), (3, 5), (5, 6), (6, 7), (7, None)]
    for (a, b) in range_list.iter().zip(range_list.iter().skip(1)) {
        // Check if prev (a) has a value or not by destructuring
        if let Some(u) = a {
            let prev = *u;  // assign to prev is value is present
            if start == -1e9 as i32 {
                start = prev;  // assign prev as start
            }
            // Check if prev (a) has a value or not by destructuring
            if let Some(v) = b {
                // Reaching this point means  prev and current have values
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
                // prev (a) has value, but current (b) is None
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
            // prev (a) has no value, check current (b)
            if let Some(v) = b {
                // assign to start if prev has no value but current has
                start = *v;
            }
        }
        // current becomes prev_item
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
pub fn blocks_to_array(block_list: Vec<&Block>) -> Vec<i32> {
    // Declare variables
    let mut pos_array: Vec<i32> = Vec::new();
    // Iterate over list of Block
    for block in block_list {
        // Object destructure to get start and stop values
        let Block { start, stop } = *block;
        // start less than stop means both are positive
        if start <= stop {
            for i in start..stop {
                pos_array.push(i);
            }
        // start is larger than stop, probably means -1, -N
        } else {
            let start = start + 1;
            let stop = stop + 1;
            // Iterate in reverse like range(..,..,-1) in Python
            // This will go from -1, -2, ... , -N for example
            for i in (stop..start).rev() {
                pos_array.push(i);
            }
        }
    }
    // Return position array
    pos_array
}

#[pyfunction]
/// pairwise_to_blocks(ref_seq, other_seq, debug)
/// 
/// Applies the positions of ungapped sites in the reference sequence unto the target.
/// Returns a list of Block objects
pub fn pairwise_to_blocks(ref_seq: &str, other_seq: &str, gap_char: &str, debug: bool) -> Vec<Block> {
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
    let gap_char = gap_char.chars().next().unwrap();

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
    // Zip ref_seq and other_seq chars, both skipping the first character
    // ref_seq char is a, other_seq char is b
    for (a, b) in ref_seq.chars().skip(1).zip(other_seq.chars().skip(1)) {
        // Current site are both filled.
        if a != gap_char && b != gap_char {
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

        // Current site is filled in the reference and
        // gapped in the other sequence.
        } else if a != gap_char && b == gap_char {
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

        // Current site is gapped for both.
        } else if a == gap_char && b == gap_char {
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
        
        // Current site is gapped in the reference and
        // filled in the other sequence.
        } else if a == gap_char && b != gap_char {
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
    // Pairwise alignment is finished
    // Handle last block
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
// TODO: Deal with empty result. I think this is panicking due to calling to option_array_to_blocks
pub fn remove_sites(seq: &str, block_list: Vec<&Block>, mut remove_pos_list:  Vec<usize>, gap_char: &str) -> (String, Vec<Block>) {
    // Check if remove_pos_list is empty
    if remove_pos_list.len() == 0 {
        let mut same_block_list: Vec<Block> = Vec::new();
        for block in block_list {
            let block = block.clone();
            same_block_list.push(block);
        }
        return (seq.to_string(), same_block_list)
    }
    // Declare variables
    let gap_char = gap_char.chars().next().unwrap();
    // Unrolls blocks into positional array
    let block_array: Vec<i32> = blocks_to_array(block_list);
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
    // Check if result is empty
    if new_seq_array.len() == 0 {
        let empty_string = String::new();
        let empty_vec: Vec<Block> = Vec::new();
        return (empty_string, empty_vec)
    }
    // Reconstruct string
    let new_seq: String = new_seq_array.iter().collect();
    // Reconstruct blocks
    let new_block_list = option_array_to_blocks(abs_pos_array);
    // Return new_seq and new_block_list
    (new_seq, new_block_list)
}

// Register python functions to PyO3
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