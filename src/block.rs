use pyo3::prelude::*;
use pyo3::{PyObjectProtocol, exceptions};

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
        if stop < start {
            return Err(exceptions::ValueError::py_err("stop cannot be less than start"))
        }
        obj.init(|_| {
            Block { start, stop }
        })
    }

    #[staticmethod]
    fn from_block_str(data_str: &str) -> PyResult<Vec<Block>> {
        from_block_str(data_str)
    }

    #[staticmethod]
    fn to_block_str(block_list: Vec<&Block>) -> PyResult<String> {
        to_block_str(block_list)
    }
}

// Block::new method for Rust
impl Block {
    pub fn new(start: i32, stop: i32) -> Block {
        Block { start, stop }
    }
    pub fn check_new(start: i32, stop: i32) -> Result<Block, &'static str> {
        if stop < start {
            return Err("stop cannot be less than start")
        }
        Ok(Block::new(start, stop))
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
/// from_block_str(data_str)
/// 
/// Converts a block string into a list of Blocks.
fn from_block_str(data_str: &str) -> PyResult<Vec<Block>> {
    // Declare variables
    let mut block_list: Vec<Block> = Vec::new();

    // // Check if '_' exists
    // // Split str at '_' and get last substr
    // let coords_str = match data_str.rfind('_') {
    //     Some(i) => {
    //         let (_, string) = data_str.split_at(i);
    //         string.trim_start_matches('_')
    //     },
    //     None => data_str
    // };
    let coords_str = data_str;

    // TODO: Use regexp to check if coords_str is ^(\d+\:\d+\;*)+$

    // Split substr by ';'
    // For each split
    for start_stop in coords_str.split(';').collect::<Vec<&str>>() {
        // split again by ':' and convert to int
        let sep_idx = match start_stop.rfind(':') {
            Some(x) => x,
            None => return Err(exceptions::ValueError::py_err("separator \":\" not found"))
        };
        let (start, stop) = start_stop.split_at(sep_idx);
        let start = match start.parse::<i32>() {
            Ok(i) => i,
            Err(error) => return Err(exceptions::ValueError::py_err(format!("cannot convert start value from &str to i32: {:?} ", error))),
        };
        let stop = match stop.trim_start_matches(':').parse::<i32>() {
            Ok(i) => i,
            Err(error) => return Err(exceptions::ValueError::py_err(format!("cannot convert end value from &str to i32: {:?} ", error))),
        };
        // Create block and push to blocklist
        let block = match Block::check_new(start, stop) {
            Ok(x) => x,
            Err(x) => return Err(exceptions::ValueError::py_err(x)),
        };
        block_list.push(block)
    }
    // Returns block_list
    Ok(block_list)
}

#[pyfunction]
/// to_block_str(block_list)
/// 
/// Converts a list of Blocks into a block string.
fn to_block_str(block_list: Vec<&Block>) -> PyResult<String> {
    let mut block_str_vec: Vec<String> = Vec::new();
    for block in block_list {
        let substr = format!("{}:{}", block.start, block.stop);
        block_str_vec.push(substr);
    }
    Ok(block_str_vec.join(";"))
}

#[pyfunction]
/// array_to_blocks(range_list)
/// 
/// Converts an explicit list of positions into a list of blocks.
/// Returns a list of Block objects.
pub fn array_to_blocks(range_list: Vec<i32>) -> PyResult<Vec<Block>> {
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
                let block = match Block::check_new(start, prev + 1) {
                    Ok(x) => x,
                    Err(x) => return Err(exceptions::ValueError::py_err(x)),
                };
                block_list.push(block);
                start = current;
            }
        } else if (prev >= 0) & (current < 0) {
            let block = match Block::check_new(start, prev + 1) {
                Ok(x) => x,
                Err(x) => return Err(exceptions::ValueError::py_err(x)),
            };
            block_list.push(block);
            start = current;
        } else if (prev < 0) & (current < 0) {
            if prev - 1 != current {
                let block = match Block::check_new(start, prev - 1) {
                    Ok(x) => x,
                    Err(x) => return Err(exceptions::ValueError::py_err(x)),
                };
                block_list.push(block);
                start = current;
            }
        } else if (prev < 0) & (current >= 0) {
            let block = match Block::check_new(start, prev - 1) {
                Ok(x) => x,
                Err(x) => return Err(exceptions::ValueError::py_err(x)),
            };
            block_list.push(block);
            start = current;
        }
        // prev gets the value of current
        prev = current;
    }

    // End of the array
    // push last remaining block
    if prev >= 0 {
        let block = match Block::check_new(start, prev + 1) {
            Ok(x) => x,
            Err(x) => return Err(exceptions::ValueError::py_err(x)),
        };
        block_list.push(block);
    } else {
        let block = match Block::check_new(start, prev - 1) {
            Ok(x) => x,
            Err(x) => return Err(exceptions::ValueError::py_err(x)),
        };
        block_list.push(block);
    }

    // Returns block_list
    Ok(block_list)
}

/// Converts an explicit list of Option into a list of blocks.
/// Returns a list of Block objects.
fn option_array_to_blocks(range_list: Vec<Option<i32>>) -> Result<Vec<Block>, &'static str> {
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
                        let block = match Block::check_new(start, prev + 1) {
                            Ok(x) => x,
                            Err(x) => return Err(x),
                        };
                        block_list.push(block);
                        start = current;
                    }
                } else if (prev >= 0) && (current < 0) {
                    let block = match Block::check_new(start, prev + 1) {
                        Ok(x) => x,
                        Err(x) => return Err(x),
                    };
                    block_list.push(block);
                    start = current;
                } else if (prev < 0) && (current < 0) {
                    if prev - 1 != current {
                        let block = match Block::check_new(start, prev - 1) {
                            Ok(x) => x,
                            Err(x) => return Err(x),
                        };
                        block_list.push(block);
                        start = current;
                    }
                } else if (prev < 0) && (current >= 0) {
                    let block = match Block::check_new(start, prev - 1) {
                        Ok(x) => x,
                        Err(x) => return Err(x),
                    };
                    block_list.push(block);
                    start = current;
                }
            } else {
                // prev (a) has value, but current (b) is None
                if prev >= 0 {
                    let block = match Block::check_new(start, prev + 1) {
                        Ok(x) => x,
                        Err(x) => return Err(x),
                    };
                    block_list.push(block);
                    start = -1e9 as i32;
                } else {
                    let block = match Block::check_new(start, prev - 1) {
                        Ok(x) => x,
                        Err(x) => return Err(x),
                    };
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
            let block = match Block::check_new(start, prev + 1) {
                Ok(x) => x,
                Err(x) => return Err(x),
            };
            block_list.push(block);
        } else {
            let block = match Block::check_new(start, prev - 1) {
                Ok(x) => x,
                Err(x) => return Err(x),
            };
            block_list.push(block);
        }
    }
    
    Ok(block_list)
}

#[pyfunction]
/// blocks_to_array(block_list)
/// 
/// Converts a list of Block objects into an explicit listing of positions.
/// Returns a list of integers.
pub fn blocks_to_array(block_list: Vec<&Block>) -> PyResult<Vec<i32>> {
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
    Ok(pos_array)
}

#[pyfunction]
/// pairwise_to_blocks(ref_seq, other_seq, debug)
/// 
/// Applies the positions of ungapped sites in the reference sequence unto the target.
/// Returns a list of Block objects
pub fn pairwise_to_blocks(ref_seq: &str, other_seq: &str, gap_char: &str, debug: bool) -> PyResult<Vec<Block>> {
    // Check if sequence lengths are the same
    // TODO: Change into an assert
    if ref_seq.len() != other_seq.len() {
        return Err(exceptions::ValueError::py_err("sequence lengths are not the same"))
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
                let block = match Block::check_new(-1, -1 - gap_cnt) {
                    Ok(x) => x,
                    Err(x) => return Err(exceptions::ValueError::py_err(x)),
                };
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
                let block = match Block::check_new(start, seq_cnt) {
                    Ok(x) => x,
                    Err(x) => return Err(exceptions::ValueError::py_err(x)),
                };
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
                let block = match Block::check_new(-1, -1 - gap_cnt) {
                    Ok(x) => x,
                    Err(x) => return Err(exceptions::ValueError::py_err(x)),
                };
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
                let block = match Block::check_new(start, seq_cnt) {
                    Ok(x) => x,
                    Err(x) => return Err(exceptions::ValueError::py_err(x)),
                };
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
                let block = match Block::check_new(-1, -1 - gap_cnt) {
                    Ok(x) => x,
                    Err(x) => return Err(exceptions::ValueError::py_err(x)),
                };
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
                let block = match Block::check_new(start, seq_cnt) {
                    Ok(x) => x,
                    Err(x) => return Err(exceptions::ValueError::py_err(x)),
                };
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
        let block = match Block::check_new(-1, -1 - gap_cnt) {
            Ok(x) => x,
            Err(x) => return Err(exceptions::ValueError::py_err(x)),
        };
        block_list.push(block);
    } else {
        if gap_cnt == 0 && {prev_ref != gap_char && prev_seq != gap_char} {
            if debug == true {
                print!("| {} ", seq_cnt);
            }
            let block = match Block::check_new(start, seq_cnt) {
                Ok(x) => x,
                Err(x) => return Err(exceptions::ValueError::py_err(x)),
            };
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
    Ok(block_list)
}

#[pyfunction]
/// remove_sites(seq, block_list, remove_pos_list, gap_char)
/// 
/// Removes parts of the sequence using a list of positions and updated associated block list.
// TODO: Deal with empty result. I think this is panicking due to calling to option_array_to_blocks
pub fn remove_sites(seq: &str, block_list: Vec<&Block>, mut remove_pos_list:  Vec<usize>, gap_char: &str) -> PyResult<(String, Vec<Block>)> {
    // Check if remove_pos_list is empty
    if remove_pos_list.len() == 0 {
        let mut same_block_list: Vec<Block> = Vec::new();
        for block in block_list {
            let block = block.clone();
            same_block_list.push(block);
        }
        return Ok((seq.to_string(), same_block_list))
    }
    // Declare variables
    let gap_char = gap_char.chars().next().unwrap();
    // Unrolls blocks into positional array
    let block_array: Vec<i32> = match blocks_to_array(block_list) {
        Ok(x) => x,
        Err(x) => return Err(exceptions::ValueError::py_err(x)),
    };
    // absolute position array (block)
    let mut abs_pos_array: Vec<Option<i32>> = Vec::new();

    // When a seq is filled (not gap), put the corresponding block array value
    let mut x = 0;
    let mut new_seq_array: Vec<char> = seq.chars().collect();
    // Check if new_seq_array and abs_pos_array are the same length
    if block_array.len() != new_seq_array.len() {
        return Err(exceptions::ValueError::py_err("sequence length and block range are not equal"))
    }
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
        // Check if i is less than array size
        if i >= abs_pos_array.len() {
            return Err(exceptions::IndexError::py_err("position out of range"))
        }
        abs_pos_array.remove(i);
        new_seq_array.remove(i);
    }
    // Check if result is empty
    if new_seq_array.len() == 0 {
        let empty_string = String::new();
        let empty_vec: Vec<Block> = Vec::new();
        return Ok((empty_string, empty_vec))
    }
    // Reconstruct string
    let new_seq: String = new_seq_array.iter().collect();
    // Reconstruct blocks
    let new_block_list = match option_array_to_blocks(abs_pos_array) {
        Ok(x) => x,
        Err(x) => return Err(exceptions::ValueError::py_err(x)),
    };
    // Return new_seq and new_block_list
    Ok((new_seq, new_block_list))
}

// Register python functions to PyO3
#[pymodinit]
fn block(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_function!(from_block_str))?;
    m.add_function(wrap_function!(to_block_str))?;
    m.add_function(wrap_function!(array_to_blocks))?;
    m.add_function(wrap_function!(blocks_to_array))?;
    m.add_function(wrap_function!(pairwise_to_blocks))?;
    m.add_function(wrap_function!(remove_sites))?;

    // Add Block class
    m.add_class::<Block>()?;

    Ok(())
}