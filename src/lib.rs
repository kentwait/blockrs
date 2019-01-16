#![feature(specialization)]

#[macro_use]
extern crate pyo3;

use pyo3::prelude::*;

use std::str;
use std::io::Write;

#[pyclass]
struct Block {
    start: i32,
    stop: i32,
}

#[pymethods]
impl Block {
    #[new]
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

#[pyfunction]
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
fn blocks_to_array(block_list: Vec<Block>) -> Vec<i32> {
    let mut pos_array: Vec<i32> = Vec::new();
    for block in block_list {
        let Block { start, stop } = block;
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
fn remove_sites(seq: &str, block_list: Vec<Block>, mut remove_pos_list:  Vec<usize>, gap_char: &str) -> (String, Vec<Block>) {
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

// // Wraps array_to_blocks in order to be exportable
// fn py_array_to_blocks(_py: Python, range_list: Vec<i32>) -> PyResult<Vec<Block>> {
//     let out = array_to_blocks(range_list);
//     Ok(out)
// }

// // Wraps blocks_to_array in order to be exportable
// fn py_blocks_to_array(_py: Python, block_list: Vec<Block>) -> PyResult<Vec<i32>> {
//     let out = blocks_to_array(block_list);
//     Ok(out)
// }

// // Wraps pairwise_to_blocks in order to be exportable
// fn py_pairwise_to_blocks(_py: Python, ref_seq: &str, other_seq: &str, debug: bool) -> PyResult<Vec<Block>> {
//     let out = pairwise_to_blocks(ref_seq, other_seq, debug);
//     Ok(out)
// }

// // Wraps remove_sites in order to be exportable
// fn py_remove_sites(_py: Python, seq: &str, block_list: Vec<Block>, remove_pos_list: Vec<usize>, gap_char: &str) -> PyResult<(String, Vec<Block>)> {
//     let out = remove_sites(seq, block_list, remove_pos_list, gap_char);
//     Ok(out)
// }

// py_module_initializer!(blockcodec, initblockcodec, PyInit_blockcodec, |py, m| { 
//     m.add_function(py, "array_to_blocks", py_fn!(py, 
//         py_array_to_blocks(range_list: Vec<i32>)))?;
//     m.add_function(py, "blocks_to_array", py_fn!(py, 
//         py_blocks_to_array(block_list: Vec<Block>)))?;
//     m.add_function(py, "pairwise_to_blocks", py_fn!(py, 
//         py_pairwise_to_blocks(ref_seq: &str, other_seq: &str, debug: bool)))?;
//     m.add_function(py, "remove_sites", py_fn!(py, 
//         py_remove_sites(seq: &str, block_list: Vec<Block>, remove_pos_list: Vec<usize>, gap_char: &str)))?;
//     Ok(())
// });

#[pymodinit]
fn blockcodec(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_function!(array_to_blocks)).unwrap();
    m.add_function(wrap_function!(blocks_to_array)).unwrap();
    m.add_function(wrap_function!(pairwise_to_blocks)).unwrap();
    m.add_function(wrap_function!(remove_sites)).unwrap();

    Ok(())
}

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
