#[macro_use]
extern crate cpython;

use cpython::{Python, PyResult};

use std::str;
use std::io::Write;

fn array_to_blocks(range_list: Vec<i32>) -> Vec<(i32, i32, i32)> {
    let block_list: Vec<(i32, i32, i32)> = Vec::new();
    let mut start = block_list[0];
    let mut prev = block_list[0];

    for current in block_list.iter().skip(1) {
        if prev >= 0 & current >= 0 {
            if prev + 1 != current {
                block_list.push((start, prev + 1, 1));
                start = current;
            }
        } else if prev >= 0 & current < 0 {
            block_list.push((start, prev + 1, 1));
            start = current;
        } else if prev < 0 & current < 0 {
            if prev - 1 != current {
                block_list.push((start, prev - 1, -1));
                start = current;
            }
        } else if prev < 0 & current >= 0 {
            block_list.push((start, prev - 1, -1));
            start = current;
        }
        prev = current;
    }

    if prev >= 0 {
        block_list.push((start, prev + 1, 1));
    } else {
        block_list.push((start, prev - 1, -1));
    }

    block_list
}

fn pairwise_to_blocks(ref_seq: &str, other_seq: &str, debug: bool) -> Vec<(i32, i32, i32)> {
    // Check if sequence lengths are the same
    // TODO: Change into an assert
    if ref_seq.len() != other_seq.len() {
        println!("sequence lengths are not the same");
    }

    // Declare variables
    // let ref_seq: &str = ref_seq.as_str();
    // let other_seq: &str = other_seq.as_str();
    let mut block_list: Vec<(i32, i32, i32)> = Vec::new();
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
                let block = (-1, -1 - gap_cnt, -1);
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
                let block = (start, seq_cnt, 0);
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
                let block = (-1, -1 - gap_cnt, -1);
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
                let block = (start, seq_cnt, 0);
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
                let block = (-1, -1 - gap_cnt, -1);
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
                let block = (start, seq_cnt, 0);
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
        let block = (-1, -1 - gap_cnt, -1);
        block_list.push(block);
    } else {
        if gap_cnt == 0 && {prev_ref != gap_char && prev_seq != gap_char} {
            if debug == true {
                print!("| {} ", seq_cnt);
            }
            let block = (start, seq_cnt, 0);
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

// Wraps array_to_blocks in order to be exportable
fn py_array_to_blocks(_py: Python, range_list: Vec<i32>) -> PyResult<Vec<(i32, i32, i32)>> {
    let out = array_to_blocks(range_list);
    Ok(out)
}

// Wraps pairwise_to_blocks in order to be exportable
fn py_pairwise_to_blocks(_py: Python, ref_seq: &str, other_seq: &str, debug: bool) -> PyResult<Vec<(i32, i32, i32)>> {
    let out = pairwise_to_blocks(ref_seq, other_seq, debug);
    Ok(out)
}

py_module_initializer!(blockcodec, initblockcodec, PyInit_blockcodec, |py, m| { 
    m.add(py, "array_to_blocks", py_fn!(py, 
        py_array_to_blocks(range_list: Vec<i32>)))?;
    m.add(py, "pairwise_to_blocks", py_fn!(py, 
        py_pairwise_to_blocks(ref_seq: &str, other_seq: &str, debug: bool)))?;

    Ok(())
});

#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
