use pyo3::prelude::*;

use std::fs::File;
use std::io::prelude::*;

#[pyclass]
/// Seq(seq_id, description, sequence)
struct Seq {
    seq_id: String,
    description: String,
    sequence: String,
}

#[pymethods]
impl Seq {
    #[new]
    /// Creates a new Block object from start and stop values.
    fn __new__(obj: &PyRawObject, seq_id: &str, description: &str, sequence: &str) -> PyResult<()> {
        obj.init(|_| {
            Seq::new(seq_id, description, sequence)
        })
    }
}

impl Seq {
    fn new(seq_id: &str, description: &str, sequence: &str) -> Seq {
        Seq {
            seq_id: seq_id.to_string(),
            description: description.to_string(),
            sequence: sequence.to_string(),
        }
    }
}

#[pyfunction]
// TODO: Change output to Some<Vec<Seq>> to account for error
fn read_fasta(path: &str) -> Vec<Seq> {
    // Open path in read-only mode
    // Returns io::Result<File>
    let mut file = match File::open(path) {
        Ok(file) => file,
        Err(x) => panic!("couldn't read {}: {}", path, x),
    };

    // Load all contents into data
    let mut data = String::new();
    match file.read_to_string(&mut data) {
        Ok(_) => (),
        Err(x) => panic!("couldn't read {}: {}", path, x),       
    };
    
    let mut seq_list: Vec<Seq> = Vec::new();
    let mut seq_id: &str = "";
    let mut desc: &str = "";
    let mut sequence: Vec<String> = Vec::new();
    for line in data.lines() {
        if line.starts_with(">") {
            if sequence.len() > 0 {
                // Create a Seq struct and push to seq_list
                let seq = sequence.concat();
                let s = Seq::new(seq_id, desc, &seq);
                seq_list.push(s);

                // Clear contents
                sequence.clear();
            }
            // Process the ID line
            // Separate the ID field from the description field
            // Check whether a description field is present
            let line_contents: Vec<&str> = line.trim_end()
                                               .trim_start_matches('>')
                                               .splitn(2, ' ')
                                               .collect();
            if line_contents.len() == 2 {
                desc = line_contents[1];
            }
            seq_id = line_contents[0]
        } else {
            sequence.push(line.trim_end().to_string());
        }
    }
    if sequence.len() > 0 {
        // Create a Seq struct and push to seq_list
        let seq = sequence.concat();
        let s = Seq::new(seq_id, desc, &seq);
        seq_list.push(s);

        // Clear contents
    }
    
    seq_list 
}

// // Wraps read_fasta in order to be exportable
// fn py_read_fasta(_py: Python, path: &str) -> PyResult<Vec<Seq>> {
//     let out = read_fasta(path);
//     Ok(out)
// }

#[pyfunction]
fn write_fastq(sequences: Vec<&Seq>, path: &str, linewidth: i32) -> i32 {
    let mut file = match File::create(path) {
        Ok(file) => file,
        Err(x) => panic!("couldn't create {}: {}", path, x),
    };
    let mut str_list: Vec<String> = Vec::new();
    for seq in sequences {
        {
            let mut line = format!(">{}", seq.seq_id);
            
            if seq.description.len() > 0 {
                line = format!("{} {}", line, seq.description);
            }
            str_list.push(line)
        }

        if linewidth == -1 {
            str_list.push(seq.sequence.to_string()); 
        } else if linewidth > 0 {
            let sub_seqs = seq.sequence.as_bytes()
                .chunks(linewidth as usize)
                .map(|s| unsafe { ::std::str::from_utf8_unchecked(s).to_string() })
                .collect::<Vec<_>>();
            str_list.extend(sub_seqs);
        } else {
            panic!("line width must be > 0 or -1")
        }
    }
    let data = str_list.join("\n");
    match file.write_all(data.as_bytes()) {
        Ok(_) => (),
        Err(x) => panic!("couldn't write to {}: {}", path, x)
    }

    1
}

// // Wraps write_fastq in order to be exportable
// fn py_write_fastq(_py: Python, sequences: Vec<PyDict>, path: &str, linewidth: i32) -> PyResult<i32> {
//     let mut sequences_: Vec<Seq> = Vec::new();
//     for dict in sequences {
//         let s = Seq {
//             seq_id: match dict.get_item(_py, "seq_id") {
//                 Some(v) => format!("{}", v),
//                 None => panic!("seq_id key not found!"),
//             },
//             description: match dict.get_item(_py, "description") {
//                 Some(v) => format!("{}", v),
//                 None => panic!("description key not found!"),
//             },
//             sequence: match dict.get_item(_py, "sequence") {
//                 Some(v) => format!("{}", v),
//                 None => panic!("sequence key not found!"),
//             },
//         };
//         sequences_.push(s);
//     }
//     let out = write_fastq(sequences_, path, linewidth);
//     Ok(out)
// }

#[pymodinit]
fn fasta(_py: Python, m: &PyModule) -> PyResult<()> {
    m.add_function(wrap_function!(read_fasta)).unwrap();
    m.add_function(wrap_function!(write_fastq)).unwrap();

    // Add Block class
    m.add_class::<Seq>()?;

    Ok(())
}

// #[cfg(test)]
// mod tests {
//     #[test]
//     fn it_works() {
//         assert_eq!(2 + 2, 4);
//     }
// }
