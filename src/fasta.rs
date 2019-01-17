use pyo3::prelude::*;
use pyo3::PyObjectProtocol;

use std::fs::File;
use std::io::prelude::*;

#[pyclass]
/// Seq(id, description, sequence)
/// 
/// Seq represents one sequence in a fasta file. 
/// It is composed of id, description, and sequence properties which correspond to the following fields in the fasta file: 
/// >{id} {description}
/// {sequence}
pub struct Seq {

    #[prop(get, set)]
    pub id: String,

    #[prop(get, set)]
    pub description: String,

    #[prop(get, set)]
    pub sequence: String,
}

#[pymethods]
impl Seq {
    #[new]
    /// Creates a new Block object from start and stop values.
    fn __new__(obj: &PyRawObject, id: &str, description: &str, sequence: &str) -> PyResult<()> {
        obj.init(|_| {
            Seq::new(id, description, sequence)
        })
    }
}

// Seq::new method for Rust
impl Seq {
    /// Creates a new Block object from start and stop values.
    pub fn new(id: &str, description: &str, sequence: &str) -> Seq {
        Seq {
            id: id.to_string(),
            description: description.to_string(),
            sequence: sequence.to_string(),
        }
    }
}

// Customizes __repr__ and __str__ of PyObjectProtocol trait
#[pyproto]
impl PyObjectProtocol for Seq {
        fn __repr__(&self) -> PyResult<String> {
        Ok(format!("Seq(id=\"{id}\", description=\"{desc}\", sequence=\"{seq}\")", 
                  id=self.id, 
                  desc=self.description,
                  seq=self.sequence))
    }

    fn __str__(&self) -> PyResult<String> {
        match self.description.len() {
            0 => Ok(format!(">{id}\n{seq}", 
                  id=self.id, 
                  seq=self.sequence)),
            _ => Ok(format!(">{id} {desc}\n{seq}", 
                  id=self.id, 
                  desc=self.description,
                  seq=self.sequence))
        }
    }
}


#[pyfunction]
// TODO: Change output to Some<Vec<Seq>> to account for error
/// read_fasta(path)
/// 
/// Reads a fasta file into a list of one or more Seq objects.
pub fn read_fasta(path: &str) -> Vec<Seq> {
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
    let mut id: &str = "";
    let mut desc: &str = "";
    let mut sequence: Vec<String> = Vec::new();
    for line in data.lines() {
        if line.starts_with(">") {
            if sequence.len() > 0 {
                // Create a Seq struct and push to seq_list
                let seq = sequence.concat();
                let s = Seq::new(id, desc, &seq);
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
            id = line_contents[0]
        } else {
            sequence.push(line.trim_end().to_string());
        }
    }
    if sequence.len() > 0 {
        // Create a Seq struct and push to seq_list
        let seq = sequence.concat();
        let s = Seq::new(id, desc, &seq);
        seq_list.push(s);

        // Clear contents
    }
    
    seq_list 
}

#[pyfunction]
/// write_fastq(sequences, path, linewidth)
/// 
/// Writes a fasta file given a list of Seq objects.
pub fn write_fastq(sequences: Vec<&Seq>, path: &str, linewidth: i32) -> i32 {
    let mut file = match File::create(path) {
        Ok(file) => file,
        Err(x) => panic!("couldn't create {}: {}", path, x),
    };
    let mut str_list: Vec<String> = Vec::new();
    for seq in sequences {
        {
            let mut line = format!(">{}", seq.id);
            
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
