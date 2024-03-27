use crate::errors::WGAError;
use crate::parser::chain::{ChainHeader, ChainReader, ChainRecord};
use crate::parser::cigar::{
    parse_cigar_to_blocks, parse_cigar_to_chain, parse_cigar_to_insert, parse_maf_seq_to_chain,
};
use crate::parser::common::{AlignRecord, Strand};
use crate::parser::maf::{MAFReader, MAFRecord, MAFSLine, MAFWriter};
use crate::parser::paf::PAFReader;
use crate::utils::reverse_complement;
use noodles::sam::header::record::value::map;
use noodles::sam::header::record::value::map::header::SortOrder;
use noodles::sam::record::ReadName;
use noodles::sam::{
    self as sam,
    header::record::value::{
        map::{Program, ReferenceSequence},
        Map,
    },
};
use rayon::prelude::*;
use rust_htslib::faidx;
use std::io::{Read, Write, BufWriter};
use std::{fs::File, fs::read_dir, fs::read_to_string};
use std::num::NonZeroUsize;
use std::collections::{HashMap, HashSet};
use tempfile::tempdir;

/// Convert a MAF Reader to output a PAF file
pub fn maf2paf<R: Read + Send>(
    mafreader: &mut MAFReader<R>,
    writer: &mut dyn Write,
) -> Result<(), WGAError> {
    // init csv writer for deserializing
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer);

    // multi-threading
    let pafrecords = mafreader
        .records()
        .par_bridge()
        .map(|record| -> Result<_, WGAError> {
            let mafrecord = record?;
            mafrecord.convert2paf()
        })
        .collect::<Result<Vec<_>, WGAError>>()?;
    for pafrec in pafrecords {
        wtr.serialize(pafrec)?;
    }
    wtr.flush()?;
    Ok(())
}

/// Convert a MAF Reader to output a Chain file
pub fn maf2chain<R: Read + Send>(
    mafreader: &mut MAFReader<R>,
    writer: &mut Box<dyn Write>,
) -> Result<(), WGAError> {
    // iterate over records and give a self-increasing chain-id
    for (id, record) in mafreader.records().enumerate() {
        let record = record?;

        // transform record to Chain Header
        let mut header = ChainHeader::from(&record);

        // set chain id
        header.chain_id = id;

        // write header without newline
        writer.write_all(format!("{}", header).as_bytes())?;

        // nom the cigar string and write to file
        parse_maf_seq_to_chain(&record, writer)?;

        // additional newline for standard chain format
        writer.write_all(b"\n\n")?;
    }
    writer.flush()?;
    Ok(())
}

/// Convert a MAF Reader to output a MSA file
pub fn maf2msa<R: Read + Send>(
    mafreader: &mut MAFReader<R>,
    writer: &mut Box<dyn Write>,
) -> Result<(), WGAError> {

    // Temporary directory of aligned sequences
    let tmp_dir = tempdir()?;

    // Map of sequence ids to writers to the corresponding temp file
    let mut name_to_writer: HashMap<String, Box<BufWriter<dyn Write>>> = HashMap::new();
    
    // Keep track of names seen
    let mut all_names: HashSet<String> = HashSet::new();

    // Keep track of alignment size in case we come across a new sequence
    let mut total_cols: usize = 0;

    // iterate over records and append each sequence to a temporary
    // fasta file according to the seq name
    for record in mafreader.records() {
        let record = record?;
        let mut block_cols: usize = 0;
        let mut block_names: HashSet<String> = HashSet::new();

        // Add each alignment entry to the corresponding temp file
        for block_entry in record.slines {
            block_cols = block_entry.seq.len();

            // Only include first copy of multi-copy paralogs
            // This may introduce some bias, but the alternative is to drop the entire block.
            // Dropping the whol block would reduce noise at the cost of potentially
            // tossing out useful info. 
            if block_names.contains(&block_entry.name) {
                log::warn!(
                    "{} is present multiple times in a single block, only using first occurence.", 
                    &block_entry.name);
                continue;
            }
            block_names.insert(block_entry.name.to_string());

            if !name_to_writer.contains_key(&block_entry.name) {
                let file_path = tmp_dir
                    .path()
                    .join(format!("seq-{}.fa", name_to_writer.len()));
                let file = File::create(file_path)?;
                name_to_writer
                    .insert(block_entry.name.to_string(), Box::new(BufWriter::new(file)));
                name_to_writer
                    .get_mut(&block_entry.name)
                    .unwrap()
                    .write(format!(">{}\n", block_entry.name).as_bytes())?;

                // If we've already written another block which didn't include this one,
                // we need to pad it so that the columns match up
                if total_cols > 0 {
                    name_to_writer
                        .get_mut(&block_entry.name)
                        .unwrap()
                        .write("-".repeat(total_cols).as_bytes())?;
                }
            } 

            // Write the sequence
            name_to_writer
                .get_mut(&block_entry.name)
                .unwrap()
                .write(block_entry.seq.as_bytes())?;
        }

        // Add gaps to missing sequence names
        for missing_name in all_names.difference(&block_names) {
            name_to_writer
                .get_mut(missing_name)
                .unwrap()
                .write("-".repeat(block_cols).as_bytes())?;
        }

        total_cols += block_cols;
        all_names.extend(block_names);
    }

    // Add newline and flush all buffers
    name_to_writer.iter_mut()
        .for_each(|(&ref _name, bw)| {bw.write("\n".as_bytes()).unwrap(); bw.flush().unwrap();});
    drop(name_to_writer);

    // Write all files to final MSA file
    for path in read_dir(&tmp_dir).unwrap() {
        writer.write_all(read_to_string(path.unwrap().path())?.as_bytes())?;
    }

    // additional newline for standard chain format
    writer.flush()?;
    Ok(())
}

pub fn maf2sam<R: Read + Send>(
    _mafreader: &mut MAFReader<R>,
    writer: &mut Box<dyn Write>,
) -> Result<(), WGAError> {
    let mut sam_writer = sam::Writer::new(writer);
    let mut header = Map::<map::Header>::default();
    *header.sort_order_mut() = Some(SortOrder::Unsorted);
    let header = sam::Header::builder()
        .set_header(header)
        .add_reference_sequence(
            "sq0".parse()?,
            Map::<ReferenceSequence>::new(NonZeroUsize::try_from(8)?),
        )
        .add_reference_sequence(
            "sq1".parse()?,
            Map::<ReferenceSequence>::new(NonZeroUsize::try_from(13)?),
        )
        .add_reference_sequence(
            "sq2".parse()?,
            Map::<ReferenceSequence>::new(NonZeroUsize::try_from(21)?),
        )
        .add_program("noodles-sam", Map::<Program>::default())
        .add_comment("an example SAM written by noodles-sam")
        .build();

    sam_writer.write_header(&header)?;
    let mut record = sam::alignment::Record::default();
    let read_name: ReadName = "sq2".parse()?;
    *record.read_name_mut() = Some(read_name);
    sam_writer.write_record(&header, &record)?;
    Ok(())
}

/// Convert a PAF Reader to output a Blocks file
pub fn paf2blocks<R: Read + Send>(
    pafreader: &mut PAFReader<R>,
    writer: &mut dyn Write,
) -> Result<(), WGAError> {
    // init writer and csv writer for deserializing
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);

    // iterate over records
    for record in pafreader.records() {
        let record = record?;
        // nom the cigar string and write to file
        parse_cigar_to_blocks(&record, &mut wtr)?;
    }
    wtr.flush()?;
    Ok(())
}

/// Convert a PAF Reader to output a Chain file
pub fn paf2chain<R: Read + Send>(
    pafreader: &mut PAFReader<R>,
    writer: &mut Box<dyn Write>,
) -> Result<(), WGAError> {
    // iterate over records and give a self-increasing chain-id
    for (id, record) in pafreader.records().enumerate() {
        let record = record?;

        // transform record to Chain Header
        let mut header = ChainHeader::from(&record);

        // set chain id
        header.chain_id = id;

        // write header without newline
        writer.write_all(format!("{}", header).as_bytes())?;

        // nom the cigar string and write to file
        parse_cigar_to_chain(&record, writer)?;

        // additional newline for standard chain format
        writer.write_all(b"\n\n")?;
    }
    writer.flush()?;
    Ok(())
}

/// Convert a PAF Reader to output a MAF file
pub fn paf2maf<R: Read + Send>(
    pafreader: &mut PAFReader<R>,
    writer: &mut dyn Write,
    t_fa_path: &str,
    q_fa_path: &str,
) -> Result<(), WGAError> {
    // get the target and query fasta reader
    let t_reader = faidx::Reader::from_path(t_fa_path)?;
    let q_reader = faidx::Reader::from_path(q_fa_path)?;

    // init a MAFWriter
    let mut mafwtr = MAFWriter::new(writer);

    // write header
    let header = format!(
        "#maf version=1.6 convert_from=paf t_seq_path={} q_seq_path={}",
        t_fa_path, q_fa_path
    );
    mafwtr.write_header(header)?;

    for pafrec in pafreader.records() {
        let pafrec = pafrec?;
        // get mapq as score
        let score = pafrec.mapq;
        // get target info
        let t_name = &pafrec.target_name;
        let t_start = pafrec.target_start;
        let t_end = pafrec.target_end - 1;
        let t_strand = pafrec.target_strand();
        let t_alilen = pafrec.target_end - pafrec.target_start;
        let t_size = pafrec.target_length;
        // get query info
        let q_name = &pafrec.query_name;
        let q_strand = pafrec.query_strand();
        let q_size = pafrec.query_length;
        let q_alilen = pafrec.query_end - pafrec.query_start;
        // NOTE: if negative strand, we should convert the start position
        let q_start = match q_strand {
            Strand::Positive => pafrec.query_start,
            Strand::Negative => q_size - pafrec.query_end,
        };

        // get seqs from indexed fasta files
        let mut whole_t_seq =
            t_reader.fetch_seq_string(t_name, t_start as usize, t_end as usize)?;
        let mut whole_q_seq = q_reader.fetch_seq_string(
            q_name,
            pafrec.query_start as usize,
            (pafrec.query_end - 1) as usize,
        )?;

        // reverse complement the query sequence if it is on the negative strand
        match q_strand {
            Strand::Positive => {}
            Strand::Negative => {
                whole_q_seq = reverse_complement(&whole_q_seq)?;
            }
        }
        // nom the cigar string and insert the `-` to sequence
        parse_cigar_to_insert(&pafrec, &mut whole_t_seq, &mut whole_q_seq)?;
        // get s-lines
        let t_sline = MAFSLine {
            mode: 's',
            name: t_name.to_string(),
            start: t_start,
            align_size: t_alilen,
            strand: t_strand,
            size: t_size,
            seq: whole_t_seq,
        };
        let q_sline = MAFSLine {
            mode: 's',
            name: q_name.to_string(),
            start: q_start,
            align_size: q_alilen,
            strand: q_strand,
            size: q_size,
            seq: whole_q_seq,
        };
        // get maf record
        let mafrec = MAFRecord {
            score,
            slines: vec![t_sline, q_sline],
        };
        // write maf record
        mafwtr.write_record(&mafrec)?;
    }
    Ok(())
}

/// Convert a Chain Reader to output a MAF file
pub fn chain2maf<R: Read + Send>(
    chainreader: &mut ChainReader<R>,
    writer: &mut dyn Write,
    t_fa_path: &str,
    q_fa_path: &str,
) -> Result<(), WGAError> {
    // get the target and query fasta reader
    let t_reader = faidx::Reader::from_path(t_fa_path)?;
    let q_reader = faidx::Reader::from_path(q_fa_path)?;

    // init a MAFWriter
    let mut mafwtr = MAFWriter::new(writer);

    // write header
    let header = format!(
        "#maf version=1.6 convert_from=chain t_seq_path={} q_seq_path={}",
        t_fa_path, q_fa_path
    );
    mafwtr.write_header(header)?;

    for chainrec in chainreader.records()? {
        let chainrec = chainrec?;
        // 255 as score
        let score = 255;
        // get target info
        let t_name = chainrec.target_name();
        let t_start = chainrec.target_start();
        let t_end = chainrec.target_end() - 1;
        let t_strand = chainrec.target_strand();
        let t_alilen = chainrec.target_end() - chainrec.target_start();
        let t_size = chainrec.target_length();
        // get query info
        let q_name = chainrec.query_name();
        let q_strand = chainrec.query_strand();
        let q_size = chainrec.query_length();
        let q_alilen = chainrec.query_end() - chainrec.query_start();
        // NOTE: if negative strand, we should convert the start position
        let q_start = match q_strand {
            Strand::Positive => chainrec.query_start(),
            Strand::Negative => q_size - chainrec.query_end(),
        };

        // get seqs from indexed fasta files
        let mut whole_t_seq =
            t_reader.fetch_seq_string(t_name, t_start as usize, t_end as usize)?;
        let mut whole_q_seq = q_reader.fetch_seq_string(
            q_name,
            chainrec.query_start() as usize,
            (chainrec.query_end() - 1) as usize,
        )?;

        // reverse complement the query sequence if it is on the negative strand
        match q_strand {
            Strand::Positive => {}
            Strand::Negative => {
                whole_q_seq = reverse_complement(&whole_q_seq)?;
            }
        }
        // read chain dataline and insert the `-` to sequence
        parse_chain_to_insert(&chainrec, &mut whole_t_seq, &mut whole_q_seq)?;
        // get s-lines
        let t_sline = MAFSLine {
            mode: 's',
            name: t_name.to_string(),
            start: t_start,
            align_size: t_alilen,
            strand: t_strand,
            size: t_size,
            seq: whole_t_seq,
        };
        let q_sline = MAFSLine {
            mode: 's',
            name: q_name.to_string(),
            start: q_start,
            align_size: q_alilen,
            strand: q_strand,
            size: q_size,
            seq: whole_q_seq,
        };
        // get maf record
        let mafrec = MAFRecord {
            score,
            slines: vec![t_sline, q_sline],
        };
        // write maf record
        mafwtr.write_record(&mafrec)?;
    }
    Ok(())
}

/// Parse the Chain Data Lines to insert the `-` to sequence
fn parse_chain_to_insert(
    rec: &ChainRecord,
    t_seq: &mut String,
    q_seq: &mut String,
) -> Result<(), WGAError> {
    let mut current_offset = 0;
    for dataline in &rec.lines {
        let ins_len = dataline.target_diff;
        let del_len = dataline.query_diff;
        current_offset += dataline.size;
        match ins_len {
            0 => {}
            _ => {
                let ins_str = "-".repeat(ins_len as usize);
                t_seq.insert_str(current_offset as usize, &ins_str);
                current_offset += ins_len;
            }
        }
        match del_len {
            0 => {}
            _ => {
                let del_str = "-".repeat(del_len as usize);
                q_seq.insert_str(current_offset as usize, &del_str);
                current_offset += del_len;
            }
        }
    }
    Ok(())
}

/// Convert a Chain Reader to output a PAF file
pub fn chain2paf<R: Read + Send>(
    chainreader: &mut ChainReader<R>,
    writer: &mut dyn Write,
) -> Result<(), WGAError> {
    // init csv writer for deserializing
    let mut wtr = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .has_headers(false)
        .from_writer(writer);

    // multi-threading
    let pafrecords = chainreader
        .records()?
        .par_bridge()
        .map(|record| -> Result<_, WGAError> {
            let chainrecord = record?;
            chainrecord.convert2paf()
        })
        .collect::<Result<Vec<_>, WGAError>>()?;
    // if we should sort pafrecords?
    for pafrec in pafrecords {
        wtr.serialize(pafrec)?;
    }
    wtr.flush()?;
    Ok(())
}
