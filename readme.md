# XXX pipeline

## Description

## Requirement

How to install the required tools below?? `conda? docker?`

`Python version >= 3.7`

`bedtools`

`pandas`

`snakemake`

`fastqc`

`trimgalore`

`fastx-multx`

`samtools`

`bwa-mem`

`bowtie2`

## Input

## Output

general description of the pipeline done.

----

# SOP

----

## Build reference genome

### Object
Build a reference genome only having exon 4~7 from ABL1 gene

### How to build the reference genome
1. get hg19.fa and correspond gtf - refGene_05_14.2020.txt

   - `md5sum refGene_05_14.2020.txt` is `7d87863aa84702725c44af5302708f47`
   - `md5sum hg19.fa` is `d6851f9f4537ff4e9beb5b7a08b89230`
2. `python3 src/extract_exons.py -t NM_005157.4 --exon_nums 4 5 6 7`
3. `python3 src/squish_fa.py -i output/NM_005157.fa -o NM_005157_exons_4_to_7`

   - We will obtain `NM_005157_exons_4_to_7.fa` in `output` folder
4. Build indexes for your alignment tools `bwa or bowtie2`

### Requirement

`python3` >= 3.7
`bedtools`
`pandas`
`bwa`/`bowtie2`

### Python scripts Usage
Help

`python3 extract_exons.py -h`

example

`python3 extract_exons.py -t NM_xxxx --exon_nums 4 5 6 7`

### Output
- `NM_xxxx.fa`
- `NM_xxxx.bed`
- `NM_xxxx.fa.bak`
 
`NM_xxxx.fa` has sequences of the exons chosen by users

After running `extract_exons.py` I also ran `squish.py` to merge the `header/sequences` from the fa file.

#### Note
`NM_xxxx.bed` file has exon info which looks like the table below.

When you look at those loci in IGV, +1 should be added to the exon_start.


| chr | exon_start | exon_end  | exon_num |
|:----|:-----------|:----------|:---------|
| chr9| 133738149  | 133738422 | 4        |
| chr9| 133747515  | 133747600 | 5        |
| chr9| 133748246  | 133748424 | 6        |
| chr9| 133750254  | 133750439 | 7        |
