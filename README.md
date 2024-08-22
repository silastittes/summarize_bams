# summarize_bams
Produce data about sites of interest from a list of bam files.

```bash
usage: summarize_bams.py [-h] -o OUTFILE chromosome_position_file bamlist

Summarize BAM files

positional arguments:
  chromosome_position_file
                        List of chromosome and position
  bamlist               List of bam files to summarize

options:
  -h, --help            show this help message and exit
  -o OUTFILE, --outfile OUTFILE
                        Output file
```

# Minimal working example
```bash
python summarize_bams.py -o summarize_bams_out.txt example/example_sites.txt example/example_bams.txt```
```
Note: bam files must be indexed.

# Output
```
ID       chromosome  position  balance  count  depth  base  base_qual  mapping_qual  duplicate  paired  proper  rc     secondary  supplementary  unmapped  qc_fail  read1  read2
example  chr1        36        1.0      1      2      T     27         20            False      True    True    False  False      False          False     False    True   False
example  chr1        36        1.0      1      2      T     27         40            False      True    True    False  False      False          False     False    True   False
```

Hopefully this output is self-explanatory. Goal is to output as much information as possible about the sites in the chromosome_position_file in a way that is easy to parse and analyze.

One happy little unit test too
```bash
pytest -q summarize_bams.py
```
