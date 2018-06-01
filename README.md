# Scop3D v3

 * [Project Description](#project-description)
 * [Usage](#usage)
 * [System requirements](#system-requirements)
 * [Previous version](#previous-version)
 * [Project Support](#project-support)

----

## Project Description

Scop3D (Sequence Conservation On Protein 3D structure) allows the visualization of sequence variation of a protein on its structure.

Scop3D is composed of two parts. The first part focuses on the analysis of sequence conservation or SNP variants, while the second part involves structural annotation. There is also an intermediate blast step available.

**New in version 3:** 
* added analysis of sequence SNP variants based on Uniprot ID, or DNA sequence(s)
* refreshed command line interface

### Sequence analysis

**Sequence conservation analysis:**
The starting point for the sequence conservation analysis is a FASTA-file which contains all sequence variants of interest. It is also possible to download sequence variants using UniprotID - the sequence will be blasted and cutoff value applied. These variants are aligned with MUSCLE. From the resulting multiple sequence alignment, the consensus sequence is constructed by retaining the most frequently found residue across all variant sequences for each position. This residue should be found more than the conservation threshold which is set be the user. Two matrices are created. The first matrix annotates the consensus sequence with the absolute abundance of each amino acid per position. In the second matrix, these numbers are converted into relative abundances in percent. Finally, a sequence logo is created with WebLogo.

**Sequence SNP variants analysis:**
The SNP variants analysis is based either on Uniprot and Ensembl data or on DNA sequence(s). In the first case the variant data is downloaded from the Uniprot and Ensembl databases using UniprotID. In the second case user-supplied fasta file with multiple sequences is used or, if using single DNA sequence, sequence is BLASTed against NCBI database with identity cutoff to obtain multiple sequences to be analyzed. Sequences are then aligned and compared to spot SNPs, as well as translated and aligned to obtain consensus protein sequence. In both cases using gathered info a matrix is created, annotating the protein sequence with the frequency, variants and type of variant. Based on that, the sequence is coloured accordingly. 

### Structure annotation

In this part, a PDB-file is needed. This file can be obtained from the PDB or an own PDB-file can be used. If no protein structure is at hand, scop3D provides the option to retrieve a homologous structure by performing a BLAST search of the consensus sequence against the PDB. The outcome of the structure annotation is two PDB-files with altered B-values and a set of figures. The output depends wheater the sequence conservation or SNP variants analysis was performed as the first step.

For conservation analysis in the first PDB-file, the B-values are replaced by the percent variation (defined as one hundred minus the percentage conservation) at that position. In the second PDB-file, the B-values are replaced by the entropy factor.

For SNP variants analysis in the first PDB-file, the B-values are replaced by the SNP frequency at that position. In the second PDB-file, the B-values are replaced by the SNP type.

### Citation
 * [Vermeire and Vermaere et al.: Proteomics. 2015 Apr;15(8):1448-52.](http://www.ncbi.nlm.nih.gov/pubmed/25641949)
 * If you use scop3D as part of a paper, please include the reference above.

----

## Usage

Scop3D is available either as web tool or as standalone command line tool.

### Web interface

Fully functional web interface is available under [http://example.com](TBD)

### Command line

Couple of typical flows are shown below. All files are read and written to `~/outputdir`. You can run `./scop3d -h` to see all available options.

**Sequence conservation analysis**

```
Analyze the sequence entropy and conservation for protein sequences from sequences.fasta:
# ./scop3d sequence -ent -w ~/outputdir -f sequences.fasta
Amongst result files in the outputdir, a consensus.fasta is created, that can be used as a BLAST input:
# ./scop3d blast ~/outputdir/consensus.fasta
After you analyze the blast output, you can choose the PDB ID, on which you want the consensus sequence entropy and conservation superimposed, e.g. 4F15:
# ./scop3d structure -ent -w ~/outputdir -i 4F15
```

```
Analyze the sequence entropy and conservation for protein A8K2U0 from Uniprot with 50% identity cutoff:
# ./scop3d sequence -ent -w ~/outputdir -u A8K2U0 -t 50
Superimpose sequence entropy and conservation on PDB from file 4F15.pdb using only chains A and C:
# ./scop3d structure -ent -w ~/outputdir -f 4F15.pdb -c A,C
```

**Sequence SNP variants analysis**

```
Analyze the SNPs for protein P06858 from Uniprot:
# ./scop3d sequence -var -w ~/outputdir -u P06858
# ./scop3d blast ~/outputdir/uniprot.fasta
# ./scop3d structure -var -w ~/outputdir -i 4ACQ
```

```
Analyze the SNPs for DNA sequences from sequences.fasta:
# ./scop3d sequence -var -w ~/outputdir -f sequences.fasta
# ./scop3d blast ~/outputdir/uniprot.fasta
# ./scop3d structure -var -w ~/outputdir -i 4ACQ
```

```
Analyze the SNPs for DNA sequence from dnasequence.fasta with 70% identity cutoff:
# ./scop3d sequence -var -w ~/outputdir -s dnasequence.fasta -t 70
# ./scop3d structure -var -w ~/outputdir -f 4F15.pdb
```


### Remarks

Please use -w to bundle all output/input files in one directory. The subsequent runs will search for input files in this directory.

----

## System requirements

* Python >= 2.7
* Perl >= 5
* Biopython
* Bioperl
* Ensembl Perl API
* Weblogo
* Muscle
* Emboss

----

## Previous version

The documentation for previous version of the tool (including binaries for windows, linux, mac) can be found in  [https://github.com/vibbits/scop3d/tree/master/scop3D.v2](scop3D.v2) directory

----

## Project Support

The scop3D project is grateful for the support by:

* [VIB Bioinformatics Core](http://www.bits.vib.be)
* [Compomics](http://www.compomics.com)
* [VIB](http://www.vib.be)
* [Ghent University](http://www.ugent.be)

[Go to top of page](#scop3d-v3)
