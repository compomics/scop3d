# Scop3D

 * [Project Description](#project-description)
 * [Downloads](#downloads)
 * [Usage](#usage)
 * [Project Support](#project-support)

----

## Project Description

Scop3D (Sequence Conservation On Protein 3D structure) allows the visualization of sequence variation of a protein on its structure.

Scop3D is composed of two parts. The first part focuses on the analysis of sequence conservation, while the second part involves structural annotation.

### Sequence analysis

The starting point for the sequence analysis is a FASTA-file which contains all sequence variants of interest. These variants are aligned with MUSCLE. From the resulting multiple sequence alignment, the consensus sequence is constructed by retaining the most frequently found residue across all variant sequences for each position. This residue should be found more than the conservation threshold which is set be the user. Two matrices are created. The first matrix annotates the consensus sequence with the absolute abundance of each amino acid per position. In the second matrix, these numbers are converted into relative abundances in percent. Finally, a sequence logo is created with WebLogo.

### Structure annotation

In this part, a PDB-file is needed. This file can be obtained from the PDB or an own PDB-file can be used. If no protein structure is at hand, scop3D provides the option to retrieve a homologous structure by performing a BLAST search of the consensus sequence against the PDB. The outcome of the structure annotation is two PDB-files with altered B-values and a set of figures. In the first PDB-file, the B-values are replaced by the percent variation (defined as one hundred minus the percentage conservation) at that position. In the second PDB-file, the B-values are replaced by the entropy factor.

### Citation
 * [Vermeire and Vermaere et al.: Proteomics. 2015 Apr;15(8):1448-52.](http://www.ncbi.nlm.nih.gov/pubmed/25641949)
 * If you use scop3D as part of a paper, please include the reference above.

[Go to top of page](#scop3d)

----

## Downloads

  * [Scop3D for windows](http://genesis.ugent.be/colims/scop3D-Windows-Version)
  * [Scop3D for mac](http://genesis.ugent.be/downloadredirect.php?toolname=scop3D-mac)
  * [Scop3D for linux](http://genesis.ugent.be/downloadredirect.php?toolname=scop3D-linux)

[Go to top of page](#scop3d)

----

## Usage
See the [wiki](https://github.com/compomics/thermo-msf-parser/wiki) for additional information on how to use *Thermo MSF Parser*.

[Go to top of page](#scop3d)

----

## Project Support

The scop3D project is grateful for the support by:

| Compomics | VIB | Ghent University|
|:--:|:--:|:--:|
| [![compomics](http://genesis.ugent.be/public_data/image/compomics.png)](http://www.compomics.com) | [![vib](http://genesis.ugent.be/public_data/image/vib.png)](http://www.vib.be) | [![ugent](http://genesis.ugent.be/public_data/image/ugent.png)](http://www.ugent.be/en) |

[Go to top of page](#scop3d)

----

| IntelliJ | Netbeans | Java | Maven |
|:--:|:--:|:--:|:--:|
| [![intellij](https://www.jetbrains.com/idea/docs/logo_intellij_idea.png)](https://www.jetbrains.com/idea/) | [![netbeans](https://netbeans.org/images_www/visual-guidelines/NB-logo-single.jpg)](https://netbeans.org/) | [![java](http://genesis.ugent.be/public_data/image/java.png)](http://java.com/en/) | [![maven](http://genesis.ugent.be/public_data/image/maven.png)](http://maven.apache.org/) |

[Go to top of page](#thermo-msf-parser)
