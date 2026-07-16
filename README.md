# CRT-Evolutionary-Modelling
This repository focuses on the evolutionary modeling of the Chloroquine Resistance Transporter (CRT) gene across diverse organisms. By reconstructing phylogenetic trees, this project aims to elucidate the evolutionary trajectory of the CRT gene, which plays a critical role in mediating antimalarial drug resistance.

## 📁 Repository Structure
### Up Stream Processing on HPC 
1. **Generate_MSA.sh:** Shell script that automates directory setup, downloads sequence data, and executes the Multiple Sequence Alignment (MSA) pipeline.
2. **species.txt:** A filtered list of organisms containing core *Plasmodium* strains and diverse evolutionary protists alongside their Accession IDs.
3. **msa/:** Output directory containing all sequence downloads and processing outputs.
4. **species/:** Subdirectory containing the individual FASTA sequences downloaded for each target protist.
5. **crt_sequences.faa:** Concatenated unaligned protein sequences (FASTA amino acid format).
6. **crt_sequences.aln:** The final Multiple Sequence Alignment file in FASTA format.
7. **crt_alignment.html:** Visual, color-coded HTML report of the finalized alignment.

### 📊 2. Downstream Analysis & Modeling (R)
The downstream phylogenetic modeling and visualization were conducted in R:

1. **msa.Rproj:** R Project environment configuration file to ensure seamless path management and working directory reproducibility.
2. **msa.R:** Comprehensive R script containing the workflow for distance-based matrix computations, maximum parsimony optimization, and maximum likelihood tree construction.
3. **Visualization:** Output plots and figures used to interpret the evolutionary relationships.
   
## 🚀 Workflow
### 🖥️ 1. Upstream Processing (HPC)
The upstream pipeline was executed on the ACE High-Performance Computing (HPC) cluster:
1. **Acquisition:** Programmatic retrieval of individual protein FASTA sequences from NCBI based on the curated `species.txt` list.
2. **Concatenation:** Consolidation of individual sequence files into a single, comprehensive FASTA database (`crt_sequences.faa`).
3. **MSA Generation:** Alignment of the concatenated database to generate the Multiple Sequence Alignment (MSA) using Clustal Omega.

