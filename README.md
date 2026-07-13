# CRT-Evolutionary-Modelling
This repository focuses on the evolutionary modeling of the Chloroquine Resistance Transporter (CRT) gene across diverse organisms. By reconstructing phylogenetic trees, this project aims to elucidate the evolutionary trajectory of the CRT gene, which plays a critical role in mediating antimalarial drug resistance.

## 📁 Repository Structure
### Up Stream Processing on HPC 
* **Generate_MSA.sh:** Shell script that automates directory setup, downloads sequence data, and executes the Multiple Sequence Alignment (MSA) pipeline.
* **species.txt:** A filtered list of organisms containing core *Plasmodium* strains and diverse evolutionary protists alongside their Accession IDs.
* **msa/:** Output directory containing all sequence downloads and processing outputs.
* **species/:** Subdirectory containing the individual FASTA sequences downloaded for each target protist.
* **crt_sequences.faa:** Concatenated unaligned protein sequences (FASTA amino acid format).
* **crt_sequences.aln:** The final Multiple Sequence Alignment file in FASTA format.
* **crt_alignment.html:** Visual, color-coded HTML report of the finalized alignment.
