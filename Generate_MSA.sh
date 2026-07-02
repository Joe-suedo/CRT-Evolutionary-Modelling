#!/bin/bash
#SBATCH --job-name=msa
#SBATCH --output=msa.out
#SBATCH --error=msa.err
#SBATCH --nodes=1
#SBATCH --mem=80G

# 1. Create the directory first
mkdir ~/msa
cd msa
mkdir species

# 2. Load Anaconda
module load anaconda/2024.10

# 3. Initialize bash
source $(conda info --base)/etc/profile.d/conda.sh

# 4. Add conda channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

# 5. Create new environment
conda create -y -n MSA

# 6. Activate the environment
conda activate MSA

# 7. Install tools
conda install bioconda::entrez-direct -y  #Entrez Direct (EDirect) - Access to NCBI's Entrez databases.
conda install bioconda::clustalo -y #Latest version of Clustal: a multiple sequence alignment program for DNA or proteins
conda install bioconda::mview -y #MView extracts and reformats the results of a sequence database search or multiple alignment. 

# 8. Download the sequences
while read -r organism sequence;do
    echo "Downloading ${organism} (${sequence})";
    efetch -db protein -id "${sequence}" -format fasta >tmp; #Download data to a temporary folder
    sed "1s/>.*/>${organism}/" tmp > "species/${organism}.faa"; #Rename the header of each file
    sleep 5
done < ../species.txt
rm -rf tmp #remove temporary folder

# 9. Concatenate all sequences into one fasta file
cat ./species/*faa > crt_sequences.faa

# 10.Perform a Multiple Sequence Alignment
clustalo -i crt_sequences.faa --seqtype=Protein --full -o crt_sequences.aln --outfmt=fa --output-order=tree-order --force

# 11.Generate visual html
mview -in fasta -html data -moltype aa -css on -bold -coloring group -colormap P1 -pagecolor "#1e1e1e" -textcolor "#ffffff" -alncolor "#2d2d2d" crt_sequences.aln > crt_alignment.html
