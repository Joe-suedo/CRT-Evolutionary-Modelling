#Tracing the Evolution of PfCRT Antimalarial Drug Resistance
#A Comparative Phylogenetic Analysis using Distance, Parsimony, & Likelihood

#Analyzing wild CRT sequences collected globally to determine the single most accurate evolutionary tree 
#showing how chloroquine resistance emerged.


# 1.INSTALLING & LOADING REQUIRED PACKAGES
install.packages(pkgs = "adegenet", dependencies = TRUE)
install.packages(pkgs = "phangorn", dependencies = TRUE)
library(stats)
library(ade4)
library(ape)
library(adegenet)
library(phangorn)

# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
#2. DISTANCE BASED PHLOGENETICS WITH APE(Naive baseline)
#.Computing pairwise distances between individuals, 
#Representing these distances on a tree and Evaluating relevance of this representation.

#a. Import global PfCRT sequence isolates
msa <- ape::read.FASTA(file = "crt_sequences.aln", type = "AA")
msa_matrix <- as.matrix(x=msa) #Create matrix for input in downstream steps

#b. Calculate raw percentage variations across different species
pair_dist <- ape::dist.aa(x=msa_matrix, pairwise.deletion = TRUE, scaled = TRUE)
matrix_pair_dist <- as.matrix(x=pair_dist) #Create a matrix of pairwise distances

#c. Paint a geographic/mutational divergence matrix to spot  clusters
df_matrix_pair_dist <- as.data.frame(matrix_pair_dist) #Data frame required
ade4::table.paint(df = df_matrix_pair_dist, csize = 1, clegend=1)

#d. Build a distance-based tree (NJ algorithm)
NJ_tree <- phangorn::NJ(pair_dist)
ape::plot.phylo(x = NJ_tree, plot = TRUE, main = "Distance tree Plot(NJ alogorithm)")

#e. Assessing quality of Phylogeny
#I.Compute cophenetic distance
cophenetic_dist <- stats::cophenetic(x =NJ_tree)
#II. Convert cophenetic distance into a distance matrix; to match the pair_dist object
cophenetic_dist <- stats::as.dist(cophenetic_dist) 
vector_cophenetic_dist <- as.vector(cophenetic_dist)
vector_pair_dist <- as.vector(pair_dist)
graphics::plot(x = vector_pair_dist, y= vector_cophenetic_dist, main = "QC of distance-based tree", 
               ylab = "Cophenetic(in the tree) Distance", xlab = "Pairwise(Original/True) distance",
               graphics::abline(a = 0, b = 1, col = "red", lty = 2))

# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*

#3. MAXIMUM PARSIMONY PHYLOGENIES(Simplest assumption)
#Seeks trees which minimize the total number of changes (substitutions) from ancestors to descendents.

#a. Reading data into R ; Create phydat 
msa2 <- phangorn::read.phyDat(file = "crt_sequences.aln", format = "fasta", type = "AA")

#c. Minimum number of mutation steps (substitutions) required to explain the NJ phylo tree in 1(c) above.
NJ_tree_score <- phangorn::parsimony(tree = NJ_tree, data = msa2)
cat("NJ_tree_score(Total Mutation Steps):", NJ_tree_score, "\n")

#d. Search for a better tree by modifying the branches to find the layout with the absolute lowest mutation score
lowest_mut_tree <- phangorn::optim.parsimony(tree =NJ_tree, data = msa2) 
cat("Tree with lowest mutation (Targeting lowest steps):", phangorn::parsimony(tree = lowest_mut_tree, data = msa2), "\n")

# e. Plot the Parsimony tree
ape::plot.phylo(x =lowest_mut_tree, plot = TRUE, main = "Phylogenetic Parsimony Tree(Lowest mutations)", use.edge.length = TRUE,
                align.tip.label = TRUE)

# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
#4. MAXIMUM LIKELIHOOD PHYLOGENIES
#. Browses a space of possible tree topologies looking for the ’best’ tree; 
#. more flexibility in that any model of sequence evolution can be taken into account
#. Also uses Phangorn and considered to be the statistically best option

#a. Initialize a biochemical evolutionary framework (Using WAG model + Gamma rates)
#.  pml computes the likelihood of a phylogenetic tree given a sequence alignment and a model. 
pml_object <- phangorn::pml(tree = NJ_tree, data = msa2, k = 4, model = "WAG") #WAG model for Amino acid data

#b. Optimize molecular parameters 
#.  optim.pml optimizes the different model parameters.
#   optNni(topology), optBf(base frequencies), optQ(rate matrix), optInv(proportion of variable size),
#   optGamma(gamma rate parameter), optEdge(edge lengths), optRate(overall rate), optRooted(edge lengths)
opt_pml_object <- phangorn::optim.pml(object = pml_object, optNni = TRUE, optBf = TRUE, optQ =  TRUE,
                           optGamma = TRUE, optEdge = TRUE)

# c. Statistical Verdict: Direct Model Comparison via Likelihood Ratio Test
# Anova:Compute analysis of variance (or deviance) tables for one or more fitted model objects.
stats::anova(pml_object)
stats::anova(opt_pml_object)

#e. Final Presentation - Visualizing the  best evolutionary tree
plot(opt_pml_object, main = "Maximum Likelihood Tree")
 
# *-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*--*-*-*-*-*-*-*-*-**-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
#5. Which of the 3 topologies best explains the data

#a.Compute lengths, then force the object from a list into a flat numeric vector
edge_weights_lowest_mut_tree <- ape::compute.brlen(lowest_mut_tree) #computes branch lengths of a tree using different methods

#b.Translate all the trees into likelihoods
#. pml computes the likelihood of a phylogenetic tree given a sequence alignment and a model
distance_likelihood  <- phangorn::pml(tree = NJ_tree, data = msa2, model = "WAG")
parsimony_likelihood <- phangorn::pml(tree = edge_weights_lowest_mut_tree, data = msa2, model = "WAG")
likelihood <- opt_pml_object 

#b. Run the Shimodaira-Hasegawa (SH) test across all topologies
# This tests the null hypothesis that all three trees fit the data equally well
topology_test <- phangorn::SH.test(distance_likelihood, parsimony_likelihood, likelihood)
print(topology_test)
