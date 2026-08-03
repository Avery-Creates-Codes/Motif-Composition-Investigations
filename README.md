# Motif-Composition-Investigations

cleaningFIMO.py: 
This script is used to remove FIMO's redundant putative hits. This script takes a Tomtom.tsv (generated from running MEME file against itself) and a FIMO.gff file as input and uses Tomtom similarities based on a user-specified threshold to group similar motifs. Motifs from the same group on the same strand that overlap are collapsed into the motif with the highest FIMO score.

Track1.ipynb: 
Track 1 is used as a baseline for clustering accessible chromatin regions (ACRs) based on motif count and binary matrices.

Track2.ipynb: 
Track 2 implements a TF-IDF motif representation for clustering ACRs.

Track3.ipynb: 
Track 3 evaluates different matrix representations that incorporate motif order, strand-awareness, and interval-awareness.
***Note - The computational cost to run the clustering and analyses on all 34 matrices is very high. It is recommended to edit the matrices = {} dictionary containing the 34 matrices down to the specific matrix representations desired by commenting out the undesired matrices from the dictionary. 

Track4.ipynb: 
Track 4 clusters ACRs based on a TF-IDF representation of k-merized sequences using k-mer lengths of 3, 5, and 7.

Cluster Evaluation CSVs:
The cluster evaluation CSVs summarize the clustering evaluation results, including the clustering method, number of clusters, largest cluster size, smallest cluster size, Silhouette score, Calinski-Harabasz Index, Davies-Bouldin Index, and Adjusted Rand Index (ARI) stability.

Biological CSVs:
The biological CSVs summarize the functional composition of clustering results. The CSVs include the clustering method, the Chi-square p-value for testing independence between functional annotation and the cluster assignment, the cluster size, Shannon Diversity Index, and the top three functional categories within each cluster along with their respective percentages.
