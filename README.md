# Motif-Composition-Investigations

cleaningFIMO.py: 
This script is used to remove FIMO's redundant putative hits. This script takes a Tomtom.tsv (generated from running MEME file against itself) and a FIMO.gff file as input and uses Tomtom similarities based on a user-specified threshold to group similar motifs. Motifs from the same group on the same strand that overlap are collapsed into the first motif that appears in the overalpping coordinates.

fimoGFF_Exploration.ipynb:
This notebook is used to explore the cleaned FIMO file in regards to the distribution of FIMO hits amongst the ACRs as well as to determine the number of ACRs with zero motif hits.

Track1.ipynb: 
Track 1 is used as a baseline for clustering accessible chromatin regions (ACRs) based on motif count and binary matrices.

Track2.ipynb: 
Track 2 implements a TF-IDF motif representation for clustering ACRs.

Track3.ipynb: 
Track 3 evaluates different matrix representations that incorporate motif order, strand-awareness, and interval-awareness.
***Note - The computational cost to run the clustering and analyses on all 34 matrices is very high. It is recommended to edit the matrices = {} dictionary containing the 34 matrices down to the specific matrix representations desired by commenting out the undesired matrices from the dictionary. 

Track4.ipynb: 
Track 4 clusters ACRs based on a TF-IDF representation of k-merized sequences using k-mer lengths of 3, 5, and 7.

Clean_VS_Redundant_FIMO.ipynb:
Compares the evaluation metrics of clusters between the cleaned FIMO file and the original redundant FIMO file.

ZeroHitExclusionClustering.ipynb:
Compares evaluation metrics between clustering all ACRs (including ACRs that contain zero motif hits) and clustering only ACRs that contain motif hits.

Cluster Evaluation CSVs:
The cluster evaluation CSVs summarize the clustering evaluation results, including the clustering method, number of clusters, largest cluster size, smallest cluster size, Silhouette score, Calinski-Harabasz Index, Davies-Bouldin Index, and Adjusted Rand Index (ARI) stability.
