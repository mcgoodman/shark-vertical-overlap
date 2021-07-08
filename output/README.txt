
**median_depths.csv**
daily median depths, by individual, for each species

**depth_binned.csv**
summary of daily median depth distribution for each species, 
aggregated into 10 meter depth bins

**distance_matrices.rds**
Distance matrices (Bhattacharyya, Schoener's D, and Euclidian) used 
in hierarchical clustering analysis

**similarity_matrices.rds**
Similarity matrices (Bhattacharya & Schoener's D) summarizing niche 
similarity for each pair of species

**bhattacharyya_heatmap.png**
Heatmap showing depth distributions of each species in the water
column (up to 100 meters), and hierarchical clustering tree based
on Bhattacharrya's dissimilarity. Analagous plots are provided for
Schoener's D and Euclidian distance.

**bhattacharyya_heatmap_1500m.png**
Heatmap showing depth distributions of each species in the water
column (up to 1510 meters), and hierarchical clustering tree based
on Bhattacharrya's dissimilarity. Analagous plots are provided for
Schoener's D and Euclidian distance.

**bhattacharrya_matrix.png**
Pairwise heatmap displaying Bhattacharyya's coefficient for each
pair of species. Higher values indicate more niche overlap.

**cluster_number_bhattacharrya.png**
Plot of within-cluster sum of squares as a function of the number of 
clusters, used to select the number of clusters plotted in the 
heatmap/dendrogram plots. Analagous plots are provided for Schoener's 
D and Euclidian distance.