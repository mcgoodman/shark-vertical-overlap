**Steps: Clustering of vertical niches**

1. Depth data were binned into 10 meter intervals for each individual of each species, summarizing the proportion of time spent ("proportion time at depth", $\text{p}_i$) by from 0-10 meters to 1840-1850 meters. For each species with five or more individuals, these proportions were averaged across individuals to obtain one distribution of proportion time at depth.

3. For each pair of species, Bhattacharyya's coefficient was computed using the proportion time at depth. For two species $x$ and $y$ and $n$ 10-meter depth bins: 
   $$
   B_{x, y} = \sum_{i = 1}^{n} \sqrt{\text{p}_i(x) \times \text{p}_i(y)}
   $$
   Bhattacharyya's coefficient ranges from 0 to 1 and measures the similarity between two discrete probability distributions, where 0 indicates no overlap between distributions, and 1 indicates identical depth distributions.

   Bhattacharyya's coefficient was converted to a dissimilarity by subtracting the value of the coefficient from 1:
   $$
   \tilde{B}_{x, y} = 1 -  B_{x, y} = 1 - \sum_{i = 1}^{n} \sqrt{\text{p}_i(x) \times \text{p}_i(y)}
   $$

4. Hierarchical clustering analysis was performed using Ward's minimum variance method.

4. Following cluster analysis, species were assigned to 4 discrete clusters. The number of clusters was chosen by plotting the within-cluster some of squares against the number of clusters: 

   <img src=".\output\cluster_number_bhattacharyya.png" alt="cluster_number_bhattacharyya" style="zoom:50%;" />