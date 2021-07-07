**Steps: Clustering of vertical niches**

1. Daily median depths for each species / individual were computed, and the data were subsetted to species with 5 or more observed individuals.

2. Depth data were binned into 10 meter intervals, summarizing the number of median daily depths in each interval from 0-10 meters to 1500-1510 meters (the maximum depth observed for all species). For each species, the proportion of daily median depths in each depth bin $i$ ("proportion time at depth", $\text{p}_i$) was computed.

3. For each pair of species, Bhattacharyya's coefficient was computed using the proportion time at depth. For two species $x$ and $y$ and 151 10-meter depth bins: 
   $$
   B_{x, y} = \sum_{i = 1}^{151} \sqrt{\text{p}_i(x) \times \text{p}_i(y)}
   $$
   Bhattacharyya's coefficient ranges from 0 to 1 and measures the similarity between two discrete probability distributions, where 0 indicates no overlap between distributions, and 1 indicates identical depth distributions.

   Bhattacharyya's coefficient was converted to a dissimilarity by subtracting the value of the coefficient from 1:
   $$
   \tilde{B}_{x, y} = 1 -  B_{x, y} = 1 - \sum_{i = 1}^{151} \sqrt{\text{p}_i(x) \times \text{p}_i(y)}
   $$

4. Hierarchical clustering analysis was performed using Ward's minimum variance method.

5. Following cluster analysis, species were assigned to 6 discrete clusters. The number of clusters was chosen by plotting the within-cluster some of squares against the number of clusters.