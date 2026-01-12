## Shannon Entropy of Conformational States

To quantify the conformational diversity sampled by each peptide tag, we computed the Shannon entropy of cluster populations obtained from RMSD-based clustering of molecular dynamics trajectories.

Peptide-tag trajectories from all nine chains and three independent trials were concatenated into a single ensemble trajectory. Conformations were then grouped into clusters using a fixed RMSD cutoff. For a given system, the probability of observing cluster *i* was defined as:

```math
p_i = \frac{N_i}{\sum_{j} N_j}
```

where \(N_i\) is the number of trajectory frames assigned to cluster *i*, and the denominator represents the total number of frames across all clusters.

The Shannon entropy was calculated as:

```math
H = -\sum_i p_i \log_2 p_i
```

where $\mathbf{H}$ is the Shannon entropy, and the summation runs over all identified clusters.

### Interpretation

In the context of molecular dynamics simulations, Shannon entropy provides a compact, ensemble-level measure of conformational freedom:

- Higher entropy indicates that the peptide samples many distinct conformations with comparable probabilities, reflecting increased structural heterogeneity and flexibility.
- Lower entropy indicates that the peptide occupies a smaller number of dominant conformational states, consistent with greater compactness or structural stabilization.

When interpreted alongside cumulative distribution functions (CDFs) of cluster sizes and end-to-end distance distributions, this metric enables direct comparison of sequence-dependent conformational landscapes within the polymer-brush environment.
