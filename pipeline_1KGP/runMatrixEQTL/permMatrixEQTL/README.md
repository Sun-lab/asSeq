1. Get an original gene fit*
2. Get coefficients estimates for the SNP using minimum p-value using lm model
3. Create several parametric bootstraps for different effect sizes: 0, 0.1,…0.9,1,2,…10 – since I used it for all genes I wanted to give estimates for less significant associations too.
4. For each simulated dataset run 1000 permutations
Timing: in this setup each permutation fit (of a 100 genes takes less then 1 minute)

*(MatrixEQTL is efficient when using multiple genes at a time, so I actually run them in batches by 100 both in original fit and in permutation setup)
