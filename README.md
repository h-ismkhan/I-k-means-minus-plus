# k-means-+

C++ Implementation of:

I-k-means−+: An iterative clustering algorithm based on an enhanced version of the k-means, Pattern Recognition Volume 79, July 2018, Pages 402-413

by: Hassan Ismkhan

Abstract

The k-means tries to minimize the sum of the squared Euclidean distance from the mean (SSEDM) of each cluster as its objective function. Although this algorithm is effective, it is too sensitive to initial centers. So, many approaches in the literature have focused on determining suitable initial centers. However, selecting suitable initial centers is not always possible, especially when the number of clusters is increased. This paper proposes an iterative approach to improve quality of the solution produced by the k-means. This approach tries to iteratively improve the quality of solution of the k-means by removing one cluster (minus), dividing another one (plus), and applying re-clustering again, in each iteration. This method called iterative k-means minus–plus (I-k-means−+). The I-k-means−+ is speeded up using some methods to determine which cluster should be removed, which one should be divided, and how to accelerate the re-clustering process. Results of experiments show that I-k-means−+ can outperform k-means++, to be known one of the accurate version of the k-means, in terms of minimizing SSEDM. For some instances, the accuracy of I-k-means−+ is about 2 times higher than both the k-means and k-means++, while it is faster than k-means++, and has the reasonable runtime, in comparison with the k-means.
