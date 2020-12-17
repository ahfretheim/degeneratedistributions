# degeneratedistributions
A new statistical method for measuring political polarization and breakdowns in the flow of ideas in society. Basically, degenerate distributions show you the degree to which factions (degeneracy) have formed in the larger population, and precisely where within society they tend to occur.

The theoretical set of all possible groups and their political opinions from the population is a series of binomial distributions, one for each distinct number of persons in a group (n), with p equal to the percentage of people in the population who support an idea. If groups were truly random, that is, not dependent upon peoples allegiances or attributes such as class, geography and education that would tend to determine their allegiances, we would expect actual groups of people that exist in the real world to follow this distribution. Degeneracy is a measure of the number of theoretical groups that are missing from this theoretical distribution of perfectly distributed ideas.

The mathematical definition is as follows:

D = sum(S*d[k, n] - A[k, n]) over all (see [2] in the implementation section for details on "all") k & n

Where S is the supremum (precise statement of use of the concept in this Python module to follow) of the distribution height within the data. Distribution height is defined as follows: d[k, n] is the density function of the Binomial distribution of probability p and repeat count n, with respect to having precisely k "yes's". A[k, n] is the actual number of groups in the data of size n with k yes's.

Specific to the implementation of this definition by this code:

0. The language is Python 3. Required libarires include SciPy, csv & math. Later versions will use Rpy2 & ggplot2 or matplotlib for visualization functionalities, but this is not yet implemented in the present version. See bullet (3) for formatting requirements of the data file.

1. The Supremum is determined as the largest ratio of A[k, n]/d[k,n] from the data that has a z-score with less than or equal to 2. The reason for the z-score cleaning is to eliminate outliers. For a sufficiently large dataset, this cleaning step could be removed and a simple maximum used, although at this time I do not know precisely how large would be large enough to achieve that goal. The odd nature of this statistic makes it a difficult question to consider.

2.The code only sums up nonmodal [k, n]'s, which both speeds up the performance and assists with the Supremum outlier cleansing mentioned in bullet 1 above. Trivially, the D contribution for the supremum would always be zero, as A[k,n] = S*d[k,n] for the supremum, by definition. We extend this trivially cancelling to all possible supremums. Lastly, any k's and n's that NEVER occur in your data are not calculated - a similar exception to the definitional logic is described in bullet 5 pertaining to the comparative chi-squared test. 

3. The degenerateBinomial class takes a .csv file with your group-based voter data, as well as optional inputs indicating the name of the "yes" column and the name of the "Group ID" column in the file. There is also an optional input for the true population level of support for the proposition, which is crucial if your dataset is not fully representative of the population. It returns a total degeneracy of the inputted data, and a chi-squared statistic indicating the p-value of your observed groups belonging to an ideal distribution. This is implemented so you can compare the calculated degeneracy with a more common statistical approach, for the purpose of reconciling with the existing literature.

4. The degenerateBionomial class also cotains a python list of the term-wise contributions - one of the advantages of this method over others is it can show you not only how different the true distribution is from the Binomial and how unlikely to be Binomial, but where the breakdown in the similarities occurs. This, in turn, allows the studious sociologist to determine the full landscape of polarization and where the breakdowns are occurring. Later versions (probably the next version) will take a "ilk" column for describing the type of group the individual observed groups belong to, which will assist with regressoin analysis of the nature of the observed degeneracy. Later versions will also include graphing and visualization functionalities, and more thorough summaries.

5. The chi-squared functionality omits, per the reccomendation of the SciPy-statistics literature, any cases where either the observed or actual number of occurrences is too small for the chi-squared statistic to be correctly calculated

CURRENT VERSION: 0.2 BETA
