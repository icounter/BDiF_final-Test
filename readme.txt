
sic idea here Jarque–Bera test:
In statistics, the Jarque–Bera test is a goodness-of-fit test of whether sample data have the skewness and kurtosis matching a normal distribution.
So my idea here is divide this file into several parts. We do this test for every part. This, I believe, would be efficient and accurate when the number of data goes to larger and larger.
For every thread, we first calculate the mean, the variance, since both in the function of skewness and kurtosis need mean and variance. After calculate the skewness and kurtosiss, we use the formula to calculate the statistic JB. As we know, this statistic follow chi-squared distribution. Then we compare the value with the distribution table to decide whether this part is normally distributed or not.

