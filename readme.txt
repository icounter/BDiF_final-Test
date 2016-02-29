This is the basic idea of my work:
I use MPI to first put data into several processes(use MPI_file_read) and how much content should each process is decided by (filesize/number of preocesses).I store the content into a char array. Then I defined a struct to read the char array line by line for every 10000 lines. If the data's price or volume is smaller than 0, I define them as a noise. Also, if the dates are different from the most date in every 10000 lines, I define them as noise data. I put noise and valid both to two vectors.And output them use MPI_allgather Method and MPI_file_write_at_all method.
basic idea here normality test:
I combine these two files together and write the normality test in the .hpp file.
This can reduce the work for read data again.
Also in statistics, the Jarque–Bera test is a goodness-of-fit test of whether sample data have the skewness and kurtosis matching a normal distribution.
For every process, we first calculate the mean, the variance, since both in the function of skewness and kurtosis need mean and variance. After calculate the skewness and kurtosiss, we use the formula to calculate the statistic JB. As we know, this statistic follow chi-squared distribution. Then we compare the value with the distribution table to decide whether this part is normally distributed or not.

