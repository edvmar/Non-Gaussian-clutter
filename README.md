# Non-Gaussian-clutter
Repository for project in the course MVE385 at Chalmers in collaboration with SAAB.

Students:

Edvin Martinson (edvmar), Erik Sahlin (Eklektikern) & Samuel Winqvist (samwin1234) 

Examiner at Chalmers: 

Alexey Geints

SAAB supervisors: 

Mika Persson & Andréas Sundström


# Special mathworks packages needed 

Statistics and Machine Learning Toolbox


# Important note: 

In the current task 2 implementation, the sampling is done by using a large list of tabulated values to get the inverse of the CDF (line 15, Task2/Sampling.m). This solution is very data ineffecient and we reccomend any potential further developer to change this to a solution using a function defined by interpolated data as the inverse CDF. It should not require a considerable effort to implement this change. 
