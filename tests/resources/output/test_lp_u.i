MCNP XS LEAK SPHERE TEST INPUT VRT
1  0          -1    IMP:N=1
2  1 -0.0708  +1 -2 +4 +5     IMP:N=1 
3  0          +2    IMP:N=0 
4  0  -4  FILL=1 IMP:N=1
5  0  -5  FILL=2 IMP:N=1
6  1 -0.07 -7 6 U=1 IMP:N=1
7  1 -0.07 -9 8 U=2 IMP:N=1

1 S 0 0 0 5
2 S 0 0 0 50
3 S 0 0 0 60
4 S 10 10 10 10
5 S -10 -10 -10 10
6 S 10 10 10 5
7 S 10 10 10 11
8 S -10 -10 -10 1
9 S -10 -10 -10 11
C

M1
C hydrogen H-1
       1001.31c          1.0     $        AB(%)     
C
SDEF POS 0 0 0 PAR=N ERG=d1 $ position, particle type, energy
SI1 H 1e-6 0.1 1 10 14      $ histogram boundaries
SP1 D 0    1   1 1  1       $ probabilities for each bin
C
MODE N
LOST 1e3
c
PRDMP  2J  -1 $ Flag to print the mctal
C
FC34 DPA production
F34:N 2
FM34 1 1 444
SD34 1
NPS 1e4
