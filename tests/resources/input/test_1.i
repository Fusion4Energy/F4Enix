Problem fv5  -  15 x 15 assembly
c
c    CELLS
10   100  0.06925613  -1      u=1      imp:n=1   $ fuel
20   200  0.042910     1 -2   u=1      imp:n=1   $ clad
30   300  0.100059     2      u=1      imp:n=1   $ water
40   0                -3    u=2 fill=1  lat=1    imp:n=1
50   0                -4  #80   fill=2    imp:n=1   $ lattice region
60   300  0.100059     4 -5            imp:n=1   $ outer water
70   400  0.083770     5 -6            imp:n=1   $ iron
80   300  -1       2 -3 5 4  fill=3 imp:n=1 
90   400 5 6     u=3 imp:n=1

c    SURFACES
1    cz   0.44
2    cz   0.49
3    RPP    -.7   .7     -.7   .7   0. 0.
4    RPP  -10.5 10.5   -10.5 10.5   0. 0.
5    RPP  -13.0 13.0   -13.0 13.0   0. 0.
*6   RPP  -13.5 13.5   -13.5 13.5   0. 0.

c    DATA
c
kcode   1000    1.0  10  50 
hsrc  15 -10.5 10.5   15 -10.5 10.5   1 -9e9 9e9
sdef   x=d1  y=d2  z=0.0
si1   -10.5  10.5
sp1      0     1
si2   -10.5  10.5
sp2      0     1
c
m100   92238 2.2380e-2   92235 8.2213e-4  8016 4.6054e-2  $ fuel
m200   40000 4.2910e-2                                    $ clad
m300    1001 6.6706e-2    8016 3.3353e-2                  $ water
mt300   lwtr
m400   26000 0.083770
