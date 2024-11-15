MCNP XS LEAK SPHERE TEST INPUT VRT
1  0          -1    IMP:N=1   IMP:P=1
2  1 -0.302  +1 -2     IMP:N=1   IMP:P=1
3  0          +2    IMP:N=0   IMP:P=0

1 S 0 0 0 5
2 S 0 0 0 50
3 S 0 0 0 60
C

C
SDEF POS 0 0 0 PAR=N ERG=d1 $ position, particle type, energy
SI1 H 1e-6 0.1 1 10 14      $ histogram boundaries
SP1 D 0    1   1 1  1       $ probabilities for each bin
C
MODE N P
PHYS:P J 1  
c
PRDMP  2J  -1 $ Flag to print the mctal
C
FC124 Neutron Flux mesh
FMESH124:N  ORIGIN=-50 -50 -50
            IMESH=50 IINTS=10
            JMESH=50 JINTS=10
            KMESH=50 KINTS=10