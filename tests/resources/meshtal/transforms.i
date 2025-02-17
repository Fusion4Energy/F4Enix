Title
1 1 -1 -1 imp:n=1
2 0 1 imp:n=0

1 so 1

M1 1001.31c 1.0
nps 1
sdef pos = 0 0 0
TR96 570.0000 50.0000 390.0000 
     1.0000 0.0000 0.0000 
     0.0000 -1.0000 0.0000 
     0.0000 0.0000 1.0000
C
*TR1    0  0  0
        20.0000    70.0000  90
       110.0000    20.0000  90
       90        90         0
TR2    1  1  1
*TR4  10 10 10
       10 100 90
       80 10 90
       90 90 0
FMESH2024:n  origin = -3 -3 -3
     imesh 3 iints 3 
     jmesh 3 jints 3 
     kmesh 3 kints 3
     tr=96
FMESH2124:n  origin = -3 -3 -3
     imesh 3 iints 3 
     jmesh 3 jints 3 
     kmesh 3 kints 3
     tr=1
FMESH2224:n  origin = -3 -3 -3
     imesh 3 iints 3 
     jmesh 3 jints 3 
     kmesh 3 kints 3
     tr=2
FMESH2324:n  origin = -3 -3 -3
     imesh 3 iints 3 
     jmesh 3 jints 3 
     kmesh 3 kints 3
     tr=4
FMESH2424:n  origin = -3 -3 -3
     imesh 3 iints 3 
     jmesh 3 jints 3 
     kmesh 3 kints 3
     out=cuv
