TEST CYL FMESH
1  400 -0.946 -1 IMP:n=1 IMP:p=1 IMP:E=1
2  0          +1 IMP:n=0 IMP:p=0 IMP:E=0

1 RPP -50 50 -50 50 -50 50

C ***********************************
C *  WATER
C *  MASS DENSITY [G/CC] - 0.946
C *  VOLUME FRACTION [%] - 100
C *
C ***********************************
C
C
M401     1001.21C 6.33910E-002 $ H  1 AMOUNT(%)  2.0000 AB(%) 99.99
         1002.21C 7.29080E-006 $ H  2 AMOUNT(%)  2.0000 AB(%)  0.01
         8016.21C 3.16221E-002 $ O 16 AMOUNT(%)  1.0000 AB(%) 99.76
C
M400     1001.21C 6.33910E-002 $ H  1 AMOUNT(%)  2.0000 AB(%) 99.99
         1002.21C 7.29080E-006 $ H  2 AMOUNT(%)  2.0000 AB(%)  0.01
         8016.21C 3.16221E-002 $ O 16 AMOUNT(%)  1.0000 AB(%) 99.76
C *
C *  T.A.D. = 9.50204E-002
C *  EFF.DENSITY = 9.46000E-001
C
sdef pos=0 0 0 erg=1
MODE N P
nps 1
c
C
FMESH4:N  GEOM=cyl ORIGIN=-50 -50 -50
           AXS=1 0 0 VEC=0 0 1
           IMESH=50  IINTS=3
           JMESH=50  JINTS=10
           KMESH=1  KINTS=10
           OUT=CF
C
FMESH14:N  GEOM=cyl ORIGIN=-50 -50 -50
           AXS=0 1 0 VEC=1 0 0
           IMESH=50  IINTS=3
           JMESH=50  JINTS=10
           KMESH=1  KINTS=10
           OUT=CF
C
FMESH24:N   GEOM=cyl ORIGIN=-50 -50 -50
           AXS=0 0 1 VEC=1 0 0
           IMESH=50  IINTS=3
           JMESH=50  JINTS=10
           KMESH=1  KINTS=10
           OUT=CF
C
FMESH124:N GEOM=XYZ ORIGIN=-50 -50 -50
           IMESH=50  IINTS=3
           JMESH=50  JINTS=3
           KMESH=50  KINTS=3
           OUT=CF
C