%compiles mex files
%other options:  for CFLAGS -fPIC -fopenmp -ftree-vectorize
%                for LDFLAGS -fopenmp
%to see results of vectorization use -ftree-loop-optimized=log.out etc.

mex relaxVEL3Dper3_mex.c            CFLAGS="\$CFLAGS -O3 -std=c99" LDFLAGS="\$LDFLAGS "
mex dirac_interp_mex.c              CFLAGS="\$CFLAGS -O3 -std=c11" LDFLAGS="\$LDFLAGS "
mex residualvel3Dper3_mex.c         CFLAGS="\$CFLAGS -O3 -std=c99" LDFLAGS="\$LDFLAGS "
mex AonVecVEL3Dper2_mex.c           CFLAGS="\$CFLAGS -O3 -std=c99" LDFLAGS="\$LDFLAGS "
mex relaxPRESSUREprodRB3Dper2_mex.c CFLAGS="\$CFLAGS -O3 -std=c99" LDFLAGS="\$LDFLAGS "
mex restricthto2h3DVper2_mex.c      CFLAGS="\$CFLAGS -O3 -std=c99" LDFLAGS="\$LDFLAGS "