#!/usr/bin/env bash

FFTW3FLAGS=" "

echo ""
echo "Verification test"
echo ""
for i in nosplit_transpose nosplit_fft_physical_in_x nosplit_fft_physical_in_z
do
   mpif90 -cpp -O3 -march=native ${i}.f90 -I../mod -L../ ${FFTW3FLAGS} -ldecomp2d -o ${i}
   mpirun -n 4 ./${i}
done

echo ""
echo "Performance test"
echo ""
for i in nosplit
do
   mpif90 -cpp -O3 -march=native ${i}.f90 -I../mod -L../ ${FFTW3FLAGS} -ldecomp2d -o ${i}
   mpirun -n 4 ./${i}
done
