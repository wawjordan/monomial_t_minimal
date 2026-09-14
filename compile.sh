#!./bin/bash
. clean_local.sh
export GFORTRAN_DEBUG_FLAGS="-g3 -Wall -Wextra -fimplicit-none -fcheck=all -ffpe-trap=zero,underflow,overflow,invalid --coverage -fbacktrace -fstack-protector-all -fstack-check -Wimplicit-procedure -Wno-unused-parameter -Wno-unused-variable -Wno-unused-dummy-argument"
gfortran $GFORTRAN_DEBUG_FLAGS monomial_t_minimal.F90
