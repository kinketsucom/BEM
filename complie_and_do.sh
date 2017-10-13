cd `dirname $0`
gfortran ./bem.f90 -o do_bem -L. -llapack -lblas;
./do_bem;
