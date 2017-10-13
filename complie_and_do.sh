cd `dirname $0`
gfortran ./bem.f90 -o do_bem -L. -llapack -lblas;
./do_bem;
# gfortran ./test.f90 -o do_test -L. -llapack -lblas;
# ./do_test;
