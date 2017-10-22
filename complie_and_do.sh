cd `dirname $0`
gfortran ./newbem.f90 -o do_new_bem -L. -llapack -lblas;
echo "#Fin#:compile"
./do_new_bem;
echo "#Fin#:execute"
# gfortran ./bem.f90 -o do_bem -L. -llapack -lblas;
# gnuplot
# gnuplot <<- EOF
# "";
# splot 'data.dat';
# EOF
# gfortran ./test.f90 -o do_test -L. -llapack -lblas;
# ./do_test;
