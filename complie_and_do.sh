cd `dirname $0`
gfortran ./bem.f90 -o do_bem -L. -llapack -lblas;
echo "#Fin#:compile"
./do_bem;
echo "#Fin#:execute"
# gnuplot
# gnuplot <<- EOF
# "";
# splot 'data.dat';
# EOF
# gfortran ./test.f90 -o do_test -L. -llapack -lblas;
# ./do_test;
