!     DGESV Example Program Text
!     NAG Copyright 2005.
!     .. Parameters ..
  Integer nin, nout
  Parameter (nin=5, nout=6)
  Integer nmax
  Parameter (nmax=8)
  Integer lda
  Parameter (lda=nmax)
!     .. Local Scalars ..
  Integer i, ifail, info, j, n
!     .. Local Arrays ..
  Double Precision a(lda, nmax), b(nmax)
  Integer ipiv(nmax)
!     .. External Subroutines ..
  External dgesv, x04caf
!     .. Executable Statements ..
  Write (nout, *) 'DGESV Example Program Results'
  Write (nout, *)
!     Skip heading in data file
  Read (nin, *)
  Read (nin, *) n
  If (n<=nmax) Then
!
!        Read A and B from data file
!
    Read (nin, *)((a(i,j),j=1,n), i=1, n)
    Read (nin, *)(b(i), i=1, n)
!
!        Solve the equations Ax = b for x
!
    Call dgesv(n, 1, a, lda, ipiv, b, n, info)
!
    If (info==0) Then
!
!           Print solution
!
      Write (nout, *) 'Solution'
      Write (nout, 100)(b(i), i=1, n)
!
!           Print details of factorization
!
      Write (nout, *)
      ifail = 0
      Call x04caf('General', ' ', n, n, a, lda, 'Details of factorization', &
        ifail)
!
!           Print pivot indices
!
      Write (nout, *)
      Write (nout, *) 'Pivot indices'
      Write (nout, 110)(ipiv(i), i=1, n)
!
    Else
      Write (nout, 120) 'The (', info, ',', info, ')', &
        ' element of the factor U is zero'
    End If
  Else
    Write (nout, *) 'NMAX too small'
  End If
  Stop
!
100 Format ((3X,7F11.4))
110 Format ((3X,7I11))
120 Format (1X, A, I3, A, I3, A, A)
End Program
