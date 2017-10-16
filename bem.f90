! boundary element method
program bem
  !宣言部
  integer N
  real R
  REAL, PARAMETER :: pi = 3.14159
	double precision,allocatable :: end_point(:,:)
	double precision,allocatable :: ans(:)
	double precision,allocatable :: b_vec(:)
	double precision,allocatable :: P(:,:)
	double precision,allocatable :: U(:,:)
	double precision,allocatable :: W(:,:)
	double precision,dimension(1,1) :: t_vec(2,1)=(0,0)
	double precision,dimension(1,1) :: n_vec(2,1)=(0,0)
	double precision,dimension(1,1) :: m_1_vec(2,1)=(0,0)
	double precision,dimension(1,1) :: m_2_vec(2,1)=(0,0)
	double precision x,y,r_1,r_2,h,theta
  Integer,allocatable :: ipiv(:)
  !ここからはテストケース
	double precision,allocatable :: test_ans(:)
  double precision input_x,input_y,length,ans_u,hankei

  !実行部
  ! print *, "Please enter integers:R="
  ! read *, R
  ! めんどいから今は10固定
  R = 10
  ! print *, "Please enter integers:N="
  ! read *, N
  N= 100

  print *, "N:",N,"Radius:",R

  allocate( end_point(N,2) )
  allocate( P(N,2) )
  allocate( U(N,N) )
  allocate( W(N,N) )
  allocate( ans(N) )
  allocate( test_ans(N) )
  allocate( ipiv(N) )
  allocate( b_vec(N) )

  do i=1, N
    ! 分割端点の設定
    end_point(i,1) = R*cos(2*pi*(i-1)/N)
    end_point(i,2) = R*sin(2*pi*(i-1)/N)
  enddo
  do i =1,N
    ! 分割代表点の設定
    if (i/=N) then
      P(i,1) = (end_point(i,1)+end_point(i+1,1))/2
      P(i,2) = (end_point(i,2)+end_point(i+1,2))/2
    else if (i==N) then
      P(i,1) = (end_point(i,1)+end_point(1,1)) / 2
      P(i,2) = (end_point(i,2)+end_point(1,2)) / 2
    end if
    ans(i) = P(i,1)**3 - 3*P(i,1)*(P(i,2)**2) !u(Pj)のこと
    test_ans(i) = ( 3*(P(i,1)**3) - 9*P(i,1)*(P(i,2))**2 )/R !微分u(Pj)の理論値的なあれ
  enddo

! Uijをつくるゾーーイいまはきったなくてもいい
  do i = 1,N
    do j = 1,N
      ! r_jをつくってみた
      if (j/=N) then
        m_1_vec(1,1) = P(i,1)- end_point(j,1)
        m_1_vec(2,1) = P(i,2)- end_point(j,2)
        m_2_vec(1,1) = P(i,1)- end_point(j+1,1)
        m_2_vec(2,1) = P(i,2)- end_point(j+1,2)
        h = sqrt( (end_point(j,1)- end_point(j+1,1))**2 + (end_point(j,2)- end_point(j+1,2))**2)
        t_vec(1,1) = (end_point(j+1,1)- end_point(j,1))/h
        t_vec(2,1) = (end_point(j+1,2)- end_point(j,2))/h
      else if (j==N) then
        m_1_vec(1,1) = P(i,1)- end_point(j,1)
        m_1_vec(2,1) = P(i,2)- end_point(j,2)
        m_2_vec(1,1) = P(i,1)- end_point(1,1)
        m_2_vec(2,1) = P(i,2)- end_point(1,2)
        r_1 = sqrt( m_1_vec(1,1)**2 + m_1_vec(2,1)**2)
        r_2 = sqrt( m_2_vec(1,1)**2 + m_2_vec(2,1)**2)
        h = sqrt((end_point(j,1)- end_point(1,1))**2 + (end_point(j,2)- end_point(1,2))**2)
        t_vec(1,1) = (end_point(1,1)- end_point(j,1))/h
        t_vec(2,1) = (end_point(1,2)- end_point(j,2))/h
      end if

      !n_vecの計算
      n_vec(1,1) = -t_vec(2,1)
      n_vec(2,1) = t_vec(1,1)
      !r_1とr_2の計算
      r_1 = sqrt( m_1_vec(1,1)**2 + m_1_vec(2,1)**2)
      r_2 = sqrt( m_2_vec(1,1)**2 + m_2_vec(2,1)**2)
      !thetaの計算
      if( abs(m_2_vec(1,1)*n_vec(1,1)+m_2_vec(2,1)*n_vec(2,1)) < 0.0001 .and. &
        abs( m_1_vec(1,1)*n_vec(1,1)+m_1_vec(2,1)*n_vec(2,1)) < 0.0001 ) then
        theta = abs(atan2( m_2_vec(1,1)*t_vec(1,1)+m_2_vec(2,1)*t_vec(2,1), m_2_vec(1,1)*n_vec(1,1)+m_2_vec(2,1)*n_vec(2,1))  &
              - atan2( m_1_vec(1,1)*t_vec(1,1)+m_1_vec(2,1)*t_vec(2,1), m_1_vec(1,1)*n_vec(1,1)+m_1_vec(2,1)*n_vec(2,1) ))
      else
        theta = atan2( m_2_vec(1,1)*t_vec(1,1)+m_2_vec(2,1)*t_vec(2,1) , m_2_vec(1,1)*n_vec(1,1)+m_2_vec(2,1)*n_vec(2,1)) &
              - atan2( m_1_vec(1,1)*t_vec(1,1)+m_1_vec(2,1)*t_vec(2,1), m_1_vec(1,1)*n_vec(1,1)+m_1_vec(2,1)*n_vec(2,1))
      end if
      !Uの計算
      U(i,j) =( (m_2_vec(1,1)*t_vec(1,1)+m_2_vec(2,1)*t_vec(2,1))*log(r_2) &
                            - (m_1_vec(1,1)*t_vec(1,1)+m_1_vec(2,1)*t_vec(2,1))*log(r_1) &
                            + h - (m_1_vec(1,1)*n_vec(1,1)+m_1_vec(2,1)*n_vec(2,1))*theta ) /(2*pi)

      W(i,j) = theta / (2*pi)

    end do
  end do

  !Wuの計算
  do i = 1,N
    do j = 1,N
      b_vec(i) = b_vec(i) + W(i,j)*ans(j)
    end do
  end do

  Call DGESV(N, 1, U, N, ipiv, b_vec, N, info)

  print *,"TestSolition"
  write(*, *) test_ans(1:4)
  print *,"Solition"
  write(*, *) b_vec(1:4) !Nこの解の一部

  ! 出力するためのうんコード
  input_x = 0.0
  input_y = 0.0
  ! ファイル書き込み
  open(1, file='data.dat', status='replace')
  do i = 1,N
    do j = 1,N
        input_x = 2*R*i/N -R
        input_y = 2*R*j/N -R
        ans_u = 0.0
        length = sqrt( input_x**2 + input_y**2 )
        if(length >= R) then
          ans_u = 0.0
        else
          !　内点計算
          do k = 1,N
              hankei = sqrt( (input_x - end_point(k,1))**2 + (input_y - end_point(k,2))**2 )
              ans_u = ans_u -log(hankei)*b_vec(k)/(2*pi)
              ans_u = ans_u - ( (input_x-end_point(k,1))*end_point(k,1)/R + (input_y - end_point(k,2))*end_point(k,2)/R ) &
                          / (2*pi*hankei**2) * ans(k)
          enddo
        endif
        ans_u = ans_u * 2*pi*R / N
        write (1,*) input_x," ",input_y," ",ans_u
    end do
  enddo

  close(1)

!



end program bem
