! boundary element method
program bem
  !宣言部
  integer N
  double precision R
  double precision, PARAMETER :: pi = 3.14159265359
	double precision,allocatable :: end_point(:,:)
	double precision,allocatable :: ans(:)
	double precision,allocatable :: b_vec(:)
	double precision,allocatable :: P(:,:)
	double precision,allocatable :: U(:,:)
	double precision,allocatable :: W(:,:)
  double precision,allocatable :: U2(:,:)
  double precision,allocatable :: W2(:,:)
  double precision,allocatable :: theta1(:,:)
  double precision,allocatable :: theta2(:,:)
	double precision,dimension(1,1) :: t_vec(2,1)=(0,0)
	double precision,dimension(1,1) :: n_vec(2,1)=(0,0)
	double precision,dimension(1,1) :: m_1_vec(2,1)=(0,0)
	double precision,dimension(1,1) :: m_2_vec(2,1)=(0,0)
	double precision x,y,r_1,r_2,h,theta
  Integer,allocatable :: ipiv(:)
  double precision,allocatable :: point_x1_vec(:,:)
  double precision,allocatable :: point_x2_vec(:,:)

  !ここからはテストケース
	double precision,allocatable :: test_ans(:)
  double precision,allocatable :: input(:,:)
  double precision input_x,input_y,length,ans_u,hankei
  Integer mesh_num,count
  double precision kai,err_sum
  !実行部
  ! print *, "Please enter integers:R="
  ! read *, R
  ! めんどいから今は10固定
  R = 10
  ! print *, "Please enter integers:N="
  ! read *, N
  N= 100
  mesh_num = 50


  print *, "N:",N,"Radius:",R

  allocate( end_point(N,2) )
  allocate( P(N,2) )
  allocate( U(N,N) )
  allocate( W(N,N) )
  allocate( U2(mesh_num**2,N) )
  allocate( W2(mesh_num**2,N) )
  allocate( point_x1_vec(max(mesh_num,N),2) )
  allocate( point_x2_vec(max(mesh_num,N),2) )
  allocate( theta1(N,1))
  allocate( theta2(mesh_num**2,1))
  allocate( ans(N) )
  allocate( test_ans(N) )
  allocate( ipiv(N) )
  allocate( b_vec(N) )
  allocate( input(mesh_num**2,2) )

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
    ans(i) = funcU( P(i,1), P(i,2) ) !u(Pj)のこと
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
      if( abs(m_2_vec(1,1)*n_vec(1,1)+m_2_vec(2,1)*n_vec(2,1)) < 0.0001 .or. &
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


  do i = 1,N
    err_sum = err_sum + abs(test_ans(i)-b_vec(i))
  end do

  ! print *,"TestSolition"
  ! write(*, *) test_ans(1:4)
  ! print *,"Solition"
  ! write(*, *) b_vec(1:4) !Nこの解の一部

  ! open(1, file='q_calc.dat', status='replace')
  !
  !   do i = 1,N
  !     write (1,*) i, " " , b_vec(i)
  !   end do
  ! close(1)
  ! open(1, file='q_ans.dat', status='replace')
  !
  !   do i = 1,N
  !     write (1,*) i, " " , test_ans(i)
  !   end do
  ! close(1)


  print * , "err_sum:",err_sum/N

  ! 出力するためのうんコード
  count = 0
  do i = 1,mesh_num**2
      ! input(i,1) = (R*((i/mesh_num)-1))/(mesh_num+1)*cos(2*pi*(i-1)/mesh_num)
      ! input(i,2) = (R*((i/mesh_num)-1))/(mesh_num+1)*sin(2*pi*(i-1)/mesh_num)
    !inputを内点としてとる


    count = count + 1
      input(i,1) = 2*R*(i/mesh_num +1 )/mesh_num - R
    if (mod(count,mesh_num) == 0 ) then
      count = 1
    end if
      input(i,2) = 2*R*count/mesh_num - R
        ! ans_u = 0.0
        ! length = sqrt( input_x**2 + input_y**2 )
        ! if(length >= R) then
        !   ans_u = 0.0
        ! else
        !   !　内点計算
        !   do k = 1,N
        !       hankei = sqrt( (input_x - end_point(k,1))**2 + (input_y - end_point(k,2))**2 )
        !       ans_u = ans_u -log(hankei)*b_vec(k)/(2*pi)
        !       ans_u = ans_u - ( (input_x-end_point(k,1))*end_point(k,1)/R + (input_y - end_point(k,2))*end_point(k,2)/R ) &
        !                   / (2*pi*hankei**2) * ans(k)
        !   enddo
        ! endif
        ! ans_u = ans_u * 2*pi*R / split_num
        ! write (1,*) input_x," ",input_y," ",ans_u
  end do


  !もういっぱつ
  ! U mesh,n
    do i = 1,mesh_num**2
      do j = 1,N
        ! r_jをつくってみた
        if (j/=N) then
          m_1_vec(1,1) = input(i,1) - end_point(j,1)
          m_1_vec(2,1) = input(i,2) - end_point(j,2)
          m_2_vec(1,1) = input(i,1) - end_point(j+1,1)
          m_2_vec(2,1) = input(i,2) - end_point(j+1,2)
          h = sqrt( (end_point(j,1)- end_point(j+1,1))**2 + (end_point(j,2)- end_point(j+1,2))**2 )
          t_vec(1,1) = (end_point(j+1,1)- end_point(j,1)) / h
          t_vec(2,1) = (end_point(j+1,2)- end_point(j,2)) / h

        else if (j==N) then
          m_1_vec(1,1) = input(i,1) - end_point(j,1)
          m_1_vec(2,1) = input(i,2) - end_point(j,2)
          m_2_vec(1,1) = input(i,1) - end_point(1,1)
          m_2_vec(2,1) = input(i,2) - end_point(1,2)
          r_1 = sqrt( m_1_vec(1,1)**2 + m_1_vec(2,1)**2)
          r_2 = sqrt( m_2_vec(1,1)**2 + m_2_vec(2,1)**2)

          h = sqrt((end_point(j,1)- end_point(1,1))**2 + (end_point(j,2)- end_point(1,2))**2)
          t_vec(1,1) = (end_point(1,1)- end_point(j,1)) / h
          t_vec(2,1) = (end_point(1,2)- end_point(j,2)) / h
        end if

        !n_vecの計算
        n_vec(1,1) = -t_vec(2,1)
        n_vec(2,1) = t_vec(1,1)
        !r_1とr_2の計算
        r_1 = sqrt( m_1_vec(1,1)**2 + m_1_vec(2,1)**2 )
        r_2 = sqrt( m_2_vec(1,1)**2 + m_2_vec(2,1)**2 )
        ! print * , i,j,m_1_vec
        !thetaの計算
        if( abs(m_2_vec(1,1)*n_vec(1,1)+m_2_vec(2,1)*n_vec(2,1)) < 0.0001 .or. &
          abs( m_1_vec(1,1)*n_vec(1,1)+m_1_vec(2,1)*n_vec(2,1)) < 0.0001 ) then
          theta = abs(atan2( m_2_vec(1,1)*t_vec(1,1)+m_2_vec(2,1)*t_vec(2,1), m_2_vec(1,1)*n_vec(1,1)+m_2_vec(2,1)*n_vec(2,1))  &
                - atan2( m_1_vec(1,1)*t_vec(1,1)+m_1_vec(2,1)*t_vec(2,1), m_1_vec(1,1)*n_vec(1,1)+m_1_vec(2,1)*n_vec(2,1) ))
        else
          theta = atan2( m_2_vec(1,1)*t_vec(1,1)+m_2_vec(2,1)*t_vec(2,1) , m_2_vec(1,1)*n_vec(1,1)+m_2_vec(2,1)*n_vec(2,1)) &
                - atan2( m_1_vec(1,1)*t_vec(1,1)+m_1_vec(2,1)*t_vec(2,1), m_1_vec(1,1)*n_vec(1,1)+m_1_vec(2,1)*n_vec(2,1))
        end if
        !Uの計算

        U2(i,j) =( (m_2_vec(1,1)*t_vec(1,1)+m_2_vec(2,1)*t_vec(2,1))*log(r_2) &
                              - (m_1_vec(1,1)*t_vec(1,1)+m_1_vec(2,1)*t_vec(2,1))*log(r_1) &
                              + h - (m_1_vec(1,1)*n_vec(1,1)+m_1_vec(2,1)*n_vec(2,1))*theta ) /(2*pi)
        W2(i,j) = theta / (2*pi)
        ! print * , i,j,":i,j",U2(i,j),r_1,r_2
      end do
    end do


! print *, "U2"
! print *, U2
! print *, "W2"
! print *, W2

! ファイル書き込み
open(1, file='data.dat', status='replace')
err_sum=0.0
do i = 1,mesh_num**2
      ans_u = 0.0
      length = sqrt( input(i,1)**2 + input(i,2)**2 )
      if(length < R) then
        ! 内点計算
          do j = 1,N
            ans_u = ans_u + U2(i,j)*b_vec(j) - W2(i,j)*funcU( end_point(j,1), end_point(j,2) )
          end do
          write (1,*) input(i,1)," ",input(i,2)," ",ans_u
          kai = funcU( input(i,1),input(i,2) )
          err_sum = err_sum + abs(kai-ans_u)
      end if
end do

! print *, "U2",U2(13,:)


print *,"err_sum:", err_sum/mesh_num**2
! print * , "x:y"
! do i = 1, mesh_num**2
!   print * ,input(i,1),input(i,2)
! end do

close(1)


end program bem


function funcU(x,y)
  double precision x,y
  funcU = x**3-3*x*y**2
end function
