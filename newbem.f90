! boundary element method
program bem
  !宣言部
  integer N
  double precision R
  double precision, PARAMETER :: pi = 3.14159265359
	double precision,allocatable :: end_point(:,:)
	double precision,allocatable :: ans(:)
	double precision,allocatable :: b_vec(:)
  double precision,allocatable :: q_vec(:)
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
  double precision,allocatable :: ans_vec(:),funcU_vec(:)

  !ここからはテストケース
	double precision,allocatable :: test_ans(:)
  double precision,allocatable :: input(:,:)
  double precision input_x,input_y,length,ans_u,hankei
  Integer mesh_num,count
  double precision kai,err_sum,err_ratio
  !実行部
  ! print *, "Please enter integers:R="
  ! read *, R
  R = 100
  ! print *, "Please enter integers:N="
  ! read *, N
  N= 100
  mesh_num = 50


  print *, "N:",N,"Radius:",R,"Mesh:",mesh_num

  allocate( end_point(N+1,2) )
  allocate( P(N,2) )
  allocate( U(N,N) )
  allocate( W(N,N) )
  allocate( U2(mesh_num**2,N) )
  allocate( W2(mesh_num**2,N) )
  allocate( theta1(N,N))
  allocate( theta2(mesh_num**2,N))
  allocate( ans(N) )
  allocate( test_ans(N) )
  allocate( ipiv(N) )
  allocate( b_vec(N) )
  allocate( q_vec(N) )
  allocate( input(mesh_num**2,2) )
  allocate( ans_vec(mesh_num**2) )
  allocate( funcU_vec(N) )
  do i=1, N
    ! 分割端点の設定
    end_point(i,1) = R*cos(2*pi*(i-1)/N)
    end_point(i,2) = R*sin(2*pi*(i-1)/N)
  enddo
    end_point(N+1,1) = end_point(1,1)
    end_point(N+1,2) = end_point(1,2)
  do i =1,N
    ! 分割代表点の設定
    P(i,1) = (end_point(i,1)+end_point(i+1,1))/2
    P(i,2) = (end_point(i,2)+end_point(i+1,2))/2
    ans(i) = funcU( P(i,1), P(i,2) ) !u(Pj)のこと
    test_ans(i) = derFuncU( P(i,1), P(i,2) , R)!微分u(Pj)の理論値的なあれ
    funcU_vec(i) = funcU(end_point(i,1),end_point(i,2))
  enddo

  do i = 1,N
    do j = 1,N
      theta1(i,j) = calcTheta(P(i,1),P(i,2),end_point(j,1),end_point(j,2),end_point(j+1,1),end_point(j+1,2) )
      W(i,j) = theta1(i,j) / (2*pi)
      U(i,j) = calcU( P(i,1),P(i,2),end_point(j,1),end_point(j,2),end_point(j+1,1),end_point(j+1,2),theta1(i,j) )
    end do
  end do

  !Wuの計算
  do i = 1,N
      b_vec(i) = DOT_PRODUCT(W(i,:),ans)
  end do

  Call DGESV(N, 1, U, N, ipiv, b_vec, N, info)

  do i = 1,N
    err_sum = max(err_sum,abs(test_ans(i)-b_vec(i))/test_ans(i))
    err_ratio = max(err_ratio,abs(1-b_vec(i)/test_ans(i)))
  end do

  open(1, file='q_calc.dat', status='replace')
    do i = 1,N
      write (1,*) i, " " , b_vec(i)
      ! write (*,*) i, " " , b_vec(i)
    end do
  close(1)
  open(1, file='q_ans.dat', status='replace')
    do i = 1,N
      write (1,*) i, " " , test_ans(i)
      ! write (*,*) i, " " , test_ans(i)
    end do
  close(1)

  print * , "相対誤差:",err_sum
  print *,"誤差比率最大:", err_ratio

  ! 出力するためのうんコード
  count = 0
  do i = 1,mesh_num**2
  !円形に点を取る場合はこっち
    ! input(i,1) = (R*((i/mesh_num)-1))/(mesh_num+1)*cos(2*pi*(i-1)/mesh_num)
    ! input(i,2) = (R*((i/mesh_num)-1))/(mesh_num+1)*sin(2*pi*(i-1)/mesh_num)
  end do
  ! input = getRoundAxis(mesh_num**2,R)
  !inputを内点としてとる
  do i = 1,mesh_num**2
    input(i,1) = (i/mesh_num)*2*(R)/(mesh_num-1)-R
  end do
  do i = 1,mesh_num**2
    input(i,2) = (mod(i,mesh_num)+1)*2*(R)/(mesh_num-1)-R
  end do

  open(1, file='axis.dat', status='replace')
  do i = 1,mesh_num**2
    write (1,*) input(i,1)," ",input(i,2)
  end do
  close(1)

! 内点計算
  do i = 1,mesh_num**2
    do j = 1,N
      theta2(i,j) = calcTheta(input(i,1),input(i,2),end_point(j,1),end_point(j,2),end_point(j+1,1),end_point(j+1,2) )
      W2(i,j) = theta2(i,j) / (2*pi)
      U2(i,j) = calcU( input(i,1),input(i,2),end_point(j,1),end_point(j,2),end_point(j+1,1),end_point(j+1,2),theta2(i,j) )
    end do
  end do

! ファイル書き込み
open(1, file='data.dat', status='replace')
err_sum=0.0
do i = 1,mesh_num**2
      ans_u = 0.0
      length = sqrt( input(i,1)**2 + input(i,2)**2 )
      if(length < R-1) then
        ! 内点計算
          ans_u = DOT_PRODUCT( U2(i,:),b_vec ) + DOT_PRODUCT( W2(i,:),funcU_vec )
          write (1,*) input(i,1)," ",input(i,2)," ",ans_u
          kai = funcU( input(i,1),input(i,2) )
          err_sum = max(err_sum,abs((kai-ans_u)/kai))
          err_ratio = max(err_ratio,abs(1-ans_u/kai))
          print *, kai, ans_u, ans_u/kai
      end if
end do
close(1)

print *,"相対誤差:", err_sum
print *,"誤差比率最大:", err_ratio

end program bem


function funcU(x,y)
  double precision x,y
  funcU = x**3-3*x*y**2
end function

function derFuncU(x,y,radius)
  double precision x,y,radius
  derFuncU = (3*(x**3) - 9*x*y**2)/radius
end function

function calcTheta(point_x,point_y,x1,y1,x2,y2)
  double precision h,r_1,r_2,ndot1,ndot2,tdot1,tdot2,theta
  double precision point_x,point_y,x1,y1,x2,y2
  double precision p_to_1_vec(2),p_to_2_vec(2),t_vec(2),n_vec(2)

  pi = 3.14159265359

  p_to_1_vec(1) = point_x - x1
  p_to_1_vec(2) = point_y - y1
  p_to_2_vec(1) = point_x - x2
  p_to_2_vec(2) = point_y - y2
  h = sqrt( (x1 - x2)**2 + (y1 - y2)**2 )
  t_vec(1) = (x2-x1)/h
  t_vec(2) = (y2-y1)/h
  n_vec(1) = (y1-y2)/h
  n_vec(2) = (x2-x1)/h

  ndot1 = DOT_PRODUCT (p_to_1_vec, n_vec)
  ndot2 = DOT_PRODUCT (p_to_2_vec, n_vec)
  tdot1 = DOT_PRODUCT (p_to_1_vec, t_vec)
  tdot2 = DOT_PRODUCT (p_to_2_vec, t_vec)

  if( abs(ndot2) < 0.1 .or. abs(ndot1) < 0.1 ) then
    theta = abs(atan2(ndot2,tdot2) - atan2(ndot1,tdot1) )
  else
    theta = atan2(ndot2,tdot2) - atan2(ndot1,tdot1)
  end if
    calcTheta = theta
end function

function calcU(point_x,point_y,x1,y1,x2,y2,theta)
  double precision h,r_1,r_2,ndot1,ndot2,tdot1,tdot2,theta,U,pi
  double precision point_x,point_y,x1,y1,x2,y2
  double precision p_to_1_vec(2),p_to_2_vec(2),t_vec(2),n_vec(2)

  pi = 3.14159265359
  p_to_1_vec(1) = point_x - x1
  p_to_1_vec(2) = point_y - y1
  p_to_2_vec(1) = point_x - x2
  p_to_2_vec(2) = point_y - y2
  h = sqrt( (x1 - x2)**2 + (y1 - y2)**2 )
  t_vec(1) = (x2-x1)/h
  t_vec(2) = (y2-y1)/h
  n_vec(1) = (y1-y2)/h
  n_vec(2) = (x2-x1)/h
  r_1 = sqrt( p_to_1_vec(1)**2 + p_to_1_vec(2)**2 )
  r_2 = sqrt( p_to_2_vec(1)**2 + p_to_2_vec(2)**2 )

  ndot1 = DOT_PRODUCT (p_to_1_vec, n_vec)
  ndot2 = DOT_PRODUCT (p_to_2_vec, n_vec)
  tdot1 = DOT_PRODUCT (p_to_1_vec, t_vec)
  tdot2 = DOT_PRODUCT (p_to_2_vec, t_vec)

  U = (tdot2*dlog(r_2) - tdot1*dlog(r_1) + h - (ndot1)*theta) / (2*pi)
  calcU = U
end function

! function getRoundAxis(mesh_num, radius)
!   Integer mesh_num
!   double precision radius
!   double precision input(mesh_num**2,2)
!   ! double precision getRoundAxis
!   !inputを内点としてとる
!   do i = 1,mesh_num**2
!     input(i,1) = (i/mesh_num)*2*R/(mesh_num-1)-R
!   end do
!   do i = 1,mesh_num**2
!     input(i,2) = (mod(i,mesh_num)+1)*2*R/(mesh_num-1)-R
!   end do
!
!   ! getRoundAxis(:,:) = input(:,:)
!
! end function
