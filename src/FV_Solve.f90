MODULE FV_Solve
USE decimal
USE FVTypes
USE scheme_utils
use omp_lib
USE LinearSystems
USE morse
use util

IMPLICIT NONE

   interface solve
      module procedure solve_1D_Diff
   end interface solve

  type, public :: RKTable
      integer        :: order
      REAL(KIND=dp), ALLOCATABLE  :: ct(:), At(:,:), bt(:)
      REAL(KIND=dp), ALLOCATABLE  :: c(:), A(:,:), b(:)
  end type RKTable

CONTAINS

subroutine solve_1D_diff(problem, time_scheme, algorithm, CFL)
  INTEGER                   :: time_scheme
  REAL(kind = dp)           :: CFL, Tend
  INTEGER                   :: N, M
  TYPE(CLS1DDiffusionProblem)        :: problem
  CLASS(FVDiff1DAlgorithm)       :: algorithm
  !Get variables
  Tend = problem%Tend
  N = problem%mesh%N
  M = problem%M
  !Set problem
  algorithm%problem = problem
  CALL generic_time_integration(problem%uu, time_scheme, algorithm, CFL,Tend, N,M,problem,problem%u0)
END subroutine solve_1D_diff

subroutine generic_time_integration(uu, time_scheme, algorithm, CFL, Tend, N,M,problem,u0)
  INTEGER                   :: time_scheme
  INTEGER                   :: N, M, percentage, i, j
  REAL(kind = dp)           :: dx,dt, limit, Tend, CFL
  REAL(kind = dp)           :: tt, uu(:,:), u0(:,:)
  REAL(kind=dp)             :: tiempo1, tiempo2
  CLASS(FVDiff1DAlgorithm)       :: algorithm
  TYPE(CLS1DDiffusionProblem)        :: problem
  CHARACTER(LEN=32)         :: message
  type(RKTable)             :: rktab
  REAL(kind=dp), allocatable :: ki(:), kjs(:), kjh(:), uold(:), rhs(:,:)
  REAL(kind=dp), allocatable :: kj(:,:), Cu(:,:), us(:), uh(:)
  type(Sparse) A, BB

  !Allocate variables
  ALLOCATE(ki(N*M), kjs(N*M), kjh(N*M), uold(N*M), rhs(N,M), Cu(N,M))
  ALLOCATE(us(N*M), uh(N*M))
  ALLOCATE(Kj(N*M,rktab%order))
  ! Zero variables
  tiempo1 = 0.0_dp;tiempo2 = 0.0_dp
  dx = problem%mesh%dx;Kj = 0.0_dp
  ! Set progress variables
  percentage = 0
  limit = Tend/20.0_dp

  ! Print Start Message
  message = algorithm%GetStartMessage()
  Print *, message
  ! Start timing
  tiempo1 = omp_get_wtime( )
  ! get tableau
  CALL get_time_tableau(rktab, time_scheme)
  ! Start loop
  DO WHILE (tt < Tend)
    uold = reshape(transpose(uu), shape(ki))
    dt = algorithm%update_dt(uu, CFL)
    if (tt + dt > Tend) THEN
      dt = Tend - tt
    end if
    !IMEX algorithm
    Ki = 0.0_dp
    ! i step
    do i = 1,RKTab%order
      Kjs = 0.0_dp; Kjh = 0.0_dp
      do j = 1,(i-1)
        Kjs = Kjs + RKTab%At(i,j)*Kj(:,j)
        Kjh = Kjh + RKTab%A(i,j)*Kj(:,j)
      end do
      us = uold + dt*Kjs
      uh = uold + dt*Kjh
      CALL assamble_B(BB, us,N,M,algorithm, problem)
      CALL assemble_A(A, BB,N,M,dt,dx,RKTab%A(i,i))
      !Reconstruct flux with comp weno5 see: WENO_Scheme.jl
      rhs = 0.0_dp
      CALL algorithm%update(rhs, transpose(reshape(us,[M,N])), dt)
      Cu = rhs(2:N+1,:)-rhs(1:N,:)
      Ki = -(1.0_dp/dx)*reshape(transpose(Cu),shape(Ki)) + (1.0_dp/dx**2) * BB%dot(uh)
      !Solve linear system
      CALL solve_system(A, Ki, 1)
      Kj(:,i) = Ki
    end do
    do j = 1,RKTab%order
      uold = uold + dt*RKTab%b(j)*Kj(:,j)
    end do
    uu(:,:) = transpose(reshape(uold,[M,N]))

    !Print progress
    IF (tt > limit) THEN
      percentage = percentage + 5
      limit = limit + Tend/20.0_dp
      print *, percentage, "% completed"
    END IF
    tt = tt + dt
  END DO
  print *, "completed..."

  ! End timing and print total time
  tiempo2 = omp_get_wtime( )
  PRINT*,'tiempo de CPU = ',tiempo2-tiempo1
  DEALLOCATE(kj, ki, kjs, kjh, uold, rhs, Cu, us, uh)
end subroutine generic_time_integration

! ASSEMBLE CRS Matrices
subroutine assamble_B(BB, u,N,M,alg, prob)
  REAL(kind=dp)   :: u(:)
  integer         :: N,M,boundary, nnz, cent, i, j, ir, ic
  integer, allocatable :: idr(:), idc(:)
  TYPE(CLS1DDiffusionProblem)        :: prob
  CLASS(FVDiff1DAlgorithm)       :: alg
  REAL(kind=dp), allocatable   :: uleft(:), uright(:), ul(:), uc(:), ur(:), vals(:), tmp(:,:)
  type(Sparse) BB
  logical   :: do_update

  ! Allocate variables
  ALLOCATE(uleft(M), uright(M),ul(M),ur(M),uc(M),tmp(M,M))
  nnz=M*M*(N-2)*3+M*M*2*2
  ALLOCATE(idr(nnz), idc(nnz), vals(nnz))

  !Get some data
  boundary = alg%problem%mesh%bdtype
  ! TODO: Dirichlet missing
  if (boundary == PERIODIC) then
    uleft = u(((N-1)*M+1):(N*M))    
    uright = u(1:M)
  else
    uleft = u(1:M)    !TODO: default is ZERO_flux?
    uright = u(((N-1)*M+1):(N*M))
  end if
  idr = 0; idc = 0
  vals = 0.0_dp
  cent = 1
  do i = 1,N
    do j = 1,N
      do_update = .false.
      if (i == j) then
        ul = merge(u(((i-2)*M+1):((i-1)*M)),uleft,i>1)
        uc = u(((i-1)*M+1):(i*M))
        ur = merge(u((i*M+1):((i+1)*M)),uright,i < N)
        tmp = -0.5_dp*(prob%K(ul)+2*prob%K(uc)+prob%K(ur))
        do_update = .true.
      else if (j == i+1) then
        uc=u(((i-1)*M+1):(i*M))
        ur=merge(u((i*M+1):((i+1)*M)), uright, i < N)
        tmp = 0.5_dp*(prob%K(uc)+prob%K(ur))
        do_update = .true.
      else if (j == i-1) then
        ul=merge(u(((i-2)*M+1):((i-1)*M)),uleft,i>1)
        uc=u(((i-1)*M+1):(i*M))
        tmp = 0.5*(prob%K(ul)+prob%K(uc))
        do_update = .true.
      end if
      if (do_update) then
        do ir = 1,M
          do ic = 1,M
            vals(cent) = tmp(ir,ic); idr(cent) = ((i-1)*M+ir); idc(cent) = ((j-1)*M+ic)
            cent = cent + 1
          end do
        end do
      end if
    end do
  end do
  BB = Sparse(idr,idc,vals)
  DEALLOCATE(uleft, uright, idr, idc, vals,ul,ur,uc,tmp)
end subroutine assamble_B

subroutine assemble_A(A,BB,N,M,dt,dx,rktabAii)
  INTEGER M,N,k,kk
  real(kind=dp)   :: dt, dx, rktabAii
  type(Sparse) A, BB
  A = BB
  ! first compute -dt/dx**2*RKTab%A(i,i)*BB
  A%aa = -dt/dx**2*rktabAii*A%aa
  ! now A = I-dt/dx**2*RKTab%A(i,i)*BB
  ! We need find the diagonal and then add 1
  ! since BB is tridiagonal we assume all diagonal elements are non-zero
  A%aa(1) = A%aa(1) + 1.0_dp
  do k = 2,M*N
    do kk = A%row(k),A%nzero
      if (A%column(kk) == k) then
        EXIT
      end if
    end do
    A%aa(kk) = A%aa(kk) + 1.0_dp 
  end do
end subroutine assemble_A

! Runge Kutta Tableaux for IMEX methods
subroutine get_time_tableau(rktab, time_scheme)
  type(RKTable)             :: rktab
  INTEGER                   :: time_scheme
  REAL(kind=dp)             :: gamma
  ! Get Time Integration tableau
  if (time_scheme == H_CN_222) then
    allocate(rktab%ct(2), rktab%At(2,2), rktab%bt(2))
    allocate(rktab%c(2), rktab%A(2,2), rktab%b(2))
    rktab%order = 2
    rktab%ct = [0.0_dp,1.0_dp]
    rktab%At = reshape([0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], shape(rktab%At))
    rktab%bt = [0.5_dp,0.5_dp]
    rktab%c = [0.0_dp,1.0_dp]
    rktab%A = reshape([0.0_dp, 0.5_dp, 0.0_dp, 0.5_dp], shape(rktab%A))
    rktab%b = [0.5_dp,0.5_dp]
  else if (time_scheme == H_DIRK2_222) then
    allocate(rktab%ct(2), rktab%At(2,2), rktab%bt(2))
    allocate(rktab%c(2), rktab%A(2,2), rktab%b(2))
    rktab%order = 2
    rktab%ct = [0.0_dp,1.0_dp]
    rktab%At = reshape([0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], shape(rktab%At))
    rktab%bt = [0.5_dp,0.5_dp]
    rktab%c = [0.5_dp,0.5_dp]
    rktab%A = reshape([0.5_dp, 0.0_dp, 0.0_dp, 0.5_dp], shape(rktab%A))
    rktab%b = [0.5_dp,0.5_dp]
  else if (time_scheme == H_LDIRK2_222) then
    gamma = 1.0_dp-1.0_dp/sqrt(2.0_dp)
    allocate(rktab%ct(2), rktab%At(2,2), rktab%bt(2))
    allocate(rktab%c(2), rktab%A(2,2), rktab%b(2))
    rktab%order = 2
    rktab%ct = [0.0_dp,1.0_dp]
    rktab%At = reshape([0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], shape(rktab%At))
    rktab%bt = [0.5_dp,0.5_dp]
    rktab%c = [gamma,1.0_dp - gamma]
    rktab%A = reshape([gamma, 1.0_dp-2.0_dp*gamma, 0.0_dp, gamma], shape(rktab%A))
    rktab%b = [0.5_dp,0.5_dp]
  else if (time_scheme == H_LDIRK3_222) then
    gamma = (3.0_dp+sqrt(3.0_dp))/6.0_dp
    allocate(rktab%ct(2), rktab%At(2,2), rktab%bt(2))
    allocate(rktab%c(2), rktab%A(2,2), rktab%b(2))
    rktab%order = 2
    rktab%ct = [0.0_dp,1.0_dp]
    rktab%At = reshape([0.0_dp, 1.0_dp, 0.0_dp, 0.0_dp], shape(rktab%At))
    rktab%bt = [0.5_dp,0.5_dp]
    rktab%c = [gamma,1.0_dp - gamma]
    rktab%A = reshape([gamma, 1.0_dp-2.0_dp*gamma, 0.0_dp, gamma], shape(rktab%A))
    rktab%b = [0.5_dp,0.5_dp]
  else if (time_scheme == SSP_LDIRK_332) then
    allocate(rktab%ct(3), rktab%At(3,3), rktab%bt(3))
    allocate(rktab%c(3), rktab%A(3,3), rktab%b(3))
    rktab%order = 3
    rktab%ct = [0.0_dp,0.5_dp,1.0_dp]
    rktab%At = reshape([0.0_dp, 0.5_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, 0.0_dp], shape(rktab%At))
    rktab%bt = [1.0/3.0_dp,1.0/3.0_dp,1.0/3.0_dp]
    rktab%c = [0.25_dp,0.25_dp,1.0_dp]
    rktab%A = reshape([0.25_dp,0.0_dp,1.0/3.0_dp,0.0_dp,0.25_dp,1.0/3.0_dp,0.0_dp,0.0_dp,1.0/3.0_dp], shape(rktab%A))
    rktab%b = [1.0/3.0_dp,1.0/3.0_dp,1.0/3.0_dp]
  else
    STOP 'time scheme not available...'
  end if
end subroutine get_time_tableau
END MODULE FV_Solve