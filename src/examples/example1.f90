! Test Problem 1
! Vehicular traffic problem
! Test based on example 1 of:
! Burger, Mulet, Villada, A difussion Corrected Multiclass LWR Traffic model
! with Anticipation Lengths and Reaction times, Advances in Applied Mathematics
! and mechanics, 2013
MODULE example1
USE decimal
USE rwdata
USE FVTypes
USE IMEX_scheme
USE FV_Solve
IMPLICIT NONE
PUBLIC example1_run
PRIVATE
INTEGER, PARAMETER :: nvars = 4   !Number of variables
REAL(kind = dp), PARAMETER :: phic = exp(-7.0_dp/exp(1.0_dp))
REAL(kind = dp), PARAMETER :: Vmax(nvars) = [60.0,55.0,50.0,45.0]
REAL(kind = dp), PARAMETER :: Lmin = 0.03
CONTAINS
subroutine example1_run()
  REAL(kind = dp)           :: Tend
  INTEGER                   :: N, M, bdtype
  REAL(kind = dp)           :: dx
  REAL(kind = dp)           :: CFL, L
  REAL(kind = dp), ALLOCATABLE    :: xx(:)
  REAL(kind = dp), ALLOCATABLE    :: uinit(:,:)
  INTEGER                         :: i
  CHARACTER(LEN=32)               :: name           ! File name to save plot data
  REAL(kind = dp), ALLOCATABLE    :: results(:,:)
  CHARACTER(LEN=8), ALLOCATABLE  :: names(:)
  type(Uniform1DMesh) :: mesh
  type(CLS1DDiffusionProblem) :: prob
  type(LIMEX1DAlgorithm)        :: LIMEX_Alg

  ! Zero variables
  Tend = 0.0_dp; dx = 0.0_dp; CFL = 0.0_dp
  ! Initialize variables
  Tend = 0.1_dp      ! Final Time
  CFL = 0.25_dp
  L = 10.0
  bdtype = PERIODIC
  !Run numerical schemes
  N = 100          ! Number of nodes
  M = nvars
  CALL setup_mesh(0.0_dp, L, N, M, mesh, uinit, bdtype)
  CALL prob%Initialize(mesh, uinit, M, Tend, Flux, JacF, BB)
  ! Save initial data
  name = 'test_1_ini'
  ALLOCATE(results(N, M+1), names(M+1))
  names = ['x       ', 'y1      ','y2      ', 'y3      ', 'y4      ']
  results(:,1) = mesh%x
  results(:,2:5) = uinit(:,1:nvars)
  CALL save_matrix(results, names, name, 0)
  
  ! Compute solution
  CALL solve(prob, H_CN_222, LIMEX_Alg, CFL)

  ! Save solution
  results(:,2:5) = prob%uu(:,1:nvars)
  name = 'test_1'
  CALL save_matrix(results, names, name, 0)
  
  DEALLOCATE(results, names, uinit)
end subroutine example1_run

SUBROUTINE setup_mesh(xinit, xend, N, M, mesh, uinit, bdtype)
  INTEGER, INTENT(IN)             :: N, M
  REAL(kind = dp)                 :: L, xinit, xend
  REAL(kind = dp), ALLOCATABLE    :: uinit(:,:)
  INTEGER                         :: i, j, k, ll, bdtype
  type(Uniform1DMesh)             :: mesh
  ! Nodes and weigths for Gauss-Lobatto integration
  real(kind = dp) :: nods(5) = [-0.906179845938664,-0.5384693101056831,0.0,&
  0.5384693101056831,0.906179845938664]
  real(kind = dp) :: wghts(5) = [0.23692688505618908, 0.47862867049936647, &
  0.5688888888888889, 0.47862867049936647, 0.23692688505618908]
  real(kind = dp) :: t_nodes(5), a, b, tmp1(5)
  real(kind = dp) :: ff0(M)
  t_nodes = 0.0_dp; tmp1 = 0.0_dp
  
  ff0 = 0.0_dp
    ! Allocate memory
  call mesh%Initialize(N, xinit, xend, bdtype)
  ALLOCATE(uinit(N,M))
  uinit = 0.0_dp;
  ! Initial Conditions
  DO j = 1, N
    a = mesh%x(j) - mesh%dx/2.0_dp
    b = mesh%x(j) + mesh%dx/2.0_dp
    ! affine transformation from [-1,1] -> [a,b]
    t_nodes = 0.5_dp*(b-a)*nods + 0.5_dp*(b+a)
    do k = 1, M
      do ll = 1, 5
        ff0 = f0(t_nodes(ll))
        tmp1(ll) = ff0(k)
      end do
      ! Compute average integral
      uinit(j,k) = 0.5_dp*(b-a)*dot_product(tmp1, wghts)/mesh%dx
    end do
  END DO
END SUBROUTINE

function f0(x)
  real(kind = dp), intent(in) :: x
  real(kind = dp) :: f0(nvars)
  f0 = 0.5_dp*exp(-(x-3)**2)*[0.2,0.3,0.2,0.3]
end function f0  

FUNCTION flux(phi) RESULT(ff)
  REAL(kind = dp), INTENT(IN)  :: phi(:)
  REAL(kind = dp)              :: ff(SIZE(phi,1))
  ff = 0.0_dp
  ff = VV(sum(phi))*phi*Vmax
END FUNCTION flux

FUNCTION JacF(phi) RESULT(F)
  REAL(kind = dp), INTENT(IN)  :: phi(:)
  REAL(kind = dp)              :: F(SIZE(phi,1),SIZE(phi,1))
  INTEGER                      :: M, i, j
  REAL(kind = dp)              :: Vphi, VPphi
  F = 0.0_dp
  M = SIZE(phi,1)
  Vphi = VV(sum(phi))
  VPphi = VP(sum(phi))
  do i =  1,M
    do j = 1,M
      if (i == j) then
        F(i,j)=Vmax(i)*(Vphi + phi(i)*VPphi)
      else
        F(i,j)=Vmax(i)*(phi(i)*VPphi)
      end if
    end do
  end do
END FUNCTION JacF

! Beta function for diffusion
function Beta(u)
  real(kind = dp), intent(in) :: u
  real(kind = dp) :: Beta
  INTEGER         :: i
  Beta = 0.0_dp
  if (u > phic ) then
    Beta = -VP(u)*Lmin*u/nvars*sum(Vmax)/nvars
  end if
end function Beta

! Difussion Matrix
function BB(phi) result(k)
  real(kind = dp), intent(in) :: phi(:)
  real(kind = dp) :: k(SIZE(phi,1), SIZE(phi,1)), tmp
  INTEGER         :: M, i
  k = 0.0_dp
  M = size(phi,1)
  tmp = Beta(sum(phi))
  do i = 1,M
     k(i,i) = tmp
  end do
end function BB

function VV(phi)
  real(kind = dp), intent(in) :: phi
  REAL(kind = dp) :: VV
  VV = 1.0_dp
  if (phi > phic) then
    VV = 1.0 - phi
  end if
end function VV

function VP(phi)
  real(kind = dp), intent(in) :: phi
  REAL(kind = dp) :: VP
  VP = 0.0_dp
  if (phi > phic) then
    VP = -1.0
  end if
end function VP

END MODULE example1