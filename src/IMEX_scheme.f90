! Based on

MODULE IMEX_scheme
USE decimal
USE FVTypes
USE scheme_utils
USE eno_weno

IMPLICIT NONE
  !Algorithm Type
  ! Pares Scheme using original variables
  type, public, extends(FVDiff1DAlgorithm) :: LIMEX1DAlgorithm
    real(kind = dp) :: alpha = 0.0_dp
  contains
    procedure :: update_dt => update_dt_imex
    procedure :: update => update_imex
    procedure :: getStartMessage => start_message_imex
  end type LIMEX1DAlgorithm

CONTAINS
  function start_message_imex(alg) result(message)
      CLASS(LIMEX1DAlgorithm)  :: alg
      CHARACTER(LEN=32)             :: message
      message = "Starting Lin IMEX scheme"
  end function start_message_imex
  
  ! Update time based on CFL condition
  function update_dt_imex (alg, u, CFL) result(dt)
    CLASS(LIMEX1DAlgorithm)  :: alg
    REAL(kind = dp), intent(in)  :: u(:,:), CFL
    REAL(kind = dp)              :: dt,dx
    alg%alpha = maxrho(u, alg%problem) !See scheme utils
    dx = alg%problem%mesh%dx
    dt = CFL*dx/alg%alpha
  end function update_dt_imex

  !!!!!!!!!!! Main methods, used to update solution in each time integration step
    subroutine update_imex(alg, rhs, uold, dt)
      CLASS(LIMEX1DAlgorithm)  :: alg
       real(kind = dp), intent(in) :: uold(:,:)
      real(kind = dp)             :: rhs(:,:), dt
      REAL(kind = dp)               :: dx, alpha
      REAL(kind=dp)             :: crj(4,3)     !WENO5 weigths
      INTEGER                       :: N, j,M, boundary, order, k, i
      REAL(kind = dp), ALLOCATABLE  :: uu(:,:), hh(:,:)
      REAL(kind = dp), ALLOCATABLE  :: fminus(:,:), fplus(:,:), Du(:,:)

      !==================
      N = alg%problem%mesh%N
      M = alg%problem%M
      dx = alg%problem%mesh%dx
      boundary = alg%problem%mesh%bdtype
      alpha = alg%alpha
      order = 3; k = 2  !order = k + 1

      ALLOCATE(uu(0:N+1,M))
      ALLOCATE(hh(N+1,M))
      uu = 0.0_dp; hh = 0.0_dp; crj = 0.0_dp

      ! Compute Weno coefficients
      crj = unif_crj(order)

      ! Add ghost cells
      uu(1:N,:) = uold
      if (boundary == PERIODIC) then
        uu(0,:) = uold(N,:); uu(N+1,:)=uold(1,:)
      else
        uu(0,:) = uold(1,:); uu(N+1,:)=uold(N,:)
      end if

      ! Global Lax Frierichs splitting
      ALLOCATE(fminus(-k:N+k+1,M), fplus(-k:N+k+1,M))
      fminus = 0.0_dp; fplus = 0.0_dp
      !$omp parallel
      !$omp do
      do j = 1,N
        fminus(j,:) = 0.5_dp*(alg%problem%f(uu(j,:))-alpha*uu(j,:))
        fplus(j,:) = 0.5_dp*(alg%problem%f(uu(j,:))+alpha*uu(j,:))
      end do
      !$omp end do
      !$omp end parallel
      do j = 1,(k+1)
        if (boundary == ZERO_FLUX) then
          fminus(1-j,:) = fminus(1,:); fminus(N+j,:) = fminus(N,:)
          fplus(1-j,:) = fplus(1,:); fplus(N+j,:) = fplus(N,:)
        else !PERIODIC
          fminus(1-j,:) = fminus(N-j+1,:); fminus(N+j,:) = fminus(j,:)
          fplus(1-j,:) = fplus(N-j+1,:); fplus(N+j,:) = fplus(j,:)
        end if
      end do

      ! WENO 5 reconstruction
      !$omp parallel
      !$omp do
      do j = 0,N
        do i = 1,M
          ! Diffusion
          hh(j+1,i) = sum(WENO_pm_rec(fminus(j-k+1:j+k+1,i),fplus(j-k:j+k,i),order*2-1, crj))
        end do
      end do
      !$omp end do
      !$omp end parallel

      !Compute Numeric Flux + Diffusion term
      if (boundary == ZERO_FLUX) then
        hh(1,:)=0.0_dp
        hh(N+1,:)=0.0_dp
      end if
      do j = 1,N
        rhs(j,:) = - 1/dx * (hh(j+1,:)-hh(j,:))
      end do
      DEALLOCATE(uu, hh, fminus, fplus)
  end subroutine update_imex
END MODULE IMEX_scheme
