MODULE FVTypes
use decimal
IMPLICIT NONE
  ENUM, bind(c)
        ENUMERATOR :: H_CN_222
        ENUMERATOR :: H_DIRK2_222
        ENUMERATOR :: H_LDIRK2_222
        ENUMERATOR :: H_LDIRK3_222
        ENUMERATOR :: SSP_LDIRK_332
  END ENUM
  ENUM, bind(c)
        ENUMERATOR :: ZERO_FLUX
        ENUMERATOR :: PERIODIC
        ENUMERATOR :: DIRICHLET
  END ENUM
    abstract interface
      function AbstractFlux(phi) result(f)
          import  :: dp
          real(kind = dp), intent(in) :: phi(:)
          real(kind = dp) :: f(SIZE(phi,1))
      end function AbstractFlux
      function AbstractJacF(phi) result(Jf)
          import  :: dp
          real(kind = dp), intent(in) :: phi(:)
          real(kind = dp) :: Jf(SIZE(phi,1),SIZE(phi,1))
      end function AbstractJacF
      function AbstractDiffMat(phi) result(k)
          import :: dp
          real(kind = dp), intent(in) :: phi(:)
          real(kind = dp) :: k(SIZE(phi,1), SIZE(phi,1))
      end function AbstractDiffMat
    end interface
!Mesh types
  type, public :: Uniform1DMesh
      integer        :: N
      REAL(KIND=dp), ALLOCATABLE  :: x(:)
      real(kind=dp)  :: dx
      integer        :: bdtype
      CONTAINS
      procedure :: initialize => init_uniform1d
  end type Uniform1DMesh
!Problem Type
  type, abstract, public :: AbstractCLS1DProblem
      type(Uniform1DMesh)     :: mesh           !Mesh
      integer                 :: M              !Number of equations
      REAL(KIND=dp), ALLOCATABLE  :: u0(:,:)    !Initial condition
      REAL(KIND=dp), ALLOCATABLE  :: uu(:,:)    !Current Solution
      real(kind=dp)           :: Tend           !End time
  end type AbstractCLS1DProblem

  type, public, extends(AbstractCLS1DProblem) :: CLS1DDiffusionProblem
      procedure(Abstractflux),pointer,nopass :: f
      procedure(AbstractJacF),pointer,nopass :: Jf
      procedure(AbstractDiffMat),pointer,nopass :: K
      CONTAINS
      procedure :: initialize => init_1DDiffProblem
  end type CLS1DDiffusionProblem

!Algorithms type 
  type, abstract, public :: AbstractFV1DAlgorithm
  CONTAINS
    procedure(abstractUpdateRHS), deferred :: update
    procedure(abstractCDT),deferred :: update_dt
    procedure(abstractGetStartMessage), deferred :: getStartMessage
  end type AbstractFV1DAlgorithm

  type, abstract, public, extends(AbstractFV1DAlgorithm) :: FVDiff1DAlgorithm
    TYPE(CLS1DDiffusionProblem) :: problem
  end type FVDiff1DAlgorithm

  abstract interface
      function abstractCDT (alg, u, CFL) result(dt)
        import :: dp
        import :: AbstractFV1DAlgorithm
        CLASS(AbstractFV1DAlgorithm) :: alg
        REAL(kind = dp), intent(in)  :: u(:,:), CFL
        REAL(kind = dp)              :: dt
      end function abstractCDT
      subroutine abstractUpdateRHS(alg, rhs, uold, dt)
          import :: dp
          import AbstractFV1DAlgorithm
          CLASS(AbstractFV1DAlgorithm)  :: alg
          real(kind = dp), intent(in) :: uold(:,:)
          real(kind = dp) :: rhs(:,:), dt
      end subroutine abstractUpdateRHS
      function abstractGetStartMessage(alg) result(message)
            import AbstractFV1DAlgorithm
            CLASS(AbstractFV1DAlgorithm)  :: alg
            CHARACTER(LEN=32)             :: message
      end function abstractGetStartMessage
end interface

  CONTAINS
  !Mesh types
  subroutine init_uniform1d(mesh, N, xinit, xend, bdtype)
      CLASS(Uniform1DMesh)    :: mesh
      integer           :: N, bdtype, i
      real(kind=dp)     :: xinit, xend, dx, L
      REAL(KIND=dp), ALLOCATABLE  :: xx(:)
      !========================
      L = xend-xinit
      dx = L/N
      ALLOCATE(mesh%x(N))
      mesh%x = 0.0_dp;
      mesh%x = (/(i*dx+dx/2+xinit,i=0, N-1)/)   
      mesh%N = N;mesh%dx = dx; mesh%bdtype = bdtype
  end subroutine init_uniform1d

  ! Problem type
  subroutine init_1DDiffProblem(problem, mesh, u0, M, Tend, Flux, JacF, DiffMat)
      CLASS(CLS1DDiffusionProblem)  :: problem
      type(Uniform1DMesh)           :: mesh
      REAL(KIND=dp)                 :: u0(:,:), Tend
      integer           :: M,N
      procedure(Abstractflux) :: Flux
      procedure(AbstractJacF) :: JacF
      procedure(AbstractDiffMat) :: DiffMat
      !========================
      N = mesh%N
      ALLOCATE(problem%uu(N,M), problem%u0(N,M))
      problem%M = M
      problem%uu = u0
      problem%u0 = u0
      problem%Tend = Tend
      problem%mesh = mesh
      problem%f => Flux
      problem%Jf => JacF
      problem%K => DiffMat
  end subroutine init_1DDiffProblem

END MODULE FVTypes