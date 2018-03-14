# Lin IMEX Methods for degenerate convection-diffusion systems of the form:

<a href="https://www.codecogs.com/eqnedit.php?latex=u_{t}&plus;f(u)_{x}&=(B(u)u_{x})_{x},\qquad\forall(x,t)\in\mathbb{R}^{n}\times\mathbb{R}_{&plus;}\\u(x,0)&=u_{0}(x),\qquad\forall&space;x\in\mathbb{R}^{n}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?u_{t}&plus;f(u)_{x}&=(B(u)u_{x})_{x},\qquad\forall(x,t)\in\mathbb{R}^{n}\times\mathbb{R}_{&plus;}\\u(x,0)&=u_{0}(x),\qquad\forall&space;x\in\mathbb{R}^{n}" title="u_{t}+f(u)_{x}&=(B(u)u_{x})_{x},\qquad\forall(x,t)\in\mathbb{R}^{n}\times\mathbb{R}_{+}\\u(x,0)&=u_{0}(x),\qquad\forall x\in\mathbb{R}^{n}" /></a>

The time integration of the semi discrete form is performed with Runge Kutta methods.

## Features
### Mesh: 
At the momento only Cartesian 1D uniform mesh available, using `FVMesh(N,a,b,boundary)` command. Where

`N` = Number of cells

`xinit,xend` = start and end coordinates.

`bdtype` = boundary type (ZERO_FLUX, PERIODIC)

* Problem types: System of Conservation Laws with diffusion term (`CLS1DDiffusionProblem`).

### Algorithms

#### Linearly implicit IMEX Runge-Kutta schemes

(`LIMEX1DAlgorithm`)

(Flux reconstruction uses Comp WENO5)

For IMEX Scheme available RK methods are: H-CN(2,2,2) `H_CN_222`, H-DIRK2(2,2,2) `H_DIRK2_222`, H-LDIRK2(2,2,2) `H_LDIRK2_222`, H-LDIRK3(2,2,2) `H_LDIRK3_222`, SSP-LDIRK(3,3,2) `SSP_LDIRK_332`.

* S. Boscarino, R. Bürger, P. Mulet, G. Russo, L. Villada, *Linearly implicit IMEX Runge Kutta methods for a class of degenerate convection diffusion problems*, SIAM J. Sci. Comput., 37(2), B305–B331

* S. Boscarino, P.G. LeFloch and G. Russo. *High order asymptotic-preserving methods for fully nonlinear relaxation problems*. SIAM J. Sci. Comput., 36 (2014), A377–A395.

* S. Boscarino, F. Filbet and G. Russo. *High order semi-implicit schemes for time dependent partial differential equations*. SIAM J. Sci. Comput. September 2016, Volume 68, Issue 3, pp 975–1001


## Example

```fortran
 ! Initialize variables
  Tend = 0.2_dp      ! Final Time
  CFL = 0.25_dp
  L = 10.0
  M = 4
  bdtype = PERIODIC
  !Run numerical schemes
  N = 100          ! Number of nodes
  CALL setup_problem(0.0_dp, L, N, M, mesh, uinit, bdtype)
  CALL prob%Initialize(mesh, uinit, M, Tend, Flux, JacF, BB)
  ! Save initial data
  name = 'test_1_ini'
  ALLOCATE(results(N, M+1), names(M+1))
  names = ['x       ', 'y1      ','y2      ', 'y3      ', 'y4      ']
  results(:,1) = mesh%x
  results(:,2:5) = uinit(:,1:4)
  CALL save_matrix(results, names, name, 0)
  
  ! Compute solution
  CALL solve(prob, SSP_LDIRK_332, LIMEX_Alg, CFL)

  ! Save solution
  results(:,2:5) = prob%uu(:,1:4)
  name = 'test_1'
  CALL save_matrix(results, names, name, 0)
```

# Disclamer
** Modules developed for personal use, some of them have not been tested enough !!!**
