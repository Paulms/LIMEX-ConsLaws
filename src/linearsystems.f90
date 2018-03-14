MODULE LinearSystems
  !
  ! Metodos iterativos
  !
  USE decimal
  USE morse
  !
  PRIVATE
  !
  PUBLIC solve_system
  !
CONTAINS
  !
  SUBROUTINE solve_system(A,xx,metodo)
    !
    ! meta-subroutina (general)
    !
    TYPE(sparse)              :: A
    INTEGER, INTENT(in)       :: metodo
    REAL(kind=dp)             :: xx(:)
    INTEGER                   :: ierr
    REAL(kind=dp)             :: start, finish  !Vars para medir tiempo
    !
    ierr = 0
    !
    !CALL cpu_time(start)
    SELECT CASE(metodo)
       !
    CASE(1)        !
       !PRINT *,"===== Solving: GMRES (SparceKit), prec: ILUT ====="
       CALL sol_gmres(A,xx)
       !
    CASE default
       PRINT*,'Metdodo aun no implementado ( M: sistemas S:solve_system)!! ,  metodo =', metodo
    END SELECT
    !
    !CALL cpu_time(finish)
    !PRINT '("CPU Time = ",f9.6," seconds.")',finish-start
    DEALLOCATE(A%aa,A%row,A%column,STAT=ierr)
    IF(ierr/=0) THEN
       PRINT*,'problemas liberando la memoria!! (M: Sistemas S: solucion)'
       STOP
    END IF
    !
  END SUBROUTINE solve_system
  !
   SUBROUTINE sol_gmres(A,xx)
    !
    ! Subrutina para usar GMRES de la libreria SPARSKIT2
    !
    TYPE(sparse)               :: A
    REAL(kind=dp)              :: xx(:),bb(size(xx))
    !
    INTEGER                    :: ipar(16),lfil,nwk,nrow,ierr,maxits,its
    INTEGER, ALLOCATABLE       :: iw(:),ju(:), jlu(:)
    REAL(kind=dp), ALLOCATABLE :: wk(:),alu(:)
    REAL(kind=dp)              :: fpar(16), tol,  res
    !
    ipar = 0; fpar = 0.0_dp; tol = 0.0_dp; lfil = 0; nwk = 0; nrow = 0; ierr = 0; maxits = 0
    tol = 0.0_dp; bb = 0.0_dp

    lfil   = 16
    maxits = A%nn
    nwk    = 2*lfil*A%nn+A%nn+1
    !
    bb = xx
    !
    ALLOCATE(iw(A%nn*2),ju(A%nn), alu(nwk), jlu(nwk))
    !
    iw = 0; ju = 0; alu = 0.0_dp; jlu = 0
    !
    ipar(1) = 0
    ipar(2) = 2
    ipar(3) = 1
    ipar(5) = 16
    ipar(4) = (A%nn + 3)*(ipar(5)+2) + ((ipar(5)+1)*ipar(5))/2
    !
    ALLOCATE(wk(ipar(4)))
    wk = 0.0_dp
    !
    ipar(6) = maxits
    fpar(1) = 1.0e-05_dp
    fpar(2) = 1.0e-10_dp
    !
    !   definiendo el precondicionador ILUT
    !    
    tol = 1.0e-05_dp
    !
    CALL ilut(A%nn,A%aa,A%column,A%row,lfil,tol,alu,jlu,ju,nwk,wk,iw,ierr)
    !
    !PRINT*,'ierr de ilut = ', ierr
    ipar(1) = 0
    ipar(2) = 2
    !
    its = 0
    res = 0.0_dp
    !
    main: DO
       !
       CALL gmres(A%nn,bb,xx,ipar,fpar,wk)
       !
       IF(ipar(7) /= its) its = ipar(7)
       !
       res = fpar(5)
       !
       SELECT CASE(ipar(1))
       CASE(1)
          CALL amux(A%nn, wk(ipar(8)), wk(ipar(9)), A%aa, A%column, A%row)
          CYCLE main
       CASE(2)
          CALL atmux(A%nn, wk(ipar(8)), wk(ipar(9)), A%aa, A%column, A%row)
          CYCLE main
       CASE(3,5)
          CALL lusol(A%nn,wk(ipar(8)),wk(ipar(9)),alu,jlu,ju)
          CYCLE main
       CASE(4,6)
          CALL lutsol(A%nn,wk(ipar(8)),wk(ipar(9)),alu,jlu,ju)
          CYCLE main
       CASE(0)
          !PRINT *, 'Se alcanzo la convergencia!!'
          EXIT main
       CASE(-1)
          PRINT *, 'El solver itero mas alla de lo permitido. No hay convergencia.'
          STOP
       CASE(-2)
          PRINT *, 'El solver no tiene suficiente espacio para trabajar'
          PRINT *, '(espacio minimo que se requiere es de ', ipar(4),' elementos)'
          STOP
       CASE(-3)
          PRINT *, 'Se dividio por cero en algun momento. No hay convergencia.'
          STOP
       CASE default
          PRINT *, 'Se termino el solver con el codigo de error =', ipar(1)
          STOP
       END SELECT
       !
    END DO main
    !
    PRINT*,'iterations = ', ipar(7)
    !PRINT*,'residual     = ', fpar(5)
    !
    DEALLOCATE(iw,wk,ju,alu,jlu, STAT=ierr)
    IF(ierr/=0) THEN
       PRINT*,'problemas liberando la memoria!! (M:sistemas S: sol_gmres)'
       STOP
    END IF
    !
  END SUBROUTINE sol_gmres
END MODULE LinearSystems






