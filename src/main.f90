PROGRAM LIMEX
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Finite Volume Methods: Linear IMEX
  !! Linearly implicit IMEX Runge-Kutta schemes
  !! Based on:
  !! S. Boscarino, R. Bürger, P. Mulet, G. Russo, L. Villada, Linearly implicit IMEX
  !! Runge Kutta methods for a class of degenerate convection difussion problems                               !!
  !!                                                      
  !!                                                      
  !! Author: Paul Mendez                                  
  !! Date: 2/June/2017                                    
  !!                                                      
  !! Version: 0.4                                           
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  USE example1
  IMPLICIT NONE

  ! Init threads
  !integer t
  !call omp_set_num_threads( 4 )
  !t = omp_get_max_threads()
  !write(*,*)'Number of threads:',t

  ! Run examples
  CALL example1_run()
END PROGRAM LIMEX