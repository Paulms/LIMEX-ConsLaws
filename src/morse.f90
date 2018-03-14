MODULE morse
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Save matrices on Sparse CRS format                   !!
  !!                                                      !!
  !! PROCEDURES:                                          !!
  !!                                                      !!
  !!   init_matrix    = Construct CRS Matrix              !!
  !!   print_flat    = print internal structure of Sparse !!
  !!                    matrix                            !!
  !!  search_element = search element by (i,j) indexes    !!
  !!                      A(i,j)                          !!
  !!  vec_prod = computes matrix vector product           !!
  !!                                                      !!
  !! Author: Paul Mendez                                  !!
  !! e-mail: paul.mendez@udec.cl                          !!
  !! Date: 11/Septiembre/2016                             !!
  !!                                                      !!
  !! Version: 0.5                                         !!
  !! last revision: 1/Noviembre/2016                      !!   
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  USE decimal
  IMPLICIT NONE
  PRIVATE
  PUBLIC  :: Sparse, print_flat, print_matrix

  TYPE Sparse
    ! Matrices in CRS format
    INTEGER                     :: nn, nzero
    REAL (kind=dp), ALLOCATABLE :: aa(:)              ! values
    INTEGER, ALLOCATABLE        :: row(:), column(:)  ! row columns
  CONTAINS
    PROCEDURE :: ix => search_element
    PROCEDURE :: dot => vec_prod
  END TYPE Sparse

  INTERFACE Sparse
    MODULE PROCEDURE init_matrix
  END INTERFACE

CONTAINS
  FUNCTION init_matrix(rows, columns, values)
    !======================================================
    ! Construct a CRS Matrix
    ! Variables:
    !   rows, columns, values.
    !======================================================
    TYPE(Sparse) init_matrix
    INTEGER       :: rows(:), columns(:)
    REAL(kind=dp) :: values(:)
    INTEGER       :: nr, N, element, j,jj, row, nzero
    INTEGER, allocatable :: ccols(:)
    REAL(kind=dp), allocatable  :: cvals(:)
    integer, allocatable :: colidx(:), idx
    !Check sizes
    if (size(rows,1) /= size(columns,1) .or. size(rows,1) /= size(values,1)) THEN
      print *, "Error colums, rows and values have different sizes: ", size(rows,1), size(columns,1), size(values,1)
      STOP
    end if
    nzero = size(values,1)
    nr = maxval(rows)
    ALLOCATE(init_matrix%row(nr+1),init_matrix%column(nzero), init_matrix%aa(nzero))
    init_matrix%row = 0.0_dp; init_matrix%column = 0.0_dp
    init_matrix%aa = 0.0_dp
    init_matrix%nn = nr
    init_matrix%nzero = nzero

    init_matrix%row(1) = 1
    DO row = 1, nr
      element = 0
      DO j = 1,size(rows,1)
        IF (rows(j)==row) then
          element = element + 1 
        end if
      END DO
      ALLOCATE(colidx(element))
      colidx = 0
      idx = 0
      DO j = 1,size(rows,1)
        IF (rows(j)==row) then
          idx = idx + 1
          colidx(idx) = j 
        end if
      END DO
      init_matrix%row(row+1) = init_matrix%row(row) + element
      !sort columns
      ALLOCATE(ccols(element), cvals(element))
      ccols = (/(columns(colidx(jj)), jj=1,element,1)/)
      cvals = (/(values(colidx(jj)), jj=1,element,1)/)
      DEALLOCATE(colidx)
      CALL quicksort(ccols,cvals,1,size(ccols,1))
      init_matrix%column(init_matrix%row(row):init_matrix%row(row+1)-1) = ccols
      init_matrix%aa(init_matrix%row(row):init_matrix%row(row+1)-1) = cvals
      DEALLOCATE(ccols,cvals)
    END DO
  END FUNCTION

  ! Print internal structure of Sparse Type
  SUBROUTINE print_flat(this)
    TYPE(Sparse),  INTENT(in) :: this
    PRINT *, 'aa = ', this%aa
    PRINT *, 'row = ', this%row
    PRINT *, 'column = ', this%column
  END SUBROUTINE print_flat

  ! Search element by index
  FUNCTION search_element(this, row, column) result(elemento)
    CLASS(Sparse), INTENT(in)  :: this               ! Sparse Matrix
    REAL (kind=dp)              :: elemento           ! Element A(i,j)
    INTEGER                     :: row, column, ii    
    ! check indices in bounds
    IF (row < 1 .OR. column < 1 .or. row > this%nn .or. column > this%nn) THEN
      PRINT *, "Error: indices out of bounds", this%nn
      STOP
    ELSE
      elemento = 0 ! Asumimos que el espacio está vacio                                   
      DO ii = this%row(row), (this%row(row+1)-1)
        ! Check if space is not empty
        IF (this%column(ii) == column) THEN
          ! get value       
          elemento = this%aa(ii)
          EXIT                 
        END IF
      END DO
    END IF
  END FUNCTION search_element

  ! Multiply sparse matrix by a vector
  FUNCTION vec_prod(this, vector) RESULT(output)
    CLASS(Sparse), INTENT(in)    :: this           ! Sparse matrix
    REAL (kind=dp), INTENT(in)    :: vector (:)    ! input vector
    REAL (kind=dp)                :: output (this%nn) ! output of product
    INTEGER                       :: ii, jj
    ! Check dimensions
    IF (this%nn /= SIZE(vector)) THEN
      PRINT *, "Error: incompatible dimensions, Matrix ", this%nn, "Vector: ", SIZE(vector)
      STOP
    END IF
    output = 0
    ! Compute product
    DO ii = 1, this%nn
      DO jj = this%row(ii), (this%row(ii+1)-1)
        output(ii) = output(ii) + this%aa(jj)*vector(this%column(jj))
      END DO
    END DO
  END FUNCTION vec_prod

  ! Print full matrix
  SUBROUTINE print_matrix(mat)
    type (Sparse)   :: mat
    INTEGER         :: cols, rows, i, j
    rows = mat%nn
    cols = mat%nn
    do i=1,rows
      write (*,"("//trim(str(cols))//"(F10.4))") ( mat%ix(i,j), j=1,cols )
    end do
  END SUBROUTINE

  CHARACTER(len=20) FUNCTION str(k)
  !   "Convert an integer to string."
    integer, intent(in) :: k
    write (str, *) k
    str = adjustl(str)
  END FUNCTION str

  recursive subroutine quicksort(a, b, first, last)
  integer  a(:), x, t
  real(kind=dp)  b(:), s
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     s = b(i);  b(i) = b(j);  b(j) = s
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, b, first, i-1)
  if (j+1 < last)  call quicksort(a, b, j+1, last)
end subroutine quicksort
END MODULE morse