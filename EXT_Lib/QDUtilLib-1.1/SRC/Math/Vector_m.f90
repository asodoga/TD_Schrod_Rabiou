!===============================================================================
!===============================================================================
!This file is part of QDUtil.
!
!===============================================================================
! MIT License
!
! Copyright (c) 2022 David Lauvergnat
!
! Permission is hereby granted, free of charge, to any person obtaining a copy
! of this software and associated documentation files (the "Software"), to deal
! in the Software without restriction, including without limitation the rights
! to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is
! furnished to do so, subject to the following conditions:
!
! The above copyright notice and this permission notice shall be included in all
! copies or substantial portions of the Software.
!
! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
! IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
! SOFTWARE.
!===============================================================================
!===============================================================================
MODULE QDUtil_Vector_m
  IMPLICIT NONE

  PRIVATE

  PUBLIC Sort_Vec
  INTERFACE Sort_Vec
    MODULE PROCEDURE QDUtil_Sort_RVec
  END INTERFACE

  PUBLIC :: Test_QDUtil_Vector
  CONTAINS
  !================================================================
  !   Sort a real matrix Rmat (in the same vector)
  !   subroutine
  !================================================================
  SUBROUTINE QDUtil_Sort_RVec(RVec,sort_type)
    USE QDUtil_NumParameters_m
    IMPLICIT NONE

    real(kind=Rkind), intent(inout)              :: RVec(:)
    integer,          intent(in),    optional    :: sort_type

    integer            :: sort_type_loc
    integer, parameter :: sort_type_default = 1 ! ascending sort

    real(kind=Rkind)   :: a
    integer            :: i,j

    IF (present(sort_type)) THEN
      sort_type_loc = sort_type
    ELSE
      sort_type_loc = sort_type_default
    END IF

    SELECT CASE (sort_type_loc)
    CASE (1) ! ascending
      DO i=lbound(RVec,dim=1),ubound(RVec,dim=1)
      DO j=i+1,ubound(RVec,dim=1)
       IF (RVec(i) > RVec(j)) THEN
          ! permutation
          a=RVec(i)
          RVec(i)=RVec(j)
          RVec(j)=a
        END IF
      END DO
      END DO
    CASE (-1) ! descending
      DO i=lbound(RVec,dim=1),ubound(RVec,dim=1)
      DO j=i+1,ubound(RVec,dim=1)
       IF (RVec(i) < RVec(j)) THEN
          ! permutation
          a=RVec(i)
          RVec(i)=RVec(j)
          RVec(j)=a
        END IF
      END DO
      END DO
    END SELECT

  END SUBROUTINE QDUtil_Sort_RVec

  SUBROUTINE Test_QDUtil_Vector()
    USE QDUtil_Test_m
    USE QDUtil_NumParameters_m
    USE QDUtil_RW_MatVec_m
    IMPLICIT NONE

    TYPE (test_t)                    :: test_var
    logical                          :: res_test

    integer                          :: io,ioerr
    real(kind=Rkind),    allocatable :: R1Vec(:),R2Vec(:)
    real (kind=Rkind),   parameter   :: ZeroTresh    = ONETENTH**10

    !====================================================================
    ! Tests the sorting
    !
    ! define the matrices
    R1Vec = [ZERO,ZERO,ONE,TWO,THREE,FIVE,TEN] ! sorted vector
    R2Vec = [TWO,ZERO,FIVE,ZERO,TEN,ONE,THREE] ! unsorted vector

    ! tests
    CALL Initialize_Test(test_var,test_name='Vector')

    CALL Sort_Vec(R2Vec)
    res_test = all(abs(R1Vec-R2Vec) < ZeroTresh)
    CALL Logical_Test(test_var,test1=res_test,info='ascending sort')
    IF (.NOT. res_test) THEN
      CALL Write_Vec(R1Vec,out_unit, nbcol=7, info='R1Vec')
      CALL Write_Vec(R2Vec,out_unit, nbcol=7, info='R2Vec')
    END IF

    R1Vec = [TEN,FIVE,THREE,TWO,ONE,ZERO,ZERO] ! sorted vector
    R2Vec = [TWO,ZERO,FIVE,ZERO,TEN,ONE,THREE] ! unsorted vector

    CALL Sort_Vec(R2Vec,sort_type=-1)

    res_test = all(abs(R1Vec-R2Vec) < ZeroTresh)
    CALL Logical_Test(test_var,test1=res_test,info='descending sort')
    IF (.NOT. res_test) THEN
      CALL Write_Vec(R1Vec,out_unit, nbcol=7, info='R1Vec')
      CALL Write_Vec(R2Vec,out_unit, nbcol=7, info='R2Vec')
    END IF

    CALL Flush_Test(test_var)
    !====================================================================

    ! finalize the tests
    CALL Finalize_Test(test_var)
  END SUBROUTINE Test_QDUtil_Vector
END MODULE QDUtil_Vector_m
