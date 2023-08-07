PROGRAM App_QDLib
  USE QDUtil_m
  IMPLICIT NONE

  integer                          :: i,n
  real(kind=Rkind),    allocatable :: RMat(:,:),REigVal(:),REigVec(:,:),RVec(:)
  real(kind=Rkind),    allocatable :: RVec2(:)

  TYPE(Frac_t)                     :: Frac1, Frac2
  TYPE(Frac_t),        allocatable :: tab_Frac(:)

  !  #if __LAPACK != 1
  !    write(out_unit,*) '  Lapack library is linked'
  !  #else
  !    write(out_unit,*) '  Lapack library is not linked'
  !  #endif
 !====================================================================
  ! Tests on fractions
  Frac1 = '1/-2' ! use the conversion from string to Frac_t
  write(*,*) 'Frac1: ',TO_String(Frac1) ! it give "Frac1: -1/2"
  Frac2 = -2*Frac1 ! here the result is one and it is simplified
  write(*,*) 'Frac2: ',TO_String(Frac2) ! it give "Frac2: 1"
  Frac2 = Frac1**3
  write(*,*) 'Frac2: ',TO_String(Frac2) ! it give "Frac2: -1/8"
  tab_Frac = Frac_t(1,[2,3,4])
  write(*,*) 'tab_Frac: ',(TO_String(tab_Frac(i)) // ' ',i=1,size(tab_Frac)) ! it give "tab_Frac: 1/2 1/3 1/4 "

  !====================================================================
  ! Tests for the matrix digonalization
  !
  ! define the matrices
  n = 3
  RMat =  reshape([ONE,HALF,ZERO,                             &
                   HALF,ONE,HALF,                             &
                   ZERO,HALF,ONE],shape=[n,n])

  allocate(REigVal(n))
  allocate(REigVec(n,n))


  CALL diagonalization(RMat,REigVal,REigVec)


  CALL Write_Mat(RMat,out_unit,5,info='RMat')
  write(out_unit,*)
  CALL Write_Mat(REigVec,out_unit,5,info='REigVec (in column)')
  write(out_unit,*)
  CALL Write_Vec(REigVal,out_unit,5,info='REigVal')
  write(out_unit,*)

  DO i=1,n
    write(out_unit,*) i,matmul(Rmat,REigVec(:,i))-REigVal(i)*REigVec(:,i)
  END DO


 !====================================================================
  ! Tests for identity matrix
  !
  ! define the matrices
  write(out_unit,*) 'Test identity matrix'
  n = 100
  RMat =  Identity_Mat(n)
  RVec  = [(real(i,kind=Rkind)/Pi,i=1,n)]
  RVec2 = matmul(Rmat,RVec)


  write(out_unit,*) 'RVec(1) and RVec(n)',RVec(1),RVec(n)
  write(out_unit,*) 'diff',maxval(abs((RVec-RVec2)))

END PROGRAM App_QDLib