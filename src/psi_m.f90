module psi_m
  USE NumParameters_m
  USE Basis_m, only : Basis_t
  implicit none

  private

  TYPE :: psi_t
    TYPE (Basis_t),    pointer     :: Basis
    real (kind=Rk),    allocatable :: RVec(:)
    complex (kind=Rk), allocatable :: CVec(:)
  END TYPE psi_t

   public :: psi_t,write_psi,init_psi,dealloc_psi,Calc_Norm, Calc_Norm_Grid
  ! operation on psi has to be defined: psi=psi1, psi1+psi2, psi=psi1*cte ...
contains
  SUBROUTINE init_psi(psi,Basis,cplx,grid)
  USE Basis_m

    TYPE(psi_t),    intent(inout)      :: psi
    TYPE (Basis_t), intent(in), target :: Basis
    logical,        intent(in)         :: cplx,grid

    CALL dealloc_psi(psi)

    IF (.NOT. Basis_IS_allocated(Basis)) STOP 'ERROR in init_psi: the Basis is not initialized'


    IF (Basis%nb < 1) STOP 'ERROR in init_psi: Basis%nb < 1!'

    psi%Basis => Basis
  If(grid)THEN   !grid
    IF (cplx) THEN
     IF(allocated(Basis%tab_basis))THEN
      allocate(psi%CVec(Basis%tab_basis(1)%nq*Basis%tab_basis(2)%nb))
     else
      allocate(psi%CVec(Basis%nq))
     END IF
    ELSE
      IF(allocated(Basis%tab_basis))THEN
        allocate(psi%RVec(Basis%tab_basis(1)%nq*Basis%tab_basis(2)%nb))
      ELSE
        allocate(psi%RVec(Basis%nq))
      END IF
    END IF
  ELSE ! grid
    IF (cplx) THEN
     IF(allocated(Basis%tab_basis))THEN
      allocate(psi%CVec(Basis%tab_basis(1)%nq*Basis%tab_basis(2)%nb))
     else
      allocate(psi%CVec(Basis%nb))
     END IF
    ELSE
      IF(allocated(Basis%tab_basis))THEN
        allocate(psi%RVec(Basis%tab_basis(1)%nb*Basis%tab_basis(2)%nb))
      ELSE
        allocate(psi%RVec(Basis%nb))
      END IF
    END IF
  END IF
  END SUBROUTINE init_psi

  SUBROUTINE dealloc_psi(psi)
    TYPE(psi_t), intent(inout) :: psi

    nullify(psi%Basis)

    IF (allocated(psi%RVec)) THEN
      deallocate(psi%RVec)
    END IF
    IF (allocated(psi%CVec)) THEN
      deallocate(psi%CVec)
    END IF

  END SUBROUTINE dealloc_psi

  SUBROUTINE write_psi(psi)
    TYPE(psi_t), intent(in) :: psi
    integer                 :: i

    IF (associated(psi%Basis)) THEN
      write(out_unitp,*) ' The basis is linked to psi.'
    END IF

    IF (allocated(psi%RVec)) THEN
      write(out_unitp,*) 'Writing psi (real):'
      write(out_unitp,*) psi%RVec
      write(out_unitp,*) 'END Writing psi'
    END IF
    IF (allocated(psi%CVec)) THEN
      write(out_unitp,*) 'Writing psi (complex):'
      do i=1, size(psi%CVec)
      write(out_unitp,*) i, psi%CVec(i)
      end do
      write(out_unitp,*) 'END Writting psi'
    END IF

  END SUBROUTINE write_psi


  SUBROUTINE Calc_Norm(psi, Norm)
  TYPE (psi_t),  intent(in)     :: psi
  real(kind = Rk),intent(inout) :: Norm

  !IF (allocated(psi%CVec)) THEN
  !Norm = sqrt(real(dot_product(psi%CVec,psi%CVec), kind=Rk))
  !END IF

  IF (allocated(psi%CVec)) THEN
   Norm = sqrt(real(dot_product(psi%CVec,psi%CVec), kind=Rk))
  END IF


  !IF (allocated(psi%RVec)) THEN
   !Norm = sqrt(dot_product(psi%RVec,psi%RVec), kind=Rk)
  !END IF
  !write(out_unitp,*) 'norm,psi',Norm

  END SUBROUTINE Calc_Norm




  SUBROUTINE Calc_Norm_Grid(G, Norm,Basis)
    USE UtilLib_m
    USE Basis_m
    COMPLEX(KIND=Rk),  intent(in)                :: G(:)
    COMPLEX(KIND = Rk) , ALLOCATABLE           :: G1(:,:)
    REAL(kind = Rk),intent(inout)              :: Norm
    TYPE (Basis_t), INTENT(IN)  ,target      ::  Basis
    REAL(KIND = Rk) , ALLOCATABLE              :: Norme(:)
    INTEGER                                    :: ib2

       Norm = 0

    ALLOCATE(Norme(2))
    ALLOCATE(G1(Basis%tab_basis(1)%nq,Basis%tab_basis(2)%nb))
      G1(:,:) = RESHAPE(G,SHAPE= [Basis%tab_basis(1)%nq,Basis%tab_basis(2)%nb])

     Do ib2 = 1,Basis%tab_basis(2)%nb
            Norme(ib2) = sqrt(real(dot_product(G1(:,ib2)*Basis%tab_basis(1)%W,G1(:,ib2)), kind=Rk))

        Norm = Norm+Norme(ib2)

     END DO
      Norm = SQRT(Norm)
     print*, Norm
     DEALLOCATE(G1)
      DEALLOCATE(Norme)
END SUBROUTINE Calc_Norm_Grid

end module psi_m
