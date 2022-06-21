module psi_m
  USE NumParameters_m
  USE Basis_m, only : Basis_t
  implicit none

  TYPE :: psi_t
    TYPE (Basis_t),    pointer     :: Basis
    real (kind=Rk),    allocatable :: RVec(:)
    complex (kind=Rk), allocatable :: CVec(:)
  END TYPE psi_t
   TYPE  :: psi0_t
     real(KIND=Rk)      ::Q0
     real(KIND=Rk)      ::K
     real(KIND=Rk)      ::phase
     real(KIND=Rk)      ::dQ
     integer            ::nb_GWP
     integer            :: ndim
     real(KIND=Rk)      :: Coef
     integer            ::  I_ElecChannel
     TYPE(psi0_t),allocatable   :: nd_psi0(:)

   END TYPE  psi0_t




   public :: psi_t,write_psi,init_psi,dealloc_psi,Calc_Norm, Calc_Norm_Grid,Norm_psi,Calc_dot_product
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
  If(grid)THEN
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





    SUBROUTINE Norm_psi(psi_in,Basis,grid,Norm)
        USE UtilLib_m
        USE Basis_m
        COMPLEX(KIND=Rk),  intent(in)                :: psi_in(:)
        COMPLEX(KIND = Rk) , ALLOCATABLE             :: psi_e(:,:)
        LOGICAL,intent(in)                           :: grid
        TYPE (Basis_t), intent(in)  ,target          ::  Basis
        REAL(KIND = Rk),intent(inout)                :: Norm(:)
        INTEGER                                      :: ib2



            Norm = 0
         If(grid)THEN
             ALLOCATE(psi_e(Basis%tab_basis(1)%nq,Basis%tab_basis(2)%nb))
             psi_e(:,:) = RESHAPE(psi_e,SHAPE= [Basis%tab_basis(1)%nq,Basis%tab_basis(2)%nb])
             Do ib2 = 1,Basis%tab_basis(2)%nb
                Norm(ib2) = sqrt(real(dot_product(psi_e(:,ib2)*Basis%tab_basis(1)%W,psi_e(:,ib2)), kind=Rk))


             END DO
         Else
           ALLOCATE(psi_e(Basis%tab_basis(1)%nb,Basis%tab_basis(2)%nb))
           psi_e(:,:) = RESHAPE(psi_e,SHAPE= [Basis%tab_basis(1)%nb,Basis%tab_basis(2)%nb])
             Do ib2 = 1,Basis%tab_basis(2)%nb
                Norm(ib2) =  sqrt(real(dot_product(psi_e(:,ib2),psi_e(:,ib2)), kind=Rk))
             End Do
           print*, Norm
          DEALLOCATE(psi_e)
        End If

    END SUBROUTINE Norm_psi


  SUBROUTINE Calc_Norm_Grid(G, Norm,Basis)
    USE UtilLib_m
    USE Basis_m
    COMPLEX(KIND=Rk),  intent(in)                :: G(:)
    COMPLEX(KIND = Rk) , ALLOCATABLE             :: G1(:,:)
    REAL(kind = Rk),intent(inout)                :: Norm
    TYPE (Basis_t), INTENT(IN)                   ::  Basis
    REAL(KIND = Rk) , ALLOCATABLE                :: Norme(:)
    INTEGER                                      :: ib2

       Norm = 0

    ALLOCATE(Norme(Basis%tab_basis(2)%nb))
    ALLOCATE(G1(Basis%tab_basis(1)%nq,Basis%tab_basis(2)%nb))
      G1(:,:) = RESHAPE(G,SHAPE= [Basis%tab_basis(1)%nq,Basis%tab_basis(2)%nb])

     Do ib2 = 1,Basis%tab_basis(2)%nb
            Norme(ib2) = sqrt(real(dot_product(G1(:,ib2)*Basis%tab_basis(1)%W,G1(:,ib2)), kind=Rk))
            !sqrt(real(dot_product(G1(:,ib2)*Basis%tab_basis(1)%W,G1(:,ib2)), kind=Rk))
          Norm = Norm+Norme(ib2)

     END DO

     print*, Norm
    print*, Norme
     DEALLOCATE(G1)
      DEALLOCATE(Norme)
END SUBROUTINE Calc_Norm_Grid




SUBROUTINE Calc_dot_product(G, dot_prdct,Basis,grid,yes)
  USE UtilLib_m
  USE Basis_m
  COMPLEX(KIND=Rk),  intent(in)                :: G(:)
  COMPLEX(KIND = Rk) , ALLOCATABLE             :: G1(:,:)
  REAL(kind = Rk),intent(inout)                :: dot_prdct
  LOGICAL                                      :: grid,yes
  TYPE (Basis_t), INTENT(IN)                   ::  Basis
  INTEGER                                      :: i_state


  !=====================================<psi|psi>======================================

  !===================================allocation========================================

    if(grid)THEN
        ALLOCATE(G1(Basis%tab_basis(1)%nq,Basis%tab_basis(2)%nb))
        G1(:,:) = RESHAPE(G,SHAPE= [Basis%tab_basis(1)%nq,Basis%tab_basis(2)%nb])


  !================================Calculation dot_product for each state=================
        dot_prdct= ZERO

            Do i_state = 1,Basis%tab_basis(2)%nb
              dot_prdct = dot_prdct+real(dot_product(G1(:,i_state)*Basis%tab_basis(1)%W,G1(:,i_state)), kind=Rk)

            END DO

          !  print*, '<psi|psi> =',dot_prdct
       DEALLOCATE(G1)
    else
      dot_prdct = real( dot_product(G,G),kind=Rk)

    end if
    if(yes)  write(*,*)'<psi|psi> =',dot_prdct



END SUBROUTINE Calc_dot_product


end module psi_m
