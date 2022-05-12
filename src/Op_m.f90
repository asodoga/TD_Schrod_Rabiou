module Op_m
  USE NumParameters_m
  USE Basis_m, only : Basis_t
  implicit none
  private

  TYPE :: Op_t

    TYPE (Basis_t),    pointer     :: Basis

    real (kind=Rk),    allocatable :: RMat(:,:)
  END TYPE Op_t

  public :: Op_t,write_Op,Set_Op,dealloc_Op,calc_OpPsi

contains
  SUBROUTINE alloc_Op(Op,nb)
    TYPE(Op_t),  intent(inout) :: Op
    integer,     intent(in)    :: nb

    CALL dealloc_Op(Op)

    IF (nb < 1) STOP 'ERROR in init_Op: nb < 1!'

    allocate(Op%RMat(nb,nb))

  END SUBROUTINE alloc_Op

  SUBROUTINE dealloc_Op(Op)
    TYPE(Op_t), intent(inout) :: Op


    nullify(Op%Basis)

    IF (allocated(Op%RMat)) THEN
      deallocate(Op%RMat)
    END IF

  END SUBROUTINE dealloc_Op

  SUBROUTINE write_Op(Op)
    USE UtilLib_m, only : Write_RMat

    TYPE(Op_t), intent(in) :: Op

   ! integer :: i

    IF (associated(Op%Basis)) THEN
      write(out_unitp,*) ' The basis is linked to Op.'
    END IF

    IF (allocated(Op%RMat)) THEN
      write(out_unitp,*) 'Writing Op (real):'
      write(out_unitp,*)
      CALL Write_RMat(Op%RMat,out_unitp,5,name_info='Op%Rmat')
      write(out_unitp,*) 'END Writing Op'
    END IF

  END SUBROUTINE write_Op

  SUBROUTINE Set_Op(Op,psi_in,Basis)
  USE Basis_m
  USE Molec_m

    TYPE(Op_t),     intent(inout)       :: Op
    COMPLEX (kind=Rk), INTENT(IN)       :: psi_in(:)
    TYPE (Basis_t), INTENT(IN)  ,target ::  Basis
    logical,         parameter          :: debug = .false.
    COMPLEX (kind=Rk), allocatable      :: OpPsi_gb(:),Psi_bb(:),Psi_gb(:),Psi_bb1(:,:)
    COMPLEX (kind=Rk), allocatable      ::Psi_gb1(:,:),Vpsigb(:,:),Kpsigb(:,:)
    real (kind=Rk), allocatable         :: Q(:),mat_pot_grid(:,:,:),  Vpsi_out(:),Kpsi_out(:)
    integer                             :: ib,iq,iq1,iq2,jb,ib1,j2,i2,i1
    nullify(Op%Basis)

    IF (debug) THEN
      write(out_unitp,*) 'BEGINNING set_Op'
      call write_Op(Op)
      call write_basis(Op%Basis)
      flush(out_unitp)
    END IF
    allocate(Op%RMat(Op%Basis%nb,Op%Basis%nb))


     allocate(Q(1))
     allocate(Psi_gb(Op%Basis%nq))
      allocate(OpPsi_gb(Op%Basis%nq))
     allocate(Psi_gb1(Basis%tab_basis(1)%nq,Basis%tab_basis(2)%nb))
     allocate(mat_pot_grid(Op%Basis%nq,2,2))
     Do iq=1,Basis%tab_basis(1)%nq
       Q=Basis%tab_basis(1)%x(iq)
      CALL sub_pot(mat_pot_grid,Q)
     END DO
     Psi_gb1 = reshape(psi_in,[ Basis%tab_basis(1)%nq, Basis%tab_basis(2)%nb])
    Vpsigb(:,:) = 0
    DO i2=1,Basis%tab_basis(2)%nb
    DO j2=1,Basis%tab_basis(2)%nb
       Vpsigb(:,j2) = Vpsigb(:,j2) + DOT_PRODUCT(  mat_pot_grid(:,j2,i2), Psi_gb1(:,i2) )
         DO i1=1,Basis%tab_basis(1)%nb
            Kpsigb(:,i2) = -HALF/mass *DOT_PRODUCT( Basis%tab_basis(1)%d2gb(:,i1,1,1),Psi_gb1(:,i2))
         END DO

    END DO
    END DO

    Kpsi_out(:) = reshape(Kpsigb,[Basis%tab_basis(1)%nq*Basis%tab_basis(2)%nb])
    Vpsi_out(:) = reshape(Vpsigb,[Basis%tab_basis(1)%nq*Basis%tab_basis(2)%nb])
    OpPsi_gb(:) =   Vpsi_out(:)+  Kpsi_out(:)
    OpPsi_gb(:) = OpPsi_gb(:) * Basis%tab_basis(1)%w(:)
         DO ib=1,Basis%nb
            Op%RMat(:,ib) = matmul(transpose(Basis%d0gb),OpPsi_gb)
         END DO

    CALL write_Op(Op)

    deallocate(Op%RMat)
    deallocate(Q)
    deallocate(Psi_gb)
    deallocate(OpPsi_gb)
    deallocate(Psi_gb1)
    deallocate(mat_pot_grid)

    IF (debug) THEN
      write(out_unitp,*) 'END set_Op'
      flush(out_unitp)
    END IF
  END SUBROUTINE Set_Op
  SUBROUTINE calc_OpPsi(Op,Psi,OpPsi)
    USE psi_m, ONLY : psi_t

    TYPE(Op_t),  intent(in)     :: Op
    TYPE(psi_t), intent(in)     :: Psi
    TYPE(psi_t), intent(inout)  :: OpPsi


    IF (.NOT. allocated(Op%RMat)) THEN
      STOP 'ERROR in calc_OpPsi: Op is not initialized!'
    END IF

    IF (allocated(Psi%RVec)) THEN
      OpPsi%RVec = matmul(Op%RMat,Psi%RVec)
    ELSE IF (allocated(Psi%CVec)) THEN
      OpPsi%CVec = matmul(Op%RMat,Psi%CVec)
    ELSE
      STOP 'ERROR in calc_OpPsi: Psi is not initialized!'
    END IF


  END SUBROUTINE calc_OpPsi








end module Op_m
