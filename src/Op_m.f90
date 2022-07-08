module Op_m
  USE NumParameters_m
  USE Basis_m, only : Basis_t
  implicit none
  private

  TYPE :: Op_t

    TYPE (Basis_t),    pointer     :: Basis

    real (kind=Rk),    allocatable :: RMat(:,:)
  END TYPE Op_t

  public :: Op_t,write_Op,Set_Op,dealloc_Op,calc_OpPsi,Calc_Hpsi  ,Kpsi_nD

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



  SUBROUTINE Set_Op(Op,Basis)
   USE Basis_m
   USE Molec_m

   TYPE(Op_t),     intent(inout)       :: Op
   TYPE (Basis_t), intent(in),  target :: Basis
   REAL(kind=Rk)  ,ALLOCATABLE      :: d0bgw(:,:)


   integer :: ib,iq,i2,j2

   real (kind=Rk), allocatable :: mat_pot_grid(:,:,:),OpPsi_g(:),OpPsi_ge(:,:),Q(:),V(:,:)
    REAL(kind=Rk)              :: Q1

   IF (.NOT. Basis_IS_allocated(Basis)) THEN
     STOP 'ERROR in Set_Op: the Basis is not initialized'
   END IF

   CALL alloc_Op(Op,Basis%nb)

   Op%Basis => Basis

   ! calculation of a potential on the grid
       allocate(d0bgw(Basis%nq,Basis%nb))
       allocate(Q(1))
       allocate(mat_pot_grid(Basis%tab_basis(1)%nq , Basis%tab_basis(2)%nb , Basis%tab_basis(2)%nb))
       allocate( V( Basis%tab_basis(2)%nb , Basis%tab_basis(2)%nb)  )

         Do iq=1,Basis%tab_basis(1)%nq
            Q=Basis%tab_basis(1)%x(iq)
            CALL sub_pot(mat_pot_grid(iq,:,:),Q)
         END DO

   ! ==============================END PRINT V(:,1,1)
   do iq = 1,Basis%tab_basis(1)%nq
           write(25,*) Basis%tab_basis(1)%x(iq),mat_pot_grid(iq,1,1)
   end do
      allocate(  OpPsi_g(Basis%nb) ,  OpPsi_ge( Basis%tab_basis(1)%nq,Basis%tab_basis(2)%nb ) )
      OpPsi_g(:)   = 0
      OpPsi_ge( :,: ) = 0
      ! calculation of Op|b_i>
    DO ib=1,Basis%nb
        DO i2=1,Basis%tab_basis(2)%nb
        DO j2=1,Basis%tab_basis(2)%nb
            DO iq=1,Basis%nq
             OpPsi_ge(iq,i2) = OpPsi_ge(iq,i2)+mat_pot_grid(iq,i2,j2)*Basis%tab_basis(1)%d0gb(iq,ib)
             !OpPsi_ge(iq,i2) = OpPsi_ge(iq,i2) + V(i2,j2) * Basis%tab_basis(1)%d0gb(iq,ib)
            END DO

        END DO

        END DO

      OpPsi_g(:)= reshape( OpPsi_ge,[Basis%tab_basis(1)%nq*Basis%tab_basis(2)%nb])
          DO iq=1,Basis%nq
             OpPsi_g(iq) = OpPsi_g(iq) - HALF/mass * Basis%tab_basis(1)%d2gb(iq,ib,1,1) ! -1/2mass d2./dx2 part
          END DO
      !  CALL TRANSPOS(d0bgw,Basis%tab_basis(1))
      Op%RMat(:,ib) = matmul(Basis%d0bgw,OpPsi_g)
    END DO

   CALL write_Op(Op)


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
  SUBROUTINE Calc_Hpsi(psi_in,psi_out,Basis)
        USE Basis_m
        USE Molec_m

        COMPLEX(kind=Rk),intent(in)                   :: psi_in(:)
        COMPLEX(kind=Rk),intent(inout)                :: psi_out(:)
        TYPE (Basis_t), intent(in),  target           :: Basis
        COMPLEX(kind=Rk),ALLOCATABLE                  ::VPSI(:),KPSI(:), VPSIGB_E(:,:), KPSIGB_E(:,:),PSIGB_E(:,:)
        REAL(kind=Rk), allocatable                    :: V(:,:,:), Q(:,:)
        INTEGER                                       ::iq,i2,j2,i1 ,Ndim
        IF(.not. allocated(Basis%tab_basis)) THEN
            STOP 'ERROR in Set_Op: the Basis%tab_bais is not initialized'
        END IF
        ! calculation of action a potential VPSI_E.PSI
         Ndim = size(Basis%tab_basis)-1
          ALLOCATE(VPSI(Basis%nq*Basis%tab_basis(Ndim+1)%nb))
          ALLOCATE(VPSIGB_E(Basis%nq,Basis%tab_basis(Ndim+1)%nb))
          ALLOCATE(KPSI(  Basis%nq*Basis%tab_basis(Ndim+1)%nb  ))
          ALLOCATE(KPSIGB_E(Basis%nq,Basis%tab_basis(Ndim+1)%nb))
          ALLOCATE(PSIGB_E(Basis%nq,Basis%tab_basis(Ndim+1)%nb))
          ALLOCATE(Q(Basis%nq,Ndim))
          ALLOCATE(V(Basis%nq,Basis%tab_basis(Ndim+1)%nb,Basis%tab_basis(Ndim+1)%nb))

          call Calc_Q_grid(Q,Basis)
          Do iq=1,Basis%nq
              CALL sub_pot(V(iq,:,:),Q(iq,:))
          END DO
          PSIGB_E(:,:)= reshape( psi_in,[Basis%nq,Basis%tab_basis(Ndim+1)%nb])

          VPSIGB_E(:,:) = 0
          KPSIGB_E(:,:) =0

          DO i2=1,Basis%tab_basis(Ndim+1)%nb
             DO j2=1,Basis%tab_basis(Ndim+1)%nb
                 VPSIGB_E(:,i2) = VPSIGB_E(:,i2) + V(:,i2,j2)*PSIGB_E(:,j2)
             END DO

          END DO
          VPSI(:)= reshape( VPSIGB_E,[Basis%nq*Basis%tab_basis(Ndim+1)%nb])
          call Kpsi_nD(KPsi,psi_in,Basis)
          psi_out(:) = VPSI(:)+KPSI(:)

          DEALLOCATE(VPSI)
          DEALLOCATE(VPSIGB_E)
          DEALLOCATE(KPSI)
          DEALLOCATE(KPSIGB_E)
          DEALLOCATE(PSIGB_E)
          DEALLOCATE(Q)
  END SUBROUTINE Calc_Hpsi

  SUBROUTINE Kpsi_nD(KPsi,psi,Basis)
      USE Basis_m
      USE UtilLib_m
      USE Molec_m
      TYPE(Basis_t),        intent(in),target      :: Basis
      complex (kind=Rk), intent(in) ,target        :: Psi(:)
      complex (kind=Rk), intent(inout),target      :: KPsi(:)
      complex (kind=Rk), pointer                   :: GGB(:,:,:)
      real (kind=Rk),    pointer                   :: d2gg(:,:)
      complex (kind=Rk), pointer                   :: KGGB(:,:,:)
      logical,           parameter                 :: debug = .false.
      integer                                      :: iq,i1,i3,inb
      integer , allocatable                        :: Iq1(:),Iq2(:),Iq3(:),ndim


      IF (debug) THEN
        write(out_unitp,*) 'BEGINNING Kpsi'
        flush(out_unitp)
      END IF
      ndim = size(Basis%tab_basis)-1

      allocate(Iq3(ndim))
      allocate(Iq2(ndim))
      allocate(Iq1(ndim))

      Kpsi(:) = CZERO
      DO inb = 1,ndim

        Iq1(inb) = Product(Basis%tab_basis(1:inb-1)%nq)
        Iq2(inb) = Basis%tab_basis(inb)%nq
        Iq3(inb) = Product(Basis%tab_basis(ndim:inb+1:-1)%nq)
        KGGB(1:Iq1(inb),1:Iq2(inb),1:iq3(inb))=> KPsi
        GGB(1:Iq1(inb),1:Iq2(inb),1:iq3(inb))  => Psi
        d2gg(1:Basis%tab_basis(inb)%nq,1:Basis%tab_basis(inb)%nq)=>Basis%tab_basis(inb)%d2gg

        DO i3=1,ubound(GGB,dim=3)
           DO i1=1,ubound(GGB,dim=1)

              KGGB(i1,:,i3) =  KGGB(i1,:,i3) -(HALF/mass)*matmul(d2gg,GGB(i1,:,i3))

           END DO
        END DO
      END DO

       IF (debug) THEN
       	write(out_unitp,*) 'END KPsi_nD'
       	flush(out_unitp)
       END IF

  END SUBROUTINE  Kpsi_nD

end module Op_m
