module Op_m
   !$ USE omp_lib
   USE QDUtil_m
   USE Basis_m, only: Basis_t
   USE Molec_m
   Use Psi_m
   implicit none
   private

   TYPE :: Op_t

      TYPE(Basis_t), pointer           :: Basis
      real(kind=Rkind), allocatable    :: RMat(:, :)
      complex(kind=Rkind), allocatable :: CMat(:, :)
      real(kind=Rkind), allocatable    :: Scal_pot(:, :, :)

   END TYPE Op_t

   public :: Op_t, write_Op, Set_Op, dealloc_Op, calc_OpPsi, Calc_Hpsi, Kpsi_nD, Make_Mat_H
   public  :: test_op,Calc_Scalar_Pot, test_openmp_op,calc_tab_Iq

contains
   SUBROUTINE alloc_Op(Op, nb)
      TYPE(Op_t), intent(inout) :: Op
      integer, intent(in)    :: nb

      CALL dealloc_Op(Op)

      IF (nb < 1) STOP 'ERROR in init_Op: nb < 1!'

      allocate (Op%RMat(nb, nb))

   END SUBROUTINE alloc_Op

   SUBROUTINE dealloc_Op(Op)
      TYPE(Op_t), intent(inout) :: Op

      nullify (Op%Basis)

      IF (allocated(Op%RMat)) THEN
         deallocate (Op%RMat)
      END IF

   END SUBROUTINE dealloc_Op

   SUBROUTINE write_Op(Op)
      USE QDUtil_m

      TYPE(Op_t), intent(in) :: Op

      ! integer :: i

      IF (associated(Op%Basis)) THEN
         write (out_unit, *) ' The basis is linked to Op.'
      END IF

      IF (allocated(Op%RMat)) THEN
         write (out_unit, *) 'Writing Op (real):'
         write (out_unit, *)
         CALL Write_VecMat(Op%RMat, out_unit, 5,  info='Op%Rmat')
         write (out_unit, *) 'END Writing Op'
      END IF

   END SUBROUTINE write_Op

   SUBROUTINE Make_Mat_H(Basis,H)
      Use QDUtil_m
      USE Basis_m
      USE Psi_m

      TYPE(Basis_t),intent(in)            :: Basis
      TYPE(Op_t), intent(inout)           :: H

      !logical, parameter                  :: debug = .true.
      logical,         parameter          :: debug = .false.
      integer, allocatable                :: Tab_iq(:, :)
      TYPE(psi_t)                         :: psi,Hpsi
      integer                             :: nb,nsurf,ib,jb,ndim

      IF (debug) THEN
         write (out_unit, *) 'BEGINNING Make_Mat_H'
         call write_basis(Basis)
         flush (out_unit)
      END IF

      call  init_psi(psi, Basis,  cplx=.true., grid=.false.)
      call  init_psi(Hpsi,Basis,  cplx=.true., grid=.false.)

      ndim = size(Basis%tab_basis) - 1
      nsurf=Basis%tab_basis(ndim+1)%nb
      nb =(Basis%nb)*nsurf
      allocate(H%CMat(nb,nb))
      call Calc_tab_Iq0(Tab_Iq,Basis)
      call Set_Op(H, Basis,Tab_Iq)

      DO jb = 1, nb
         psi%CVec(:)  = CZERO
         psi%CVec(jb) = CONE
         call  calc_Oppsi(H, psi, Hpsi)

         DO ib = 1, nb
            psi%CVec(:) = CZERO
            psi%CVec(ib) = CONE
            H%CMat(ib, jb) = dot_product(psi%CVec,Hpsi%CVec)
         END DO

      END DO
      call dealloc_psi(psi)
      call dealloc_psi(Hpsi)
      deallocate(Tab_Iq)

      IF (debug) THEN
         write (out_unit, *) 'END Make_Mat_H'
         flush (out_unit)
      END IF

   END SUBROUTINE 

   SUBROUTINE Set_Op(Op, Basis,Tab_Iq)
      USE Basis_m

      TYPE(Op_t), intent(inout)         :: Op
      TYPE(Basis_t), intent(in), target :: Basis
      integer, intent(in)               :: Tab_iq(:, :)

      Integer                            :: ndim,nsurf,nq

      IF (.NOT. Basis_IS_allocated(Basis)) THEN
         STOP 'ERROR in Set_Op: the Basis is not initialized'
      END IF

      ndim = size(Basis%tab_basis) - 1
      nq =Basis%nq
      nsurf=Basis%tab_basis(ndim+1)%nb

      !CALL alloc_Op(Op, Basis%nb*Basis%tab_basis(size(Basis%tab_basis))%nb)
      !allocate(Op%Scal_pot(nq, nsurf, nsurf))
      !Op%Basis => Basis

      !call Calc_Scalar_Pot(Op%Scal_pot, Basis)
      !call Calc_tab_Iq0(Tab_Iq,Basis)
      call Calc_Scalar_Pot_openmp(Op%Scal_pot, Basis,Tab_Iq)

   END SUBROUTINE 

   SUBROUTINE Kpsi_nD(KPsi_g, Psi_g, Basis)
      USE Basis_m
      USE  QDUtil_m
      USE Molec_m
      TYPE(Basis_t), intent(in), target               :: Basis
      complex(kind=Rkind), intent(in), target         :: psi_g(:)
      complex(kind=Rkind), intent(inout), target      :: Kpsi_g(:)
      complex(kind=Rkind), pointer                    :: psi_ggb(:, :, :)
      complex(kind=Rkind), pointer                    :: d2gg(:, :)
      complex(kind=Rkind), pointer                    :: Kpsi_ggb(:, :, :)
      real(kind=Rkind), allocatable                   :: GGdef(:, :)
      logical, parameter                              :: debug = .true.
      integer                                         :: iq, i1, i3, inb, ndim
      integer, allocatable                            :: Iq1, Iq2, Iq3

      IF (debug) THEN
         !write(out_unit,*) 'BEGINNING Kpsi'
         flush (out_unit)
      END IF
      
      ndim = size(Basis%tab_basis)-1
      allocate (GGdef(ndim , ndim))
      CALL get_Qmodel_GGdef(GGdef)

      Kpsi_g(:) = CZERO
      DO inb = 1, ndim 
         
      IF (inb == 1) THEN
      
          Iq1 = 1
          Iq2 = Basis%tab_basis(1)%nq
          Iq3 = Product(Basis%tab_basis(2:ndim)%nq)*Basis%tab_basis(ndim + 1)%nb

      ELSE IF (inb == ndim) THEN

         Iq1  = Product(Basis%tab_basis(1:ndim - 1)%nq)
         Iq2  = Basis%tab_basis(ndim)%nq
         Iq3  = Basis%tab_basis(ndim + 1)%nb

      ELSE
          Iq1 = Product(Basis%tab_basis(1:inb - 1)%nq)
          Iq2 = Basis%tab_basis(inb)%nq
          Iq3 = Product(Basis%tab_basis(inb + 1:ndim)%nq)*Basis%tab_basis(ndim + 1)%nb
      END IF

         Kpsi_ggb(1:Iq1, 1:Iq2, 1:Iq3) => Kpsi_g
         psi_ggb(1:Iq1, 1:Iq2, 1:Iq3) => psi_g
         d2gg(1:Iq2, 1:Iq2) => Basis%tab_basis(inb)%d2gg

         DO i3 = 1, ubound(psi_ggb, dim=3)
            DO i1 = 1, ubound(psi_ggb, dim=1)
               !KPsi_ggb(i1, :, i3) = KPsi_ggb(i1, :, i3) - HALF*GGdef(inb, inb)*matmul(d2gg, psi_ggb(i1, :, i3))
               Kpsi_ggb(i1, :, i3) = Kpsi_ggb(i1, :, i3) - HALF*matmul(d2gg, psi_ggb(i1, :, i3))
            END DO
         END DO

      END DO
      deallocate (Iq1, Iq2, Iq3)
      IF (debug) THEN
         !        write(out_unit,*) 'END KPsi_nD'
         flush (out_unit)
      END IF
   END SUBROUTINE Kpsi_nD

SUBROUTINE Calc_Scalar_Pot(V, Basis)
   USE  QDUtil_m
   USE Basis_m
   TYPE(Basis_t), intent(in), target                :: Basis
   real(kind=Rkind), allocatable ,intent(inout)     :: V(:, :, :)

   integer, allocatable                             :: Tab_iq(:)
   integer                                          :: inb, ndim, iq,nq,nsurf,i
   real(Kind=Rkind)               , allocatable     :: Q(:)
   logical                                          :: Endloop

   ndim = size(Basis%tab_basis) - 1
   nq =Basis%nq
   nsurf=Basis%tab_basis(ndim+1)%nb

    allocate (Q(ndim),Tab_iq(ndim))
    allocate(V(nq, nsurf, nsurf))

   Call Init_tab_ind(Tab_iq, Basis%NDindexq)
   Iq = 0
   DO
      Iq = Iq + 1
     ! i =Tab_iq(2)
      CALL increase_NDindex(Tab_iq, Basis%NDindexq, Endloop)
      IF (Endloop) exit
      do inb = 1, ndim
          Q(inb) = Basis%tab_basis(inb)%x(Tab_iq(inb))
      end do
      CALL sub_Qmodel_V(V(iq, :, :), Q(:))
      !if ( (Tab_iq(2)/= i) )  write(270,*)
      !write(270,*) Q(:),V(iq, 1, 1),V(iq, 2, 2)

   END DO

   deallocate(Tab_iq,Q)
END SUBROUTINE 


SUBROUTINE test_op(Basis,psi0)
   USE  QDUtil_m
   USE Basis_m
   TYPE(Basis_t), intent(in), target             :: Basis
   TYPE(psi_t),intent(in)                        :: psi0
    TYPE(Op_t)                                   :: H
    TYPE(psi_t)                                  :: psi
    real(kind=Rkind), allocatable                :: V(:, :, :)
    complex(kind=Rkind) ,allocatable             :: CEigVal(:),CEigVec(:,:)
    real(kind=Rkind) ,allocatable                :: prob(:),vec(:)
    real(Kind=Rkind)            , allocatable    ::Pop(:)
    integer                                      :: nb,ndim,ib,nsurf

    open(unit=50, file = 'proj')
    open(unit=51, file = 'EingVec')
    open(unit=52, file = 'population')

     call  Make_Mat_H(Basis,H)
     ndim = size(Basis%tab_basis) - 1
     nb =Basis%nb*Basis%tab_basis(ndim+1)%nb
     nsurf =Basis%tab_basis(ndim+1)%nb
     allocate(CEigVal(nb),CEigVec(nb,nb),prob(nb),vec(nb),Pop(nsurf))
     call init_psi(psi, psi0%Basis, cplx=.true., grid=.false.)
     call   Write_VecMat(H%CMat,out_unit,5, info='<psi|H|psi>',Rformat='e13.4')
     call diagonalization(H%CMat,CEigVal,CEigVec)
     prob = matmul(CEigVec,psi0%CVec)
     Do ib = 1,nb
      psi%CVec(:)=CEigVec(:,ib)
      call Popu(psi,Pop)
      write(51,*) ib, CEigVal(ib)%re
      write(50,*) ib, abs(prob(ib))
      write(52,*) ib, Pop(:)
     End Do

     call Calc_Scalar_Pot(V, Basis)
     

     call dealloc_psi(psi)
END SUBROUTINE 


subroutine Popu(Psi, Pop)
   implicit none
   type(Psi_t), intent(in), target                 :: Psi
   complex(kind=Rkind), pointer                    :: Psi_bb(:, :)
   real(Kind=Rkind), intent(inout), allocatable    ::Pop(:)
   integer                                         :: inb,nsurf,ndim,nb
   real(Kind=Rkind)                                :: Norm
   ndim = size(Psi%Basis%tab_basis)-1
   nb = Psi%Basis%nb
   nsurf = Psi%Basis%tab_basis(ndim+1)%nb
   Psi_bb(1:nb, 1:nsurf) => Psi%CVec
   call Calc_Norm_OF_Psi(Psi, Norm)
   do inb = 1, nsurf
      Pop(inb) = dot_product(Psi_bb(:, inb), Psi_bb(:, inb))/Norm
   end do
end subroutine 


 SUBROUTINE Calc_Hpsi(psi_g, HPsi_g, Basis,V)
    USE Basis_m
    USE Psi_m
    USE Molec_m
    complex(kind=Rkind), intent(in), target            :: psi_g(:)
    complex(kind=Rkind), intent(inout)                 :: HPsi_g(:)
    type(Basis_t), intent(in), target                  :: Basis
    real(kind=Rkind) , intent(in)                      :: V(:, :, :)

    complex(kind=Rkind), allocatable, target           :: VPsi_g(:), KPsi_g(:)
    complex(kind=Rkind), pointer                       :: VPsi_gb(:, :), Psi_gb(:, :)
    Integer                                            :: i, j, ndim,nsurf,nq


    IF (.not. allocated(Basis%tab_basis)) THEN
       STOP 'ERROR in Set_Op: the Basis%tab_bais is not initialized'
    END IF

    ndim = SIZE(Basis%tab_basis) - 1
    nq =Basis%nq
    nsurf=Basis%tab_basis(ndim+1)%nb

    ! action potential V|Psi_g> 
    allocate(VPsi_g(nq*nsurf))
    allocate(KPsi_g(nq*nsurf))

    VPsi_g(:) = CZERO
    KPsi_g(:) = CZERO
    HPsi_g(:) = CZERO

   VPsi_gb(1:nq,1:nsurf) => VPsi_g
   Psi_gb(1:nq, 1:nsurf) => Psi_g

   DO i = 1, nsurf

      DO j = 1, nsurf

         VPsi_gb(:, i) = VPsi_gb(:, i) + V(:, i, j)*Psi_gb(:, j)

      END DO

   END DO

    ! action potential K|Psi_g>

     call Kpsi_nD(KPsi_g, Psi_g, Basis)
     HPsi_g(:) = VPsi_g(:) + KPsi_g(:)

    DEALLOCATE (VPsi_g)
    DEALLOCATE (KPsi_g)

 END SUBROUTINE 



 SUBROUTINE calc_OpPsi(Op, psi, Oppsi)

   USE psi_m
   TYPE(Op_t), intent(in)           :: Op
   TYPE(psi_t), intent(in)          :: psi
   TYPE(psi_t), intent(inout)       :: Oppsi
   TYPE(psi_t)                      :: psi_g, Oppsi_g

   IF (allocated(Op%RMat)) THEN

      Oppsi%CVec = matmul(Op%RMat, psi%CVec)

   ELSE

      IF (psi%Grid) THEN

         call init_psi(Oppsi_g, Psi%Basis, .true., .true.)
         call Calc_Hpsi(psi%CVec, Oppsi_g%CVec, psi%Basis,Op%Scal_pot)
         call GridTOBasis_nD_cplx(Oppsi%CVec, Oppsi_g%CVec, psi%Basis)
         call dealloc_psi(Oppsi_g)

      ELSE

         call init_psi(psi_g, Psi%Basis, .true., .true.)
         call init_psi(Oppsi_g, psi%Basis, .true., .true.)
         call BasisTOGrid_nD_cplx(psi_g%CVec, psi%CVec, psi%Basis)
         call Calc_Hpsi(psi_g%CVec, Oppsi_g%CVec, psi%Basis,Op%Scal_pot)
         call GridTOBasis_nD_cplx(Oppsi%CVec, Oppsi_g%CVec, psi%Basis)
         call dealloc_psi(psi_g)
         call dealloc_psi(Oppsi_g)

      END IF

   END IF

END SUBROUTINE 



SUBROUTINE Calc_tab_Iq(Tab_Iq,Basis)
   USE  QDUtil_m
   USE Basis_m
   TYPE(Basis_t), intent(in), target             :: Basis
   integer, allocatable ,intent(inout)           :: Tab_iq(:, :)
   integer, allocatable                          :: Tab_iq0(:)
   integer                                       :: ndim, iq,nq
   logical                                       :: Endloop

   ndim = size(Basis%tab_basis) - 1
   nq =Basis%nq

    allocate (Tab_iq(ndim,nq),Tab_iq0(ndim))
   Call Init_tab_ind(Tab_iq0, Basis%NDindexq)
   Iq = 0
   DO
      Iq = Iq + 1
      CALL increase_NDindex(Tab_Iq0, Basis%NDindexq, Endloop)
      IF (Endloop) exit
      Tab_iq(:,Iq) = Tab_Iq0
       !print*,iq,Tab_Iq(:,Iq)
   END DO
   deallocate(Tab_Iq0)
END SUBROUTINE 


SUBROUTINE Calc_Scalar_Pot_openmp(V, Basis,Tab_Iq)
   !$ USE omp_lib
   USE  QDUtil_m
   USE Basis_m
   TYPE(Basis_t), intent(in), target               :: Basis
   real(kind=Rkind), allocatable ,intent(inout)    :: V(:, :, :)
   integer ,intent(in)                             :: Tab_iq(:, :)

   integer                                          :: Ib, ndim, iq,nq,nsurf,maxth
   real(kind=Rkind)               , allocatable     :: Q(:)
   logical                                          :: Endloop

   ndim = size(Basis%tab_basis) - 1
   nq   =Basis%nq
   nsurf=Basis%tab_basis(ndim+1)%nb

    allocate (Q(ndim))
    allocate(V(nq, nsurf, nsurf))
    V(:,:,:) = ZERO
    Q(:) = ZERO

    maxth = 1
    !$ maxth  = omp_get_max_threads()
    write(*,*) 'nbr de coeurs :',maxth
   !$OMP   PARALLEL DEFAULT(NONE) &
   !$OMP   SHARED(Basis,Tab_Iq,maxth,V,nq,ndim,Ib) &
   !$OMP   PRIVATE(Iq,Q) &
   !$OMP   NUM_THREADS(maxth)

   !$OMP BARRIER
   !$OMP   DO SCHEDULE(STATIC)

   DO  Iq = 1,nq

      DO Ib = 1, ndim
          Q(Ib) = Basis%tab_basis(Ib)%x(Tab_Iq(Ib,Iq))
      ENDDO

      !CALL sub_pot(V(Iq, :, :), Q(:), 0)
      CALL sub_Qmodel_V(V(iq, :, :), Q(:))

   END DO
   !$OMP END DO
   !$OMP BARRIER
   !$OMP END PARALLEL

   deallocate(Q)
END SUBROUTINE 



SUBROUTINE test_openmp_op(Basis)
   USE  QDUtil_m
   USE Basis_m
   TYPE(Basis_t), intent(in), target             :: Basis
   real(kind=Rkind), allocatable                 :: V(:, :, :)
   integer, allocatable                          :: Tab_Iq(:, :)
   real(kind=Rkind)                              :: t1,t2,tps, tpsopenmp


   call Calc_tab_Iq(Tab_Iq,Basis)
     
   !call cpu_time(time=t1)
    call Calc_Scalar_Pot_openmp(V, Basis,Tab_Iq)
   ! call cpu_time(time=t2)
    !tpsopenmp = t2-t1
    !write(*,*)"temps CPU  openmp :",tpsopenmp
    deallocate(V)
    !call cpu_time(time=t1)
    !call Calc_Scalar_Pot(V, Basis)
    !call cpu_time(time=t2)
    !tps = t2-t1
    !write(*,*)"temps CPU  non parallele :",tps

END SUBROUTINE 




 
end module Op_m