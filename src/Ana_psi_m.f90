module Ana_psi_m
   USE Basis_m
   USE NDindex_m
   USE psi_m
   implicit none
   private
   public:: Population, Qpop
   public::  Calc_Av_imp_k_nD
   public :: Calc_Avg_A_nD,test_analyse_psi,Calc_AVQ_SQ_nD

contains


   SUBROUTINE test_analyse_psi(psi)
      USE  QDUtil_m
      USE Basis_m
      type(psi_t)    ,intent(in)                    :: psi

      real(kind=Rkind), allocatable                 :: V(:, :, :)
      integer, allocatable                          :: Tab_Iq(:, :)
      real(kind=Rkind) ,allocatable                 :: SQt(:),Qt(:)
      complex(kind=Rkind) ,allocatable              :: At(:)
      integer                                       :: ndim 
     

      ndim = size(psi%Basis%tab_basis) - 1
      allocate(Qt(ndim), SQt(ndim),At(ndim))
   
      call Calc_tab_Iq0(Tab_Iq,psi%Basis)
      call Calc_Avg_A_nD(psi, At)

        write(out_unit,*) 'At',At
        write(out_unit,*) '=============================================================================='
   
   END SUBROUTINE 
   
   subroutine Population(Psi, Pop)
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

   SUBROUTINE Qpop(Psi, Qp)
      USE Basis_m
      USE QDUtil_m
      type(Psi_t), intent(in), target               :: Psi
      type(Psi_t), target                           :: Psi_g
      complex(kind=Rkind), pointer                     :: psi_gb(:, :)

      real(kind=Rkind), intent(inout)                  :: Qp(:)
      real(kind=Rkind)                                 :: Norm(2)
      logical, parameter                            :: debug = .true.
      integer                                       :: iq, inb, ndim

      IF (debug) THEN
         !write(out_unit,*) 'BEGINNING Qpop'
         flush (out_unit)
      END IF
      Ndim = size(Psi%Basis%tab_basis)
      call init_psi(psi_g, psi%Basis, cplx=.TRUE., grid=.true.)
      Psi_g%CVec(:) = CZERO
      call BasisTOGrid_nD_cplx(Psi_g%CVec, Psi%CVec, Psi%Basis)
      do inb = 1, Psi%Basis%tab_basis(2)%nb
         psi_gb(1:Psi%Basis%tab_basis(1)%nq, 1:Psi%Basis%tab_basis(2)%nb) => psi_g%CVec
         Qp(inb) = dot_product(psi_gb(:, inb), Psi%Basis%tab_basis(1)%w*Psi%Basis%tab_basis(1)%x*psi_gb(:, inb))
         Norm(inb) = dot_product(psi_gb(:, inb), Psi%Basis%tab_basis(1)%w*psi_gb(:, inb))
         if (Norm(inb) /= ZERO) then
            Qp(inb) = Qp(inb)/Norm(inb)
         end if
      end do
      !Do iq = 1,Psi%Basis%tab_basis(1)%nq
      !   write(666,*)   psi_gb(iq,1) , psi_gb(iq,2)
      !End Do

      print *, Qp, Norm
      IF (debug) THEN
         !        write(out_unit,*) 'END Qpop
         flush (out_unit)
      END IF
   END SUBROUTINE Qpop


SUBROUTINE Calc_Avg_A_nD(psi, At)
   USE QDUtil_m
   type(psi_t), intent(in)                                 :: psi
   complex(kind=Rkind), intent(inout)                      :: At(:)
   
   !locals variables---------------------------------------------------
   
   real(kind=Rkind), allocatable                           :: VQ(:),VP(:),VQQ(:),VQP(:),SQt(:)
   real(kind=Rkind), allocatable                           :: CA(:),CB(:)
   integer, allocatable                                    :: Tab_Iq(:, :)
   logical, parameter                                      :: debug = .true.
   integer                                                 :: Ndim, Ib
   
   !debuging----------------------------------------------------------------
   
   IF (debug) THEN
   
      flush (out_unit)
      
   END IF
   
   Ndim = size(psi%Basis%tab_basis) - 1 
   call Calc_tab_Iq0(Tab_Iq,psi%Basis)
   

  allocate (VQ(Ndim), VQQ(Ndim),VP(Ndim),VQP(Ndim),SQt(Ndim))
  allocate (CA(Ndim), CB(Ndim))
  VQ(:) = ZERO; VP(:) = ZERO;VQQ(:) = ZERO; VQP(:) = ZERO;SQt(:) = ONE

  call Calc_AVQ_SQ_nD_scd(psi, VQ, SQt, VQQ, Tab_Iq)
  call Calc_Av_imp_k_nD(psi, VP)
  call Calc_VQP_nD(VQP, psi)
  write (out_unit, *) 'VQQ = ', VQQ
  
   

  At(:) = CZERO 
  CB(:) = ZERO
  CA(:) = ZERO
   
   Do Ib = 1, Ndim
    !CA(Ib) = ONE/(TWO*(VQQ(Ib)-VQ(Ib)*VQ(Ib)))
      CA(Ib) = SQt(Ib)**2
    CB(Ib) = (VQP(Ib)-TWO*VP(Ib)*VQ(Ib))/(TWO*(VQQ(Ib)-VQ(Ib)*VQ(Ib))) 
    At(Ib) = complex(CA(Ib),-CB(Ib))   
   End do
   
   !At(:) = At(:)
   !write (out_unit, *) 'At = ', At
   
   IF (debug) THEN
   
      flush (out_unit)
      
   END IF

 deallocate(VQQ,VQP,CB,CA,VQ,VP,SQt)   
   
END SUBROUTINE


   SUBROUTINE Calc_Av_imp_k_1D(psi0, K, nio)
      USE QDUtil_m
      type(psi_t), intent(in), target                     :: psi0
      real(kind=Rkind), intent(inout)                     :: K
      integer, intent(in)                                 :: nio
      !locals variables---------------------------------------------------------
      type(psi_t), target                                 :: psi, ikpsi
      type(psi_t), target                                 :: psi_b, ikpsi_b
      logical, parameter                                  :: debug = .true.
      complex(kind=Rkind), pointer                        :: GB(:,:,:)
      complex(kind=Rkind), pointer                        :: d1gg(:,:,:)
      complex(kind=Rkind), pointer                        :: ikpsi0(:,:,:)
      Integer, allocatable                                :: Iq1(:), Iq2(:), Iq3(:)
      Integer, allocatable                                :: Ib1(:), Ib2(:), Ib3(:)
      integer                                             :: nq,inb,Ndim,Inbe,i1,i3
      !debuging----------------------------------------------------------------
      
      IF (debug) THEN
         flush (out_unit)
      END IF
      
      Ndim = size(psi0%Basis%tab_basis)
      call init_psi(psi, psi0%Basis, cplx=.true., grid=.true.)
      call init_psi(ikpsi, psi0%Basis, cplx=.true., grid=.true.)
      call init_psi(psi_b, psi0%Basis, cplx=.true., grid=.false.)
      call init_psi(ikpsi_b, psi0%Basis, cplx=.true., grid=.false.)
      Call Calc_index(Ib1=Ib1, Ib2=Ib2, Ib3=Ib3,Iq1=Iq1, Iq2=Iq2,&
       & Iq3=Iq3, Basis=psi0%Basis)
      psi%CVec = CZERO
      ikpsi%CVec = CZERO
      ikpsi_b%CVec = CZERO
      psi_b%CVec = CZERO
      IF (psi0%Grid) then
         psi%CVec(:) = psi0%CVec(:)
      ELSE
         call BasisTOGrid_nD_cplx(psi%CVec, psi0%CVec, psi0%Basis)
      END IF
       
        GB(1:Iq1(nio),1:Iq2(nio),1:Iq3(nio))    => psi%CVec
        ikpsi0(1:Iq1(nio),1:Iq2(nio),1:Iq3(nio))=> ikpsi%CVec
        d1gg(1:Iq2(nio),1:Iq2(nio),1:1)=> psi%Basis%tab_basis(nio)%d1gg(:, :, 1)
        
        DO i3 = 1, ubound(GB, dim=3)
          DO i1 = 1, ubound(GB, dim=1)
            ikpsi0(i1, :, i3) = ikpsi0(i1, :, i3)-EYE*matmul(d1gg(:,:,1),GB(i1,:,i3))
          END DO
        END DO
      call GridTOBasis_nD_cplx(psi_b%CVec, psi%CVec, psi0%Basis)
      call GridTOBasis_nD_cplx(ikpsi_b%CVec, ikpsi%CVec, psi0%Basis)
      K = dot_product(psi_b%CVec, ikpsi_b%CVec) 
      ! write (out_unit, *) '<psi/-id_x/psi> =', K
      IF (debug) THEN
         flush (out_unit)
      END IF
      
     call dealloc_psi(psi)
     call dealloc_psi(ikpsi)
     call dealloc_psi(ikpsi_b)
     call dealloc_psi(psi_b)
     Deallocate(Ib1, Ib2, Ib3,Iq1, Iq2, Iq3)
   END SUBROUTINE

   SUBROUTINE Calc_Av_imp_k_nD(psi0, K)
      USE QDUtil_m
      type(psi_t), intent(in)                                 :: psi0
      real(kind=Rkind), intent(inout)                         :: K(:)
      !locals variables---------------------------------------------------
      type(psi_t), target                                     :: psi
      logical, parameter                                      :: debug = .true.
      integer                                                 :: Ndim, Inb
      real(kind=Rkind)                                           :: p

      !debuging----------------------------------------------------------------

      IF (debug) THEN
         flush (out_unit)
      END IF
      Ndim = size(psi0%Basis%tab_basis) - 1
      call init_psi(psi, psi0%Basis, cplx=.TRUE., grid=.true.)

      psi%CVec(:) = CZERO
      K(:) = ZERO

      IF (psi0%Grid) then
         psi%CVec(:) = psi0%CVec(:)
      ELSE
         call BasisTOGrid_nD_cplx(psi%CVec, psi0%CVec, psi0%Basis)
      END IF

      Do Inb = 1, Ndim
         call Calc_Av_imp_k_1D(psi, K(inb), inb)
      End do

      write (out_unit, *) '<psi/-id_xi/psi> =', K

      IF (debug) THEN
         flush (out_unit)
      END IF
      CALL dealloc_psi(psi)

   END SUBROUTINE

 SUBROUTINE Calc_VQP_1D(psi0, VQP, Ib)
    USE QDUtil_m
    type(psi_t), intent(in), target              :: psi0
    real(kind=Rkind), intent(inout)              :: VQP
    integer, intent(in)                          :: Ib
    !locals variables----------------------      ---------------
    type(psi_t), target                           :: psi, d1psi
    logical, parameter                            :: debug = .true.
    complex(kind=Rkind), pointer                  :: GB(:,:,:)
    complex(kind=Rkind), pointer                  :: d1gg(:,:,:)
    complex(kind=Rkind), pointer                  :: d1psi0(:,:,:)
    Integer, allocatable                          :: Ib1(:), Ib2(:), Ib3(:)
    Integer, allocatable                          :: Iq1(:), Iq2(:), Iq3(:)
    integer                                       :: nq,inb,Ndim,Inbe,i1,i3,Iq

    complex(kind=Rkind), allocatable              :: psi_gb(:, :),d1psi_gb(:, :)
    logical                                       :: Endloop_q
    real(kind=Rkind), allocatable                 :: VQPEl(:)
    real(kind=Rkind)                              :: W, X
    real(kind=Rkind), allocatable                 :: N(:)
    integer, allocatable                          :: Tab_iq(:)
    !debuging--------------------------------------------------
    
    IF (debug) THEN
       flush (out_unit)
    END IF
    
    Ndim = size(psi0%Basis%tab_basis)
    call init_psi(psi,   psi0%Basis, cplx=.true., grid=.true.)
    call init_psi(d1psi, psi0%Basis, cplx=.true., grid=.true.)
    Call Calc_index(Ib1=Ib1, Ib2=Ib2, Ib3=Ib3,Iq1=Iq1, Iq2=Iq2,&
     & Iq3=Iq3, Basis=psi0%Basis)

    psi%CVec = CZERO
    d1psi%CVec = CZERO

    IF (psi0%Grid) then
       psi%CVec(:) = psi0%CVec(:)
    ELSE
       call BasisTOGrid_nD_cplx(psi%CVec, psi0%CVec, psi0%Basis)
    END IF
     
      GB(1:Iq1(Ib),1:Iq2(Ib),1:Iq3(Ib))    => psi%CVec
      d1psi0(1:Iq1(Ib),1:Iq2(Ib),1:Iq3(Ib))=> d1psi%CVec
      d1gg(1:Iq2(Ib),1:Iq2(Ib),1:1) => psi%Basis%tab_basis(Ib)%d1gg(:, :, 1)
      
      DO i3 = 1, ubound(GB, dim=3)
        DO i1 = 1, ubound(GB, dim=1)
          d1psi0(i1, :, i3) = d1psi0(i1, :, i3) + matmul(d1gg(:,:,1),GB(i1,:,i3))
        END DO
     END DO
     
     Allocate (psi_gb(psi%Basis%nq, psi%Basis%tab_basis(size(psi%Basis%tab_basis))%nb))
     Allocate (d1psi_gb(d1psi%Basis%nq, d1psi%Basis%tab_basis(size(d1psi%Basis%tab_basis))%nb))
     Allocate (VQPEl(psi%Basis%tab_basis(size(psi%Basis%tab_basis))%nb))
     Allocate (Tab_iq(size(Psi%Basis%tab_basis) - 1),N(ndim))

     psi_gb(:, :) =reshape(psi%CVec, shape=[psi%Basis%nq, psi%Basis%tab_basis(size(psi%Basis%tab_basis))%nb])
     d1psi_gb(:, :) =reshape(d1psi%CVec, shape=[d1psi%Basis%nq, d1psi%Basis%tab_basis(size(d1psi%Basis%tab_basis))%nb])
     X = ZERO
     N(:) = ZERO

    DO inbe = 1, psi%Basis%tab_basis(size(psi%Basis%tab_basis))%nb !electronic state
     VQPel(inbe) = ZERO
     Call Init_tab_ind(Tab_iq, psi%Basis%NDindexq)
      Iq = 0
      DO
         Iq = Iq + 1
         CALL increase_NDindex(Tab_iq, psi%Basis%NDindexq, Endloop_q)
         IF (Endloop_q) exit
         W = ONE
         DO inb = 1, size(psi%Basis%tab_basis) - 1
            W = W*psi%Basis%tab_basis(inb)%w(tab_iq(inb))
         END DO

         X = psi%Basis%tab_basis(ib)%x(tab_iq(ib))
         N(inbe) = N(inbe) + conjg(psi_gb(iq, inbe))*psi_gb(iq, inbe)*W
         VQPel(inbe) = VQPel(inbe)-EYE*conjg(psi_gb(iq, inbe))*W*(psi_gb(iq, inbe)+TWO*X*d1psi_gb(iq, inbe))

      END DO
    END DO
     VQP = sum(VQPel)/(Sum(N)**2)
    

      call dealloc_psi(psi)
      call dealloc_psi(d1psi)
      deallocate(Ib1, Ib2, Ib3,Iq1, Iq2, Iq3)
      deallocate(N,VQPEl,psi_gb,d1psi_gb,Tab_iq)

    IF (debug) THEN
       flush (out_unit)
    END IF
   
 END SUBROUTINE  


 SUBROUTINE Calc_VQP_nD(VQP,psi)
  TYPE(psi_t), intent(in)                        :: psi
  real(kind=Rkind) , intent(inout)               :: VQP(:)

  !Locals variabls ----------------------------------------------------------

  integer                                         :: Ib,ndim

    ndim  = size(psi%Basis%tab_basis) - 1

    Do Ib = 1,ndim
      call Calc_VQP_1D(psi, VQP(Ib), Ib)
    End do
    
    !write (out_unit, *) 'VQP',VQP

END SUBROUTINE

SUBROUTINE Calc_AVQ_SQ_nD(psi_in, AVQ, SQ,Tab_Iq)
   USE QDUtil_m
   TYPE(Psi_t), intent(in)                       :: psi_in
   integer,intent(in)                            :: Tab_Iq(:, :)
   real(kind=Rkind), intent(inout)               :: AVQ(:), SQ(:)

   logical, parameter                            :: debug = .false.
   TYPE(Psi_t)                                   :: psi
   complex(kind=Rkind), allocatable              :: psi_gb(:, :)
   real(kind=Rkind), allocatable                 :: AVQel(:, :), SQel(:, :),Q(:)
   real(kind=Rkind)                              :: W
   real(kind=Rkind), allocatable                 :: N(:)
   integer                                       :: iq, inbe, inb,nsurf,nq,ndim,Ib


   IF (debug) THEN
      write (out_unit, *) 'Beging AVQ'
      flush (out_unit)
   END IF

   ndim = size(psi_in%Basis%tab_basis) - 1
   nq =psi_in%Basis%nq
   nsurf=psi_in%Basis%tab_basis(ndim+1)%nb
   allocate (N(nsurf),AVQel(ndim, nsurf),SQel(ndim, nsurf),Q(ndim),psi_gb(nq, nsurf))
   CALL init_psi(psi, psi_in%Basis, cplx=.true., grid=.true.)

   IF (psi_in%Grid) then
      psi%CVec(:) = psi_in%CVec(:)
   ELSE
      CALL BasisTOGrid_nD_cplx(psi%CVec, psi_in%CVec, psi_in%Basis)
   END IF

   psi_gb(:, :) = reshape(psi%CVec, shape=[nq, nsurf])
   N(:) = ZERO
   AVQel(:, :) = ZERO
   SQel(:, :)  = ZERO

   DO inbe = 1, nsurf 

      DO Iq =  1,nq
         w = ONE
         Q(:)=ZERO

         DO inb = 1, ndim
            w = w*psi%Basis%tab_basis(inb)%w(Tab_Iq(inb,Iq))
            Q(inb)= psi%Basis%tab_basis(inb)%x(Tab_Iq(inb,Iq))!+psi%Basis%tab_basis(inb)%Q0
         END DO
         
           N(inbe) = N(inbe) + conjg(psi_gb(iq, inbe))*psi_gb(iq, inbe)*W
           AVQel(:, inbe) = AVQel(:, inbe) + conjg(psi_gb(iq, inbe))*(Q(:)*w)*psi_gb(iq, inbe)
           SQel(:, inbe) = SQel(:, inbe) + conjg(psi_gb(iq, inbe))*(Q(:)*Q(:)*w)*psi_gb(iq, inbe)

      END DO


   END DO

   DO inb = 1, ndim
      AVQ(inb) = sum(AVQel(inb,:))/(Sum(N)**2)
      SQ(inb) = sum(SQel(inb, :))/(Sum(N)**2)
      SQ(inb) = TWO*(SQ(inb) - AVQ(inb)*AVQ(inb))
      SQ(inb) = sqrt(ONE/SQ(inb))
   
   END DO   

   Deallocate (Psi_gb,N,AVQel,SQel,Q)
   CALL dealloc_psi(psi)
   IF (debug) THEN
      write (out_unit, *) 'END AVQ'
      flush (out_unit)
   END IF
END SUBROUTINE

SUBROUTINE Calc_AVQ_SQ_nD_scd(psi_in, AVQ, SQ,VQQ,Tab_Iq)
   USE QDUtil_m
   TYPE(Psi_t), intent(in)                       :: psi_in
   integer, intent(in)                           :: Tab_Iq(:, :)
   real(kind=Rkind), intent(inout)               :: AVQ(:), SQ(:),VQQ(:)

   logical, parameter                            :: debug = .false.
   TYPE(Psi_t)                                   :: psi
   complex(kind=Rkind), allocatable              :: psi_gb(:, :)
   real(kind=Rkind), allocatable                 :: AVQel(:, :), SQel(:, :),Q(:)
   real(kind=Rkind)                              :: W
   real(kind=Rkind), allocatable                 :: N(:)
   integer                                       :: iq, inbe, inb,nsurf,nq,ndim,Ib


   IF (debug) THEN
      write (out_unit, *) 'Beging AVQ'
      flush (out_unit)
   END IF

   ndim = size(psi_in%Basis%tab_basis) - 1
   nq =psi_in%Basis%nq
   nsurf=psi_in%Basis%tab_basis(ndim+1)%nb
   allocate (N(nsurf),AVQel(ndim, nsurf),SQel(ndim, nsurf),Q(ndim),psi_gb(nq, nsurf))
   CALL init_psi(psi, psi_in%Basis, cplx=.true., grid=.true.)

   IF (psi_in%Grid) then
      psi%CVec(:) = psi_in%CVec(:)
   ELSE
      CALL BasisTOGrid_nD_cplx(psi%CVec, psi_in%CVec, psi_in%Basis)
   END IF

   psi_gb(:, :) = reshape(psi%CVec, shape=[nq, nsurf])
   N(:) = ZERO
   AVQel(:, :) = ZERO
   SQel(:, :)  = ZERO

   DO inbe = 1, nsurf 

      DO Iq =  1,nq
         w = ONE
         Q(:)=ZERO

         DO inb = 1, ndim
            w = w*psi%Basis%tab_basis(inb)%w(Tab_Iq(inb,Iq))
            Q(inb)= psi%Basis%tab_basis(inb)%x(Tab_Iq(inb,Iq))
         END DO
         
           N(inbe) = N(inbe) + conjg(psi_gb(iq, inbe))*psi_gb(iq, inbe)*W
           AVQel(:, inbe) = AVQel(:, inbe) + conjg(psi_gb(iq, inbe))*(Q(:)*w)*psi_gb(iq, inbe)
           SQel(:, inbe) = SQel(:, inbe) + conjg(psi_gb(iq, inbe))*(Q(:)*Q(:)*w)*psi_gb(iq, inbe)

      END DO


   END DO

   DO inb = 1, ndim
      AVQ(inb) = sum(AVQel(inb,:))/(Sum(N)**2)
      VQQ(inb) = sum(SQel(inb,:))/(Sum(N)**2)
      SQ(inb) = sum(SQel(inb, :))/(Sum(N)**2)
      SQ(inb) = sqrt(SQ(inb) - AVQ(inb)*AVQ(inb))
      SQ(inb) = ONE/(SQ(inb)*sqrt(TWO))
   
   END DO   

   Deallocate (Psi_gb,N,AVQel,SQel,Q)
   CALL dealloc_psi(psi)
   IF (debug) THEN
      write (out_unit, *) 'END AVQ'
      flush (out_unit)
   END IF
END SUBROUTINE


end module
