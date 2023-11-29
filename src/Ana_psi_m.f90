module Ana_psi_m
   USE Basis_m
   USE NDindex_m
   USE psi_m
   implicit none
   private
   public:: Population, Qpop
   public::  Calc_AVQ_nD, Calc_Av_imp_k_nD
   public :: Calc_Avg_A_nD

contains

   SUBROUTINE Calc_AVQ_1D(psi_in, AVQ1, AVQel1, SQ1, SQel1, ib)
      USE QDUtil_m
      logical, parameter                            :: debug = .false.
      TYPE(Psi_t), intent(in)                       :: psi_in
      TYPE(Psi_t)                                   :: psi
      complex(kind=Rkind), allocatable              :: psi_gb(:, :)
      logical                                       :: Endloop_q
      real(kind=Rkind), intent(inout), optional     :: AVQ1, SQ1
      real(kind=Rkind), intent(inout), optional     :: AVQel1(:), SQel1(:)

      real(kind=Rkind), allocatable                 :: AVQel(:), SQel(:)
      real(kind=Rkind)                              :: AVQ, SQ
      integer, intent(in)                           :: ib
      real(kind=Rkind)                              :: WnD, X
      real(kind=Rkind), allocatable                 :: N(:)
      integer, allocatable                          :: Tab_iq(:)
      integer                                       :: iq, inbe, inb
      IF (debug) THEN
         write (out_unit, *) 'Beging AVQ'
         flush (out_unit)
      END IF
      allocate (N(psi_in%Basis%tab_basis(size(psi_in%Basis%tab_basis))%nb))
      allocate (AVQel(psi_in%Basis%tab_basis(size(psi_in%Basis%tab_basis))%nb))
      allocate (SQel(psi_in%Basis%tab_basis(size(psi_in%Basis%tab_basis))%nb))
      CALL init_psi(psi, psi_in%Basis, cplx=.TRUE., grid=.true.)

      IF (psi_in%Grid) then
         psi%CVec(:) = psi_in%CVec(:)
      ELSE
         CALL BasisTOGrid_nD_cplx(psi%CVec, psi_in%CVec, psi_in%Basis)
      END IF

      Allocate (Psi_gb(psi%Basis%nq, psi%Basis%tab_basis(size(psi%Basis%tab_basis))%nb))
      Allocate (Tab_iq(size(Psi%Basis%tab_basis) - 1))
      psi_gb(:, :) = reshape(psi%CVec, shape=[psi%Basis%nq, psi%Basis%tab_basis(size(psi%Basis%tab_basis))%nb])

      X = ZERO
      N(:) = ZERO

      DO inbe = 1, psi%Basis%tab_basis(size(psi%Basis%tab_basis))%nb !electronic state

         AVQel(inbe) = ZERO
         Call Init_tab_ind(Tab_iq, psi%Basis%NDindexq)
         Iq = 0
         DO
            Iq = Iq + 1
            CALL increase_NDindex(Tab_iq, psi%Basis%NDindexq, Endloop_q)
            IF (Endloop_q) exit
            WnD = ONE
            DO inb = 1, size(psi%Basis%tab_basis) - 1
               WnD = WnD*psi%Basis%tab_basis(inb)%w(tab_iq(inb))
            END DO
            X = psi%Basis%tab_basis(ib)%x(tab_iq(ib))
            N(inbe) = N(inbe) + conjg(psi_gb(iq, inbe))*psi_gb(iq, inbe)*WnD
            AVQel(inbe) = AVQel(inbe) + conjg(psi_gb(iq, inbe))*(X*WnD)*psi_gb(iq, inbe)
            SQel(inbe) = SQel(inbe) + conjg(psi_gb(iq, inbe))*(X*X*WnD)*psi_gb(iq, inbe)
         END DO
      END DO
       !write (out_unit, *) 'Ib,VQQ',Ib,SQel,Sum(N)**2
      AVQ = sum(AVQel)/(Sum(N)**2)
      SQ = sum(SQel)/(Sum(N)**2)
      SQ = sqrt(SQ - AVQ*AVQ)
      SQ = ONE/(SQ*sqrt(TWO))
      DO inbe = 1, psi%Basis%tab_basis(size(Psi%Basis%tab_basis))%nb !electronic state
         if (N(inbe) /= ZERO) AVQel(inbe) = AVQel(inbe)/(N(inbe)**2)
      END DO
      if (present(AVQ1)) AVQ1 = AVQ
      if (present(AVQel1)) AVQel1 = AVQel
      if (present(SQ1)) SQ1 = SQ
      if (present(AVQel1)) SQel1 = SQel

      ! if(present(AVQ1))  print*,   'AVQ1=',    AVQ1
      ! if(present(AVQel1)) print*,  'AVQel1=',  AVQel1
      ! if(present(SQ1))    print*,  'SQ1=',     SQ1

      Deallocate (Tab_iq)
      Deallocate (Psi_gb)
      CALL dealloc_psi(psi)
      IF (debug) THEN
         write (out_unit, *) 'END AVQ'
         flush (out_unit)
      END IF
   END SUBROUTINE

   SUBROUTINE Calc_AVQ_nD(psi0, AVQ, AVQel, SQ, SQel)
      USE QDUtil_m
      type(psi_t), intent(in), target                         :: psi0
      type(psi_t), target                                     :: psi
      real(kind=Rkind), intent(inout), optional                  :: AVQ(:), SQ(:)
      real(kind=Rkind), intent(inout), optional                  :: AVQel(:, :), SQel(:, :)
      logical, parameter                                      :: debug = .true.
      integer                                                 :: Inb, Ndim

      IF (debug) THEN
         flush (out_unit)
      END IF

      Ndim = size(psi0%Basis%tab_basis)
      call init_psi(psi, psi0%Basis, cplx=.TRUE., grid=.true.)

      IF (psi0%Grid) then
         psi%CVec(:) = psi0%CVec(:)
      ELSE
         call BasisTOGrid_nD_cplx(psi%CVec, psi0%CVec, psi0%Basis)
      END IF

      DO Inb = 1, Ndim - 1

         if (present(AVQel)) call Calc_AVQ_1D(psi_in=psi, AVQel1=AVQel(Inb, :), ib=Inb)
         if (present(AVQ)) call Calc_AVQ_1D(psi_in=psi, AVQ1=AVQ(Inb), ib=Inb)

         if (present(SQel)) call Calc_AVQ_1D(psi_in=psi, SQel1=SQel(Inb, :), ib=Inb)
         if (present(SQ)) call Calc_AVQ_1D(psi_in=psi, SQ1=SQ(Inb), ib=Inb)

      END DO

      write (out_unit, *) '<psi/Q/psi> =', AVQ
      write (out_unit, *) 'SQ =', SQ
      IF (debug) THEN
         flush (out_unit)
      END IF
      CALL dealloc_psi(psi)
   END SUBROUTINE


   subroutine Population(Psi, Pop)
      implicit none
      type(Psi_t), intent(in), target                 :: Psi
      complex(kind=Rkind), pointer                       :: Psi_bb(:, :)
      real(Kind=Rkind), intent(inout), allocatable       ::Pop(:)
      integer                                         :: inb
      real(Kind=Rkind)                                   :: Norm

      Psi_bb(1:Psi%Basis%nb, 1:Psi%Basis%tab_basis(size(Psi%Basis%tab_basis))%nb) => Psi%CVec
      call Calc_Norm_OF_Psi(Psi, Norm)

      do inb = 1, Psi%Basis%tab_basis(size(Psi%Basis%tab_basis))%nb
         Pop(inb) = dot_product(Psi_bb(:, inb), Psi_bb(:, inb))/Norm
      end do

   end subroutine Population

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
     logical, parameter                                      :: debug = .true.
     integer                                                 :: Ndim, Ib
     
     !debuging----------------------------------------------------------------
     
     IF (debug) THEN
     
        flush (out_unit)
        
     END IF
     
     Ndim = size(psi%Basis%tab_basis) - 1 

    allocate (VQ(Ndim), VQQ(Ndim),VP(Ndim),VQP(Ndim),SQt(Ndim))
    allocate (CA(Ndim), CB(Ndim))
    VQ(:) = ZERO; VP(:) = ZERO;VQQ(:) = ZERO; VQP(:) = ZERO;SQt(:) = ONE

    call Calc_AVQ_nD(psi0=psi, AVQ=VQ, SQ=SQt)
    call Calc_Av_imp_k_nD(psi,VP)
    call Calc_VQQ_nD(VQQ,psi) 
    call Calc_VQP_nD(VQP,psi)

    At(:) = CZERO 
    CB(:) = ZERO
    CA(:) = ZERO
     
     Do Ib = 1, Ndim
      CA(Ib) = ONE/(FOUR*(VQQ(Ib)-VQ(Ib)*VQ(Ib))) 
      CB(Ib) = (VQP(Ib)-TWO*VP(Ib)*VQ(Ib))/(FOUR*(VQQ(Ib)-VQ(Ib)*VQ(Ib))) 
      At(Ib) = complex(CA(Ib),CB(Ib))   
     End do
     
     At(:) = TWO*At(:)
     write (out_unit, *) 'At = ', At
     
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
      call init_psi(psi, psi0%Basis, cplx=.TRUE., grid=.true.)
      call init_psi(ikpsi, psi0%Basis, cplx=.TRUE., grid=.true.)
      call init_psi(psi_b, psi0%Basis, cplx=.TRUE., grid=.false.)
      call init_psi(ikpsi_b, psi0%Basis, cplx=.TRUE., grid=.false.)
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


    SUBROUTINE Calc_VQQ_1D(VQQ,psi_in,Ib)
      USE QDUtil_m
      logical, parameter                            :: debug = .false.
      TYPE(Psi_t), intent(in)                       :: psi_in
      TYPE(Psi_t)                                   :: psi
      complex(kind=Rkind), allocatable              :: psi_gb(:, :)
      logical                                       :: Endloop_q
      real(kind=Rkind), intent(inout)               :: VQQ
      real(kind=Rkind), allocatable                 :: VQQel(:)
      integer, intent(in)                           :: ib
      real(kind=Rkind)                              :: WnD, X
      real(kind=Rkind), allocatable                 :: N(:)
      integer, allocatable                          :: Tab_iq(:)
      integer                                       :: iq, inbe, inb

      IF (debug) THEN
         write (out_unit, *) 'Beging AVQ'
         flush (out_unit)
      END IF

      allocate (N(psi_in%Basis%tab_basis(size(psi_in%Basis%tab_basis))%nb))
      allocate (VQQel(psi_in%Basis%tab_basis(size(psi_in%Basis%tab_basis))%nb))
      CALL init_psi(psi, psi_in%Basis, cplx=.true., grid=.true.)

      IF (psi_in%Grid) then
         psi%CVec(:) = psi_in%CVec(:)
      ELSE
         CALL BasisTOGrid_nD_cplx(psi%CVec, psi_in%CVec, psi_in%Basis)
      END IF
      Allocate (Psi_gb(psi%Basis%nq, psi%Basis%tab_basis(size(psi%Basis%tab_basis))%nb))
      Allocate (Tab_iq(size(Psi%Basis%tab_basis) - 1))
      psi_gb(:, :) =reshape(psi%CVec, shape=[psi%Basis%nq, psi%Basis%tab_basis(size(psi%Basis%tab_basis))%nb])
      X = ZERO
      N(:) = ZERO
      DO inbe = 1, psi%Basis%tab_basis(size(psi%Basis%tab_basis))%nb !electronic state
         VQQel(inbe) = ZERO
         Call Init_tab_ind(Tab_iq, psi%Basis%NDindexq)
         Iq = 0
         DO
            Iq = Iq + 1
            CALL increase_NDindex(Tab_iq, psi%Basis%NDindexq, Endloop_q)
            IF (Endloop_q) exit
            WnD = ONE
            DO inb = 1, size(psi%Basis%tab_basis) - 1
               WnD = WnD*psi%Basis%tab_basis(inb)%w(tab_iq(inb))
            END DO
            X = psi%Basis%tab_basis(ib)%x(tab_iq(ib))
            N(inbe) = N(inbe) + conjg(psi_gb(iq, inbe))*psi_gb(iq, inbe)*WnD
            VQQel(inbe) = VQQel(inbe) + conjg(psi_gb(iq, inbe))*(X*X*WnD)*psi_gb(iq, inbe)
         END DO
      END DO
      VQQ = sum(VQQel)/(Sum(N)**2)
   
       DO inbe = 1, psi%Basis%tab_basis(size(Psi%Basis%tab_basis))%nb !electronic state
         if (N(inbe) /= ZERO) VQQel(inbe) = VQQel(inbe)/(N(inbe)**2)
      END DO

      Deallocate (Tab_iq,N)
      Deallocate (psi_gb)
      CALL dealloc_psi(psi)

      IF (debug) THEN
         write (out_unit, *) 'END AVQ'
         flush (out_unit)
      END IF

   END SUBROUTINE


   SUBROUTINE Calc_VQQ_nD(VQQ,psi)
      TYPE(psi_t), intent(in)                                :: psi
      real(kind=Rkind) , intent(inout)                       :: VQQ(:)
     
      !Locals variabls ----------------------------------------------------------
     
      integer                                                :: Ib,ndim
     
        ndim  = size(psi%Basis%tab_basis) - 1
     
     
        Do Ib = 1,ndim
          call Calc_VQQ_1D(VQQ(Ib),psi,Ib)
        End do

        
       ! write (out_unit, *) 'VQQ',VQQ

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


end module
