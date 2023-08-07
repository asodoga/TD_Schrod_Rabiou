module Ana_psi_m
   USE Basis_m
   USE NDindex_m
   USE psi_m
   implicit none
   private
   public:: Population, Qpop,Calc_Integral_cplx
   public:: Calc_AVQ_nD0, Calc_AVQ_nD, Calc_Av_imp_k_nD

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
      AVQ = sum(AVQel)/(Sum(N)**2)
      SQ = sum(SQel)/(Sum(N)**2)
      SQ = sqrt(SQ - AVQ*AVQ)
      SQ = ONE/(SQ*sqrt(TWO))
      DO inbe = 1, Psi%Basis%tab_basis(size(Psi%Basis%tab_basis))%nb !electronic state
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

   SUBROUTINE Calc_AVQ_nD0(psi0, AVQ, AVQel, SQ, SQel)
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

   SUBROUTINE Calc_AVQ_nD(psi0, AVQ, SQ)
      USE QDUtil_m
      type(psi_t), intent(in), target                          :: psi0
      type(psi_t), target                                      :: psi
      real(kind=Rkind), intent(inout)                             :: AVQ(:), SQ(:)
      real(kind=Rkind), allocatable                               :: Q(:, :), W(:), Xe(:, :)
      real(kind=Rkind)                                            :: Norm
      logical, parameter                                       :: debug = .true.
      integer                                                  :: Inb, Ndim, Iq, inbe
      complex(kind=Rkind), pointer                                :: BB(:, :)
      real(kind=Rkind), allocatable                               :: Qte(:, :), SQe(:, :)

      IF (debug) THEN
         flush (out_unit)
      END IF

      Ndim = size(psi0%Basis%tab_basis) - 1
      call init_psi(psi, psi0%Basis, cplx=.TRUE., grid=.true.)
      call Calc_Q_grid(Q=Q, Basis=psi0%Basis, WnD=W)

      IF (psi0%Grid) then
         psi%CVec(:) = psi0%CVec(:)
      ELSE
         call BasisTOGrid_nD_cplx(psi%CVec, psi0%CVec, psi0%Basis)
      END IF
      allocate (Xe(Ndim, psi%Basis%tab_basis(Ndim + 1)%nb))
      SQ = ZERO
      Xe = ZERO
      AVQ = ZERO
      BB(1:psi%Basis%nq, 1:psi%Basis%tab_basis(Ndim + 1)%nb) => psi%CVec
      allocate (SQe(Ndim, psi%Basis%tab_basis(Ndim + 1)%nb), Qte(Ndim, psi%Basis%tab_basis(Ndim + 1)%nb))

      DO inbe = 1, psi%Basis%tab_basis(Ndim + 1)%nb
         DO Inb = 1, Ndim

            Qte(Inb, inbe) = dot_product(BB(:, inbe), W(:)*Q(:, Inb)*BB(:, inbe))

            Xe(Inb, inbe) = dot_product(BB(:, inbe), W(:)*Q(:, Inb)*Q(:, Inb)*BB(:, inbe))

            SQe(Inb, inbe) = sqrt(Xe(Inb, inbe) - Qte(Inb, inbe)*Qte(Inb, inbe))
            SQ(Inb) = ONE/(SQe(Inb, inbe)*sqrt(TWO))

         END DO
      End do

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



   SUBROUTINE Calc_Av_imp_k_1D(psi0, K, nio)
      USE QDUtil_m
      type(psi_t), intent(in), target                            :: psi0
      real(kind=Rkind), intent(inout)                            :: K
      integer, intent(in)                                        :: nio
      !locals variables---------------------------------------------------------
      type(psi_t), target                                     :: psi, ikpsi
      type(psi_t), target                                     :: psi_b, ikpsi_b
      logical, parameter                                      :: debug = .true.
      complex(kind=Rkind), pointer                               :: GB(:,:,:)
      complex(kind=Rkind), pointer                                  ::d1gg_cplx(:,:,:)
      complex(kind=Rkind), pointer                               :: ikpsi0(:,:,:)
      real(kind=Rkind)                                           :: p
      Integer, allocatable         :: Ib1(:), Ib2(:), Ib3(:),Iq1(:), Iq2(:), Iq3(:)
      integer                                                 :: nq,inb,Ndim,Inbe,i1,i3
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
        d1gg_cplx(1:Iq2(nio),1:Iq2(nio),1:1)=> psi%Basis%tab_basis(nio)%d1gg_cplx(:, :, 1)
        
        DO i3 = 1, ubound(GB, dim=3)
          DO i1 = 1, ubound(GB, dim=1)
            ikpsi0(i1, :, i3) = ikpsi0(i1, :, i3)-EYE*matmul(d1gg_cplx(:,:,1),GB(i1,:,i3))
          END DO
       END DO
      call GridTOBasis_nD_cplx(psi_b%CVec, psi%CVec, psi0%Basis)
      call GridTOBasis_nD_cplx(ikpsi_b%CVec, ikpsi%CVec, psi0%Basis)
      p = psi%Basis%tab_basis(nio)%Imp_k
      K = dot_product(psi_b%CVec, ikpsi_b%CVec) !+ p
      ! write (out_unit, *) '<psi/-id_x/psi> =', K
      IF (debug) THEN
         flush (out_unit)
      END IF
      
      CALL dealloc_psi(psi)
      CALL dealloc_psi(ikpsi)
      CALL dealloc_psi(ikpsi_b)
      CALL dealloc_psi(psi_b)
      Deallocate(Ib1, Ib2, Ib3,Iq1, Iq2, Iq3)
   END SUBROUTINE

   SUBROUTINE Calc_Av_imp_k_nD(psi0, K)
      USE QDUtil_m
      type(psi_t), intent(in)                                 :: psi0
      real(kind=Rkind), intent(inout)                            :: K(:)
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



   SUBROUTINE dnpsi_cplx(dnpsi_g, Psi_g, nio,Basis)
    USE QDUtil_m
    USE Basis_m
    TYPE(Basis_t), intent(in), target            :: Basis
    complex(kind=Rkind), intent(in), target         :: psi_g(:)
    complex(kind=Rkind), intent(inout), target      :: dnpsi_g(:)
    integer,intent(in)                           :: nio
    
    !locals variables=============================================================================
    
    complex(kind=Rkind), pointer                    :: psi_ggb(:, :, :)
    complex(kind=Rkind), pointer                       :: dngg(:, :)
    complex(kind=Rkind), pointer                    :: dnpsi_ggb(:, :, :)
    logical, parameter                           :: debug = .true.
    integer                                      :: iq, i1, i3, inb, Ndim
    integer, allocatable                         :: Iq1(:), Iq2(:), Iq3(:), Ib1(:), Ib2(:), Ib3(:)
    
    IF (debug) THEN
       write(out_unit,*) 'BEGINNING d2psi'
       flush (out_unit)
    END IF
    
    Ndim = size(Basis%tab_basis)
    call Calc_index(Ib1, Ib2, Ib3, Iq1, Iq2, Iq3, Basis)
    
    dnpsi_g(:) = CZERO
    DO inb = 1, Ndim - 1
    
       dnpsi_ggb(1:Iq1(inb), 1:Iq2(inb), 1:Iq3(inb)) => dnpsi_g
       psi_ggb(1:Iq1(inb), 1:Iq2(inb), 1:Iq3(inb))   => psi_g
       
       if(nio==1) then
           dngg(1:Iq2(inb), 1:Iq2(inb)) => Basis%tab_basis(inb)%d1gg_cplx
       else
           dngg(1:Iq2(inb), 1:Iq2(inb)) => Basis%tab_basis(inb)%d2gg_cplx
       end if
       
       DO i3 = 1, ubound(psi_ggb, dim=3)
          DO i1 = 1, ubound(psi_ggb, dim=1)
             dnpsi_ggb(i1, :, i3) = dnpsi_ggb(i1, :, i3)+ matmul(dngg, psi_ggb(i1, :, i3))
          END DO
       END DO
    END DO
    Deallocate (Ib1, Ib2, Ib3, Iq1, Iq2, Iq3)
    IF (debug) THEN
       write(out_unit,*) 'END d2psi'
       flush (out_unit)
    END IF
 END SUBROUTINE 



  SUBROUTINE Calc_Integral_cplx(psi_in, Alpha, nio)
    USE QDUtil_m
    TYPE(Psi_t), intent(in)                    :: psi_in
    integer, intent(in)                        :: nio
    complex(kind=Rkind)  ,intent(inout)           :: Alpha(:)
    
    !locals variables========================================================================================
    
    TYPE(Psi_t)                                :: psi,dnpsi
    complex(kind=Rkind), allocatable              :: psi_gb(:, :)
    complex(kind=Rkind), allocatable              :: dnpsi_gb(:, :)
    logical                                    :: Endloop_q
    logical, parameter                         :: debug = .false.
    complex(kind=Rkind), allocatable              :: Int1(:),Int2(:),Int1el(:),Int2el(:),Int0(:)
  
    real(kind=Rkind)                              :: W
    integer, allocatable                       :: Tab_iq(:)
    integer                                    :: iq, inbe, inb,ndim
    
    !debunging & allocation=============================================================================================
    
    if (debug) then
       write (out_unit, *) 'Beging Calc_Integral_cplx'
       flush (out_unit)
    end if
    ndim = size(psi_in%Basis%tab_basis)-1
    
    allocate (Int1(ndim))
    allocate (Int2(ndim))
    allocate (Int0(ndim))
    
    allocate (Int1el(psi_in%Basis%tab_basis(ndim+1)%nb))
    allocate (Int2el(psi_in%Basis%tab_basis(ndim+1)%nb))
    allocate (psi_gb(psi_in%Basis%nq, psi_in%Basis%tab_basis(ndim+1)%nb))
    !sprint*,'cc',size(psi_gb,dim=1),size(psi_gb,dim=2)
    !STOP 'cc'
    allocate (dnpsi_gb(psi_in%Basis%nq, psi_in%Basis%tab_basis(ndim+1)%nb))
    allocate (Tab_iq(ndim))
    
    call init_psi(psi, psi_in%Basis, cplx=.true., grid=.true.)
    call init_psi(dnpsi, psi_in%Basis, cplx=.true., grid=.true.)
    
    if (psi_in%Grid) then
       psi%CVec(:) = psi_in%CVec(:)
    else
       call BasisTOGrid_nD_cplx(psi%CVec, psi_in%CVec, psi_in%Basis)
    end if
    
    call dnpsi_cplx(dnpsi%CVec, psi%CVec, nio,psi%Basis)
    psi_gb(:, :)    =  reshape(psi%CVec,shape=[psi%Basis%nq, psi%Basis%tab_basis(ndim+1)%nb])
    dnpsi_gb(:, :) =   reshape(dnpsi%CVec,shape=[dnpsi%Basis%nq, dnpsi%Basis%tab_basis(ndim+1)%nb])
    
    !Evaluation of Alpha=================================================================================================    
    
    do inb= 1,ndim
    
    DO inbe = 1, psi%Basis%tab_basis(ndim+1)%nb !electronic state
       Int1el(inbe) = CZERO
       Int2el(inbe) = CZERO
       Call Init_tab_ind(Tab_iq, psi%Basis%NDindexq)
       
       Iq = 0
       DO
          Iq = Iq + 1
          call increase_NDindex(Tab_iq, psi%Basis%NDindexq, Endloop_q)
          if (Endloop_q) exit
         
             W = psi%Basis%tab_basis(inb)%w(tab_iq(inb))
             Int1el(inbe) = Int1el(inbe) + psi_gb(iq, inbe)*W*psi_gb(iq, inbe)
             Int2el(inbe) = Int2el(inbe) + psi_gb(iq, inbe)*W*dnpsi_gb(iq, inbe)
          
       END DO
    END DO
    
    Int1(inb) = sum(Int1el)
    Int2(inb) = sum(Int2el)
    Int0(inb)  = Int2(inb)/Int1(inb)
    Alpha(inb) = EYE*Int0(inb)
    
  end do ! basis  
   
    print*,  'Alpha=', Alpha
    
    Deallocate (Tab_iq)
    Deallocate (psi_gb,dnpsi_gb,Int0,Int1,Int2,Int1el,Int2el)
    call dealloc_psi(psi)
    call dealloc_psi(dnpsi)
    
    if (debug) then
       write (out_unit, *) 'END Calc_Integral_cplx'
       flush (out_unit)
    end if
 END SUBROUTINE

end module
