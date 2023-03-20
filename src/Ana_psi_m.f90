 module Ana_psi_m
    USE Basis_m
    USE NDindex_m
    USE psi_m
    implicit none
     private
     public:: Population,Average_Q,Qpop,Calc_AVQ_1D,Calc_std_dev_AVQ_1D

     contains

     subroutine Population(Psi,Pop)
            implicit none
            type (Psi_t) ,intent(in)              ,target        :: Psi
            complex (kind=Rk), pointer                           :: Psi_bb(:,:)
            real(Kind = Rk), intent(inout) ,allocatable          ::Pop(:)
            integer                                              :: inb
            real(Kind = Rk)                                      :: Norm


            Psi_bb(1:Psi%Basis%nb,1:Psi%Basis%tab_basis(size(Psi%Basis%tab_basis))%nb) => Psi%CVec
            call Calc_Norm_OF_Psi(Psi,Norm)

            do inb = 1,Psi%Basis%tab_basis(size(Psi%Basis%tab_basis))%nb
             Pop(inb) = dot_product(Psi_bb(:,inb),Psi_bb(:,inb))/Norm
            end do
          ! write(*,*) Pop

         end subroutine Population


         SUBROUTINE Calc_std_dev_AVQ_1D(Psi_in,i_1D,AVQ,SQ)
             USE UtilLib_m
             logical,         parameter      :: debug = .false.
             TYPE(Psi_t), intent(in)         :: Psi_in
             integer,intent(in)              :: i_1D ! indicating  the dimension
             TYPE(Psi_t)                     :: Psi
             real(kind=Rk),intent(inout)     :: AVQ,SQ
             real(kind=Rk)                   :: WnD,X,STD_DQ,Norm
             integer                         :: Iq

              call Calc_Norm_OF_Psi(Psi_in,Norm)
             ! write(out_unitp,*) 'Norm=',Norm
             IF (debug) THEN
                ! write(out_unitp,*) 'Beging AVQ,STD_DQ'
                 flush(out_unitp)
             END IF
             CALL init_psi(psi,   Psi_in%Basis,    cplx=.TRUE.   ,grid =.true.)
             IF(Psi_in%Grid) then
                 psi%CVec(:)= psi_in%CVec(:)
             ELSE
                 CALL BasisTOGrid_nD_cplx(Psi%CVec,Psi_in%CVec,Psi_in%Basis)
             END IF
                 AVQ     = ZERO
                 STD_DQ  = ZERO
                 X       =ZERO
              do Iq = 1,psi%Basis%tab_basis(i_1D)%nq
                  AVQ =  AVQ+conjg(psi%CVec(Iq))*psi%Basis%tab_basis(i_1D)%x(Iq)*psi%CVec(Iq)*psi%Basis%tab_basis(i_1D)%w(Iq)
                  X   =  X+conjg(psi%CVec(Iq))*(psi%Basis%tab_basis(i_1D)%x(Iq)**2)*psi%CVec(Iq)*psi%Basis%tab_basis(i_1D)%w(Iq)

              end do

             STD_DQ =sqrt(X-AVQ*AVQ)
             SQ = ONE/(STD_DQ*sqrt(TWO))
             if( (x-real(int(x),kind=Rk) )<= ONETENTH**8)   x = real(int(x),kind=Rk)
             if( (SQ-real(int(SQ),kind=Rk) )<= ONETENTH**8)   SQ = real(int(SQ),kind=Rk)
               Print *,"<psi|Q|psi> = ",AVQ ,"<psi|Q**2|psi> = ",X,"sqrt(<psi|Q**2|psi> - <psi|Q|psi> )= ", STD_DQ,'SQ=',SQ
             CALL dealloc_psi(psi)
             IF (debug) THEN
                ! write(out_unitp,*) 'END AVQ,STD_DQ'
                 flush(out_unitp)
             END IF

         END SUBROUTINE Calc_std_dev_AVQ_1D


         SUBROUTINE Calc_sx_avq(psi_in,Q,sQ)
             USE UtilLib_m
             real(kind=Rk)  ,intent(inout)   :: Q(:),sQ(:)
             TYPE(Psi_t), intent(in)         :: Psi_in
             integer                         :: ib
             logical,         parameter      :: debug = .false.



             IF (debug) THEN
                 write(out_unitp,*) 'Beging evaluation of sx and <q>'
                 flush(out_unitp)
             END IF

               Do ib = 1,size(Q)
                  CALL  Calc_std_dev_AVQ_1D(psi_in,ib,Q(ib),sQ(ib))
               End Do

             IF (debug) THEN
                 write(out_unitp,*) 'End evaluation of sx and <q>'
                 flush(out_unitp)
             END IF

         END SUBROUTINE  Calc_sx_avq






         SUBROUTINE Calc_AVQ_1D(Psi_in,i_dim, AVQ,AVQel)

             USE UtilLib_m
             logical,         parameter      :: debug = .false.
             TYPE(Psi_t), intent(in)         :: Psi_in
             integer , intent(in)            :: i_dim
             TYPE(Psi_t)                     :: Psi
             complex (kind=Rk),allocatable   :: Psi_gb(:,:)
             logical                         :: Endloop_q
             real(kind=Rk),intent(inout)     :: AVQ
             real(kind=Rk),intent(inout)     :: AVQel(:)
             real(kind=Rk)                   :: WnD,X
             real(kind=Rk),allocatable       :: N(:)
             integer,        allocatable     :: Tab_iq(:)
             integer                         :: iq,inbe,inb

             IF (debug) THEN
                 write(out_unitp,*) 'Beging AVQ'
                 flush(out_unitp)
             END IF
             allocate(N(Psi_in%Basis%tab_basis(size(Psi_in%Basis%tab_basis))%nb))
             CALL init_psi(psi,   Psi_in%Basis,    cplx=.TRUE.   ,grid =.true.)
             IF(Psi_in%Grid) then
                 psi%CVec(:)= psi_in%CVec(:)
             ELSE
                 CALL BasisTOGrid_nD_cplx(Psi%CVec,Psi_in%CVec,Psi_in%Basis)
             END IF
             Allocate(Psi_gb(Psi%Basis%nq,Psi%Basis%tab_basis(size(Psi%Basis%tab_basis))%nb))
             Allocate(Tab_iq(size(Psi%Basis%tab_basis)-1))
             Psi_gb(:,:) = reshape(Psi%CVec,shape= [psi%Basis%nq,Psi%Basis%tab_basis(size(Psi%Basis%tab_basis))%nb])
             X = ZERO
             N(:) = ZERO
             DO inbe = 1,Psi%Basis%tab_basis(size(Psi%Basis%tab_basis))%nb !electronic state
                 AVQel(inbe) = ZERO

                 Call Init_tab_ind(Tab_iq,Psi%Basis%NDindexq)
                 Iq = 0
                 DO
                     Iq = Iq+1
                     CALL increase_NDindex(Tab_iq,Psi%Basis%NDindexq,Endloop_q)
                     IF (Endloop_q) exit
                     WnD= ONE
                     DO inb = 1,size(psi%Basis%tab_basis)-1
                         WnD  =WnD* Psi%Basis%tab_basis(inb)%w(tab_iq(inb))
                     END DO
                     X = Psi%Basis%tab_basis(i_dim)%x(tab_iq(i_dim))*WnD
                     N(inbe) = N(inbe) + conjg(Psi_gb(iq,inbe))*Psi_gb(iq,inbe)*WnD
                     AVQel(inbe) = AVQel(inbe) + conjg(Psi_gb(iq,inbe))*Psi_gb(iq,inbe)*X
                 END DO

             END DO
             AVQ = sum(AVQel)

             DO inbe = 1,Psi%Basis%tab_basis(size(Psi%Basis%tab_basis))%nb !electronic state
                 if(N(inbe) /= ZERO)  AVQel(inbe) = AVQel(inbe)/N(inbe)
             END DO

             Deallocate(Tab_iq)
             Deallocate(Psi_gb)
             CALL dealloc_psi(psi)
             IF (debug) THEN
                 write(out_unitp,*) 'END AVQ'
                 flush(out_unitp)
             END IF

         END SUBROUTINE Calc_AVQ_1D





        SUBROUTINE Average_Q(Psi,Qm)
      USE Basis_m
      USE UtilLib_m
      type(Psi_t)      , intent(in) ,target        :: Psi
       type(Psi_t)     ,target                     :: Psi_g
      complex (kind=Rk), pointer                   :: Psi_ggb(:,:,:)
      real (kind=Rk),    pointer                   :: Q_g(:),W_g(:)
      real (kind=Rk),intent(inout)                 :: Qm(:)
      logical,           parameter                 :: debug = .true.
      integer                                      :: iq,i1,i3,inb,ndim
      integer , allocatable                        :: Iq1(:),Iq2(:),Iq3(:),Ib1(:),Ib2(:),Ib3(:)

      IF (debug) THEN
        !write(out_unitp,*) 'BEGINNING Average_Q'
        flush(out_unitp)
      END IF
      Ndim = size(Psi%Basis%tab_basis)
      call Calc_iqib( Ib1,Ib2,Ib3,Iq1,Iq2,Iq3,Ndim,Psi%Basis)
       Qm(:) = ZERO
       if(Psi%Grid) then
            DO inb = 1,ndim-1
              Psi_ggb(1:Iq1(inb),1:Iq2(inb),1:iq3(inb))  => Psi%CVec
              Q_g(1:Psi%Basis%tab_basis(inb)%nq) => Psi%Basis%tab_basis(inb)%x
              W_g(1:Psi%Basis%tab_basis(inb)%nq) => Psi%Basis%tab_basis(inb)%w

              DO i3=1,ubound(Psi_ggb,dim=3)
                 DO i1=1,ubound(Psi_ggb,dim=1)
                    Qm(inb) =  Qm(inb) + dot_product(Psi_ggb(i1,:,i3), W_g(:)*Q_g(:)*Psi_ggb(i1,:,i3))
                 END DO
              END DO
            END DO
       else
           CALL init_psi(psi_g,   psi%Basis,    cplx=.TRUE.   ,grid =.true.)
           Psi_g%CVec(:)= CZERO
             call BasisTOGrid_nD_cplx(Psi_g%CVec,Psi%CVec,Psi%Basis)
             DO inb = 1,ndim-1
                 Psi_ggb(1:Iq1(inb),1:Iq2(inb),1:iq3(inb))  => Psi_g%CVec
                 Q_g(1:Psi%Basis%tab_basis(inb)%nq) => Psi%Basis%tab_basis(inb)%x
                 W_g(1:Psi%Basis%tab_basis(inb)%nq) => Psi%Basis%tab_basis(inb)%w

                 DO i3=1,ubound(Psi_ggb,dim=3)
                     DO i1=1,ubound(Psi_ggb,dim=1)
                         Qm(inb) =  Qm(inb) + dot_product(Psi_ggb(i1,:,i3), W_g(:)*Q_g(:)*Psi_ggb(i1,:,i3))
                     END DO
                 END DO
             END DO
             CALL dealloc_psi(Psi_g)
           deallocate(Iq1,Iq2,Iq3,Ib1,Ib2,Ib3)
       end if
      IF (debug) THEN
       !	write(out_unitp,*) 'END Average_Q
       	flush(out_unitp)
      END IF
  END SUBROUTINE  Average_Q

         SUBROUTINE Qpop(Psi,Qp)
                 USE Basis_m
                 USE UtilLib_m
                 type(Psi_t)      , intent(in) ,target         :: Psi
                 type(Psi_t)     ,target                       :: Psi_g
                 complex(kind=Rk), pointer                     :: psi_gb(:,:)

                 real (kind=Rk),intent(inout)                  :: Qp(:)
                 real (kind=Rk)                                :: Norm(2)
                 logical,           parameter                  :: debug = .true.
                 integer                                       :: iq,inb,ndim

                 IF (debug) THEN
                    !write(out_unitp,*) 'BEGINNING Qpop'
                    flush(out_unitp)
                 END IF
                 Ndim = size(Psi%Basis%tab_basis)
                      call  init_psi(psi_g,   psi%Basis,    cplx=.TRUE.   ,grid =.true.)
                          Psi_g%CVec(:)= CZERO
                     call  BasisTOGrid_nD_cplx(Psi_g%CVec,Psi%CVec,Psi%Basis)
                   do inb = 1,Psi%Basis%tab_basis(2)%nb
                     psi_gb(1:Psi%Basis%tab_basis(1)%nq,1:Psi%Basis%tab_basis(2)%nb)  => psi_g%CVec
                     Qp(inb) = dot_product(psi_gb(:,inb), Psi%Basis%tab_basis(1)%w*Psi%Basis%tab_basis(1)%x*psi_gb(:,inb))
                     Norm(inb) = dot_product(psi_gb(:,inb), Psi%Basis%tab_basis(1)%w*psi_gb(:,inb))
                       if(Norm(inb) /=   ZERO ) then
                      Qp(inb)  = Qp(inb)/Norm(inb)
                      end if
                   end do
                  !Do iq = 1,Psi%Basis%tab_basis(1)%nq
                  !   write(666,*)   psi_gb(iq,1) , psi_gb(iq,2)
                  !End Do

                   print*,Qp ,Norm
                IF (debug) THEN
                !	write(out_unitp,*) 'END Qpop
                flush(out_unitp)
                END IF
         END SUBROUTINE  Qpop





end module
