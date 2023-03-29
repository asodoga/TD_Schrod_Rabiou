 module Ana_psi_m
    USE Basis_m
    USE NDindex_m
    USE psi_m
    implicit none
     private
     public:: Population,Qpop,Calc_std_dev_AVQ_1D,Calc_std_dev_AVQ_nD
     public:: Calc_AVQ_1D,Calc_AVQ_nD

     contains


     SUBROUTINE Calc_AVQ_1D(AVQ,SQ,GGB,Q,W,Norm)
      USE UtilLib_m
      
       complex (kind=Rk), intent(in)            :: GGB(:,:,:)
       Real (kind=Rk), intent(in)               :: Q(:),W(:),Norm
        Real (kind=Rk), intent(inout)           :: AVQ,SQ
        Real (kind=Rk)                          :: SDR_D
       logical          , parameter             :: debug = .true.
       integer                                  :: i1,i3,iq,ib
   
         IF (debug) THEN
           flush(out_unitp)
         END IF
   
         DO i3 = 1,ubound(GGB, dim=3)
         DO i1 = 1,ubound(GGB, dim=1)
   
           AVQ =  dot_product( GGB(i1,:,i3) ,Q(:)*W(:)*GGB(i1,:,i3))
           SDR_D =  dot_product( GGB(i1,:,i3) , Q(:)*Q(:)*W(:)*GGB(i1,:,i3))
   
         END DO
         END DO
   
         AVQ   = AVQ/(Norm*Norm)
         SDR_D   =  SDR_D/(Norm*Norm)
         SDR_D =sqrt( SDR_D-AVQ*AVQ)
         SQ = ONE/(SDR_D*sqrt(TWO))
   
         IF (debug) THEN
           flush(out_unitp)
         END IF
     END SUBROUTINE  Calc_AVQ_1D



     SUBROUTINE Calc_AVQ_nD(psi0,AVQ,SQ)
      USE UtilLib_m
      
      type(psi_t) , intent(in) ,target                        :: psi0
      type(psi_t)              ,target                        :: psi
      real(kind=Rk),intent(inout)                             :: AVQ(:),SQ(:)
      complex (kind=Rk), pointer                              :: GGB(:,:,:)
      real(kind=Rk),pointer                                   :: Q(:),W(:)
      logical,         parameter                              :: debug = .true.
      integer                                                 :: Inb,Ndim
      real(kind=Rk)                                           :: Norm
      integer , allocatable                                   :: Ib1(:),Ib2(:),Iq3(:),Iq1(:),Iq2(:),Ib3(:)
      IF (debug) THEN
          flush(out_unitp)
      END IF
      Ndim = size(psi0%Basis%tab_basis)
      call Calc_iqib( Ib1,Ib2,Ib3,Iq1,Iq2,Iq3,Ndim,psi0%Basis) 
      call Calc_Norm_OF_Psi(psi0,Norm)
      
       CALL init_psi(psi,   psi0%Basis,    cplx=.TRUE.   ,grid =.true.)
               IF(psi0%Grid) then
                   psi%CVec(:)= psi0%CVec(:)
               ELSE
                   CALL BasisTOGrid_nD_cplx(psi%CVec,psi0%CVec,psi0%Basis)
               END IF
     
         DO   Inb = 1,Ndim-1
         
         
               Q(1:Iq2(inb)) => psi%Basis%tab_basis(inb)%X
               W(1:Iq2(inb)) => psi%Basis%tab_basis(inb)%W
               GGB( 1:Iq1(inb),1:Iq2(inb),1:Iq3(inb)) => psi%CVec
               
               CALL Calc_AVQ_1D(AVQ(Inb),SQ(Inb),GGB,Q,W,Norm)
              
                      
          END DO
          write(out_unitp,*) '<psi/Q/psi> =',AVQ
          write(out_unitp,*)  'SQ =',SQ
      
      IF (debug) THEN
          flush(out_unitp)
      END IF
      deallocate (Iq1,Iq2,Iq3,Ib1,Ib2,Ib3)
  END SUBROUTINE Calc_AVQ_nD




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
               AVQ = AVQ/(Norm*Norm)
               X   = X/(Norm*Norm)
               STD_DQ =sqrt(X-AVQ*AVQ)
             SQ = ONE/(STD_DQ*sqrt(TWO))
              ! Print *,"<psi|Q|psi> = ",AVQ ,"<psi|Q**2|psi> = ",X,"sqrt(<psi|Q**2|psi> - <psi|Q|psi> )= ", STD_DQ,'SQ=',SQ
             CALL dealloc_psi(psi)
             IF (debug) THEN
                ! write(out_unitp,*) 'END AVQ,STD_DQ'
                 flush(out_unitp)
             END IF

         END SUBROUTINE Calc_std_dev_AVQ_1D





         SUBROUTINE Calc_std_dev_AVQ_nD(psi,AVQ,SQ)
          USE UtilLib_m
          logical,         parameter      :: debug = .false.
          TYPE(Psi_t), intent(in)         :: psi
          real(kind=Rk),intent(inout)     :: AVQ(:),SQ(:)
          integer                         :: Ib,Ndim
          IF (debug) THEN
              !flush(out_unitp)
          END IF
          Ndim =size(psi%Basis%tab_basis)
          Do Ib = 1,Ndim
            If(psi%Basis%tab_basis(Ib)%Basis_name == 'el') then
              AVQ(Ib)=ZERO
              SQ(Ib)= ONE
            else
            call Calc_std_dev_AVQ_1D(psi,Ib,AVQ(Ib),SQ(Ib))
            End if
          End Do

          Print *,'Ib=',Ib
          Print *,"<psi|Q|psi> = ",AVQ 
          Print *,'SQ=',SQ

          IF (debug) THEN
             
              flush(out_unitp)
          END IF

      END SUBROUTINE Calc_std_dev_AVQ_nD








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
