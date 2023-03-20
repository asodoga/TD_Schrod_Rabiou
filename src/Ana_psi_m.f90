 module Ana_psi_m
    USE Basis_m
    USE NDindex_m
    USE psi_m
    implicit none
     private
     public:: Population,Qpop,Calc_std_dev_AVQ_1D

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
               AVQ = AVQ/(Norm*Norm)
               X   = X/(Norm*Norm)
               STD_DQ =sqrt(X-AVQ*AVQ)
             SQ = ONE/(STD_DQ*sqrt(TWO))
               Print *,"<psi|Q|psi> = ",AVQ ,"<psi|Q**2|psi> = ",X,"sqrt(<psi|Q**2|psi> - <psi|Q|psi> )= ", STD_DQ,'SQ=',SQ
             CALL dealloc_psi(psi)
             IF (debug) THEN
                ! write(out_unitp,*) 'END AVQ,STD_DQ'
                 flush(out_unitp)
             END IF

         END SUBROUTINE Calc_std_dev_AVQ_1D


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
