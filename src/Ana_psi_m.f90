 module Ana_psi_m
    USE Basis_m
    USE NDindex_m
    USE psi_m
    implicit none
     private
     public:: Population,Qpop
     public:: Calc_AVQ_nD0,Calc_AVQ_nD 

     contains



      
     SUBROUTINE Calc_AVQ_1D(psi_in,AVQ1,AVQel1,SQ1,SQel1,ib)
      USE UtilLib_m
      logical,         parameter               :: debug = .false.
      TYPE(Psi_t), intent(in)                  :: psi_in
      TYPE(Psi_t)                              :: psi
      complex (kind=Rk),allocatable            :: psi_gb(:,:)
      logical                                  :: Endloop_q
      real(kind=Rk),intent(inout),optional     :: AVQ1,SQ1
      real(kind=Rk),intent(inout),optional     :: AVQel1(:),SQel1(:)
      
      real(kind=Rk),allocatable                :: AVQel(:),SQel(:)
      real(kind=Rk)                            :: AVQ,SQ
      integer  ,intent(in)                     :: ib
      real(kind=Rk)                            :: WnD,X
      real(kind=Rk),allocatable                :: N(:)
      integer,        allocatable              :: Tab_iq(:)
      integer                                  :: iq,inbe,inb
      IF (debug) THEN
          write(out_unitp,*) 'Beging AVQ'
          flush(out_unitp)
      END IF 
      allocate(N(psi_in%Basis%tab_basis(size(psi_in%Basis%tab_basis))%nb))
      allocate(AVQel(psi_in%Basis%tab_basis(size(psi_in%Basis%tab_basis))%nb))
      allocate(SQel(psi_in%Basis%tab_basis(size(psi_in%Basis%tab_basis))%nb))
      CALL init_psi(psi,   psi_in%Basis,    cplx=.TRUE.   ,grid =.true.)

      IF(psi_in%Grid) then
          psi%CVec(:)= psi_in%CVec(:)
      ELSE
          CALL BasisTOGrid_nD_cplx(psi%CVec,psi_in%CVec,psi_in%Basis)
      END IF

      Allocate(Psi_gb(psi%Basis%nq,psi%Basis%tab_basis(size(psi%Basis%tab_basis))%nb))
      Allocate(Tab_iq(size(Psi%Basis%tab_basis)-1))
      psi_gb(:,:) = reshape(psi%CVec,shape= [psi%Basis%nq,psi%Basis%tab_basis(size(psi%Basis%tab_basis))%nb])

      X = ZERO
      N(:) = ZERO

      DO inbe = 1,psi%Basis%tab_basis(size(psi%Basis%tab_basis))%nb !electronic state
        
          AVQel(inbe) = ZERO
          Call Init_tab_ind(Tab_iq,psi%Basis%NDindexq)
          Iq = 0
          DO
              Iq = Iq+1
              CALL increase_NDindex(Tab_iq,psi%Basis%NDindexq,Endloop_q)
              IF (Endloop_q) exit
              WnD= ONE
              DO inb = 1,size(psi%Basis%tab_basis)-1
                  WnD  =WnD*psi%Basis%tab_basis(inb)%w(tab_iq(inb))
              END DO
              X = psi%Basis%tab_basis(ib)%x(tab_iq(ib))
              N(inbe) = N(inbe) + conjg(psi_gb(iq,inbe))*psi_gb(iq,inbe)*WnD
              AVQel(inbe) = AVQel(inbe) + conjg(psi_gb(iq,inbe))*(X*WnD)*psi_gb(iq,inbe)
              SQel(inbe)  = SQel(inbe)  + conjg(psi_gb(iq,inbe))*(X*X*WnD)*psi_gb(iq,inbe)
          END DO
      END DO
      AVQ = sum(AVQel)/(Sum(N)**2)
      SQ  = sum(SQel)/(Sum(N)**2)
      SQ  = sqrt(SQ-AVQ*AVQ)
      SQ  = ONE/(SQ*sqrt(TWO))  
      DO inbe = 1,Psi%Basis%tab_basis(size(Psi%Basis%tab_basis))%nb !electronic state
          if(N(inbe) /= ZERO)  AVQel(inbe) = AVQel(inbe)/(N(inbe)**2)
      END DO
      if(present(AVQ1))    AVQ1     = AVQ
      if(present(AVQel1))  AVQel1   = AVQel
      if(present(SQ1))     SQ1      = SQ
      if(present(AVQel1))  SQel1    = SQel


     ! if(present(AVQ1))  print*,   'AVQ1=',    AVQ1
     ! if(present(AVQel1)) print*,  'AVQel1=',  AVQel1 
     ! if(present(SQ1))    print*,  'SQ1=',     SQ1     
      
      Deallocate(Tab_iq)
      Deallocate(Psi_gb)
      CALL dealloc_psi(psi)
      IF (debug) THEN
          write(out_unitp,*) 'END AVQ'
          flush(out_unitp)
      END IF
    END SUBROUTINE 


    SUBROUTINE Calc_AVQ_nD0(psi0,AVQ,AVQel,SQ,SQel)
      USE UtilLib_m
      
      type(psi_t) , intent(in) ,target                        :: psi0
      type(psi_t)              ,target                        :: psi
      real(kind=Rk),intent(inout) ,optional                   :: AVQ(:),SQ(:)
      real(kind=Rk),intent(inout)  ,optional                  :: AVQel(:,:),SQel(:,:)
      logical,         parameter                              :: debug = .true.
      integer                                                 :: Inb,Ndim
      
      
      IF (debug) THEN
          flush(out_unitp)
      END IF
      
      Ndim = size(psi0%Basis%tab_basis)     
      call  init_psi(psi,   psi0%Basis,    cplx=.TRUE.   ,grid =.true.)
       
      IF(psi0%Grid) then
        psi%CVec(:)= psi0%CVec(:)
      ELSE
        call BasisTOGrid_nD_cplx(psi%CVec,psi0%CVec,psi0%Basis)
      END IF
     
      DO  Inb = 1,Ndim-1
         
         if ( present(AVQel) ) call Calc_AVQ_1D(psi_in=psi,AVQel1=AVQel(Inb,:),ib=Inb)        
         if ( present(AVQ) )   call Calc_AVQ_1D(psi_in=psi,AVQ1=AVQ(Inb),ib=Inb)

         if ( present(SQel) ) call Calc_AVQ_1D(psi_in=psi,SQel1=SQel(Inb,:),ib=Inb)  
         if ( present(SQ) )  call Calc_AVQ_1D(psi_in=psi,SQ1=SQ(Inb),ib=Inb)
                 
          
      END DO
         
        write(out_unitp,*) '<psi/Q/psi> =',AVQ
        write(out_unitp,*) 'SQ =',SQ
      IF (debug) THEN
          flush(out_unitp)
      END IF
      CALL dealloc_psi(psi)
  END SUBROUTINE



  SUBROUTINE Calc_AVQ_nD(psi0,AVQ,SQ)
    USE UtilLib_m
    
    type(psi_t) , intent(in) ,target                        :: psi0
    type(psi_t)              ,target                        :: psi
    real(kind=Rk),intent(inout)                             :: AVQ(:),SQ(:)
    real(kind=Rk),allocatable                               :: Q(:,:),W(:),X(:)
    real(kind=Rk)                                           :: Norm
    logical,         parameter                              :: debug = .true.
    integer                                                 :: Inb,Ndim,Iq
    
    
    IF (debug) THEN
        flush(out_unitp)
    END IF
    
    Ndim = size(psi0%Basis%tab_basis)-1 
    call Calc_Norm_OF_Psi(psi0,Norm)
    allocate(X(Ndim))   
    call  init_psi(psi,   psi0%Basis,    cplx=.TRUE.   ,grid =.true.)
    call Calc_Q_grid(Q=Q,Basis=psi0%Basis,WnD=W)
     
    IF(psi0%Grid) then
      psi%CVec(:)= psi0%CVec(:)
    ELSE
      call BasisTOGrid_nD_cplx(psi%CVec,psi0%CVec,psi0%Basis)
    END IF
    SQ  = ZERO
    X   = ZERO
    AVQ = ZERO
    DO  Inb = 1,Ndim
       
       AVQ(Inb) = dot_product(psi%CVec(:),W(:)*Q(:,Inb)*psi%CVec(:)) 
       X(Inb)   = dot_product(psi%CVec(:),W(:)*Q(:,Inb)*Q(:,Inb)*psi%CVec(:)) 
       X(Inb)   = X(Inb)/(Norm*Norm)
       AVQ(Inb) = AVQ(Inb)/(Norm*Norm)
       SQ(Inb)  = sqrt(X(Inb)-AVQ(Inb)*AVQ(Inb))
       SQ(Inb)  = ONE/(SQ(Inb)*sqrt(TWO)) 
        
    END DO
       
      write(out_unitp,*) '<psi/Q/psi> =',AVQ
      write(out_unitp,*) 'SQ =',SQ
    IF (debug) THEN
        flush(out_unitp)
    END IF
    CALL dealloc_psi(psi)
END SUBROUTINE




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
