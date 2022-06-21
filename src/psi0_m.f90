module psi0_m
  USE NumParameters_m
  USE Basis_m, only : Basis_t
  USE psi_m
  implicit none

   TYPE  ::  gwp_dim_para_t
     integer                         ::nb_GWP
     integer                         :: ndim
     integer                         ::  I_ElecChannel

   END TYPE   gwp_dim_para_t


  TYPE  :: gwp_para_t
     real(KIND=Rk)      ::Q0
     real(KIND=Rk)      ::K
     real(KIND=Rk)      ::phase
     real(KIND=Rk)      ::DQ
   END TYPE    gwp_para_t


  TYPE  ::  para_psi0_t

     complex(KIND=Rk),allocatable          :: vector_Coef(:)
     type(gwp_para_t),allocatable            :: nd_wavepaket_para(:)
     type(gwp_dim_para_t),allocatable        :: nd_wavepaket_dim_para(:)

  END TYPE   para_psi0_t


   public ::  Write_psi0, Read_psi0
contains


 SUBROUTINE Read_psi0(paragwp0,nio)
USE UtilLib_m
  implicit none
   integer,             intent(in)     :: nio
  logical,             parameter                  :: debug = .true.
  type( para_psi0_t)  ,intent(inout)                :: paragwp0
  real(kind= Rk)                                  :: K,DQ,phase,Q0
  integer                                         :: I_ElecChannel,i,j,ndim,nb_GWP,err_io
  complex(kind= Rk)                               :: scalar_coef

  namelist/paragwp/  nb_GWP  , ndim
  namelist /defGWP / I_ElecChannel, scalar_coef
  namelist /defWP0/  DQ, Q0,  k, phase
  nb_GWP = 1;  ndim  = 1;
  I_ElecChannel = 2  ; scalar_coef   = (1,0);
  DQ=0.2 ; Q0=0.0 ; k=0.0 ; phase=0.


    OPEN(nio, FILE='my_nml.dat',FORM='formatted', ACTION='read')
    read(nio,nml=paragwp,IOSTAT=err_io)
    IF (err_io /= 0) THEN
      write(out_unitp,*) ' ERROR in Reed_psi0'
      print*, err_io
      stop 'problem with nameist'

    else
        print*,'paragwp read normaly'
        write(1,nml=paragwp,IOSTAT=err_io)
        allocate(paragwp0%vector_Coef(nb_GWP))
        allocate(paragwp0%nd_wavepaket_para(ndim))
        allocate(paragwp0%nd_wavepaket_dim_para(nb_GWP))
    END IF

 do i = 1,nb_GWP,1
   read(nio,nml=defGWP,IOSTAT=err_io)
   paragwp0%nd_wavepaket_dim_para(i)%nb_gwp        = nb_GWP
   paragwp0%nd_wavepaket_dim_para(i)%ndim          = ndim
   paragwp0%nd_wavepaket_dim_para(i)%I_ElecChannel = I_ElecChannel
   do j = 1,ndim
         read(nio,nml=defWP0,IOSTAT=err_io)
        paragwp0%nd_wavepaket_para(j)%Q0    = Q0
        paragwp0%nd_wavepaket_para(j)%K     = K
        paragwp0%nd_wavepaket_para(j)%phase = phase
        paragwp0%nd_wavepaket_para(j)%DQ    = DQ
        paragwp0%vector_Coef(i)             = scalar_coef
   end do
 end do


END SUBROUTINE Read_psi0

SUBROUTINE Write_psi0(paragwp0,nio)
 USE UtilLib_m

    type( para_psi0_t)  ,intent(inout)                :: paragwp0
    integer,             intent(in)                   :: nio
    integer                                           ::  i,j
    integer                                           :: ndim,nb_GWP,err_io


    namelist/paragwp/  nb_GWP  , ndim
       nb_gwp = 1
       ndim   = 1
    write(out_unitp,*) '**********************************************************'
     write(out_unitp,*) '***************Beging paragwp0***************************'
    write(out_unitp,*) '**********************************************************'
  read(nio,nml=paragwp,IOSTAT=err_io)
 do i = 1,nb_GWP,1
    write(out_unitp,*) paragwp0%nd_wavepaket_dim_para(i)%nb_gwp
    write(out_unitp,*) paragwp0%nd_wavepaket_dim_para(i)%ndim
    write(out_unitp,*)  paragwp0%nd_wavepaket_dim_para(i)%I_ElecChannel
    do j = 1,ndim
       write(out_unitp,*) paragwp0%nd_wavepaket_para(j)%Q0
       write(out_unitp,*)  paragwp0%nd_wavepaket_para(j)%K
       write(out_unitp,*)  paragwp0%nd_wavepaket_para(j)%phase
       write(out_unitp,*)  paragwp0%nd_wavepaket_para(j)%DQ
       write(out_unitp,*)   paragwp0%vector_Coef(i)
   end do
 end do
    write(out_unitp,*) '**********************************************************'
      write(out_unitp,*) '**************End paragwp0*****************************'
    write(out_unitp,*) '**********************************************************'


END SUBROUTINE Write_psi0
    SUBROUTINE GWP0(psi0,paragwp0,Basis,nio)
      USE NumParameters_m
      USE Basis_m
      TYPE(para_psi0_t),  intent(inout)                   :: paragwp0
      TYPE(psi_t),INTENT(inout)                           :: psi0
      TYPE(Basis_t)  ,INTENT(IN)                          :: Basis
      integer,             intent(in)                     :: nio
      COMPLEX(KIND=Rk), ALLOCATABLE                       :: g1(:,:),g2(:,:)
      INTEGER                                             :: IQ , IB
      REAL(KIND= Rk)                                      ::dot_prdct,Norm
      integer                                             ::  i,j,k
      integer                                             :: ndim,nb_GWP

       namelist/paragwp/  nb_GWP  , ndim
       nb_gwp = 1
       ndim   = 1
       CALL Read_psi0( paragwp0,nio)
       CALL write_psi0( paragwp0,nio)
      stop 'ici'
      read(nio,nml=paragwp)
        ALLOCATE(g1(Basis%tab_basis(1)%nq, Basis%tab_basis(2)%nb))

        ALLOCATE(g2(Basis%tab_basis(1)%nq, Basis%tab_basis(2)%nb))

        do i = 1,nb_GWP,1
            do k = 1,paragwp0%nd_wavepaket_dim_para(i)%I_ElecChannel,1
               do j = 1,ndim,1
                 g1(:,k)= g1(:,k)*exp(-((basis%tab_basis(1)%x(:)-paragwp0%nd_wavepaket_para(j)%Q0)/&
                        & paragwp0%nd_wavepaket_para(j)%DQ)**2  &
                        &+ EYE*paragwp0%nd_wavepaket_para(j)%K*(basis%tab_basis(1)%x(:)-paragwp0%nd_wavepaket_para(j)%Q0) &
                        &+EYE*paragwp0%nd_wavepaket_para(j)%phase )
               end do
            end do
            g2(:,k) = g2(:,k)+ g1(:,k)*paragwp0%vector_Coef(i)
        end do
            psi0%CVec(:)  = reshape( g2,[Basis%tab_basis(1)%nq*Basis%tab_basis(2)%nb])
            CALL Calc_dot_product(psi0%CVec,dot_prdct,Basis,.true.,.false.)
            psi0%CVec(:)  =  psi0%CVec(:)/SQRT(dot_prdct)
            CALL Calc_dot_product(psi0%CVec,dot_prdct,Basis,.true.,.true.)



      END SUBROUTINE GWP0




end module psi0_m
