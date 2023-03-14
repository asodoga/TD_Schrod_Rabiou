
      MODULE param_WP0_m
      USE UtilLib_m
      IMPLICIT NONE

        TYPE GWP1D_t

        real (kind=Rk)                  :: sigma = ONETENTH  ! width of WP0
        real (kind=Rk)                  :: Q0    = ZERO      ! position of WP0
        real (kind=Rk)                  :: imp_k = ZERO      ! impultion for WP0
        real (kind=Rk)                  :: phase = ZERO      ! phase for WP0


        END TYPE GWP1D_t

        TYPE GWP_t
          integer                     :: ndim           = 1
          complex (kind=Rk)           :: Coef           = CZERO
          integer                     ::    Elecindex   = 1

          TYPE (GWP1D_t), allocatable :: tab_GWP1D(:)

        END TYPE GWP_t

        CONTAINS

   RECURSIVE SUBROUTINE Read_GWP1D(GWP1D)
    TYPE (GWP1D_t), intent(inout)          :: GWP1D

    !------ initial WP definition -----------------------------
    !     GWP(Q)=exp[-((Q-Qeq)/sigma)2+i*imp_k*(Q-Qeq)+i*phase]
    real (kind=Rk)                         :: sigma,imp_k,Qeq,phase
    integer                                :: Rerr

    NAMELIST /defWP0/ sigma,Qeq,imp_k,phase


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Read_GWP1D'
    logical, parameter :: debug =.FALSE.
    !logical, parameter :: debug =.TRUE.
    !-----------------------------------------------------------

    sigma    = ONETENTH
    Qeq      = ZERO
    imp_k    = ZERO
    phase    = ZERO

    read(in_unitp,defWP0,iostat=Rerr)
    if (Rerr /= 0) THEN
      write(out_unitp,*) ' ERROR in ',name_sub
      write(out_unitp,*) ' problem while reading the namelist "defWP0"'
      write(out_unitp,defWP0)
      STOP 'ERROR in Read_GWP1D: problem while reading the namelist "defWP0"'
    end if
          GWP1D = GWP1D_t(sigma,Qeq,imp_k,phase)
       

  END SUBROUTINE Read_GWP1D

         SUBROUTINE Write_GWP1D(GWP1D)
           TYPE (GWP1D_t), intent(in) :: GWP1D
        
           !----- for debuging --------------------------------------------------
           character (len=*), parameter :: name_sub='Write_GWP1D'
           logical, parameter :: debug =.FALSE.
           !logical, parameter :: debug =.TRUE.
           !-----------------------------------------------------------
        
           IF (debug) write(out_unitp,*) 'Q0,sigma,imp_k,phase',           &
                  GWP1D%Q0,GWP1D%sigma,GWP1D%imp_k,GWP1D%phase
        
           write(out_unitp,*)                                             &
                'sigma=',GWP1D%sigma,'Qeq=',GWP1D%Q0 ,'imp_k =',GWP1D%imp_k,'phase=',GWP1D%phase
        
         END SUBROUTINE Write_GWP1D

         FUNCTION calc_GWP1D(GWP1D,Q)
           complex (kind=Rk)             :: calc_GWP1D
           real (kind=Rk), intent(in)    :: Q
           TYPE (GWP1D_t),    intent(in) :: GWP1D
        
           real (kind=Rk) :: ze,zk,DQ
           !----- for debuging --------------------------------------------------
           character (len=*), parameter :: name_sub='calc_GWP1D'
           logical, parameter :: debug =.FALSE.
           !logical, parameter :: debug =.TRUE.
           !-----------------------------------------------------------
           DQ = Q-GWP1D%Q0
           ze = (DQ/GWP1D%sigma)**2
           zk = mod(DQ*GWP1D%imp_k + GWP1D%phase ,TWO*pi)
        
           calc_GWP1D = exp(-ze + EYE*zk) / sqrt(sqrt(pi/TWO)*GWP1D%sigma)
        
         END FUNCTION calc_GWP1D
        
        SUBROUTINE Read_GWP(GWP,nio)
          TYPE (GWP_t), intent(inout) :: GWP
          integer,intent(in)          :: nio
          !------ initial WP definition -----------------------------
          integer                      :: i,Rerr
          integer                      ::  ndim,Elecindex

          complex (kind=Rk)            :: Coef

          NAMELIST /defGWP/ ndim,Elecindex,Coef

          !----- for debuging --------------------------------------------------
          character (len=*), parameter :: name_sub='Read_GWP'
          logical, parameter :: debug =.FALSE.
          !logical, parameter :: debug =.TRUE.
          !-----------------------------------------------------------


         ndim            = 1
         Coef            = CONE
         Elecindex       = 1 

          read(nio,nml=defGWP,iostat=Rerr)
          IF (Rerr /= 0) THEN  
            print*,'Rerr',Rerr
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' problem while reading the namelist "defGWP"'
            write(out_unitp,defGWP)
            STOP 'ERROR in Read_GWP: problem while reading the namelist "defGWP"'
          END IF
           GWP%ndim          = ndim
           GWP%Coef          = Coef
           GWP%Elecindex     = Elecindex

          allocate(GWP%tab_GWP1D(ndim))
          DO i=1,ndim
            CALL Read_GWP1D(GWP%tab_GWP1D(i))
          END DO

        END SUBROUTINE Read_GWP

       FUNCTION calc_GWP(GWP,Q)
         complex (kind=Rk)          :: calc_GWP
         real (kind=Rk), intent(in) :: Q(:)
         TYPE (GWP_t),      intent(in) :: GWP

         integer :: i
         !----- for debuging --------------------------------------------------
         character (len=*), parameter :: name_sub='calc_GWP'
         logical, parameter :: debug =.FALSE.
         !logical, parameter :: debug =.TRUE.
         !-----------------------------------------------------------

         calc_GWP = CONE
         DO i=1,size(GWP%tab_GWP1D)
           calc_GWP = calc_GWP * calc_GWP1D(GWP%tab_GWP1D(i),Q(i))
         END DO
         calc_GWP = calc_GWP * GWP%Coef

       END FUNCTION calc_GWP


      SUBROUTINE Write_GWP(GWP)
        TYPE (GWP_t), intent(in) :: GWP
      
        integer :: i
        !----- for debuging --------------------------------------------------
        character (len=*), parameter :: name_sub='Write_GWP'
        logical, parameter :: debug =.FALSE.
        !logical, parameter :: debug =.TRUE.
        !-----------------------------------------------------------
      
        write(out_unitp,*) 'ndim,Elecindex,Coef',GWP%ndim , GWP%Elecindex , GWP%Coef
        write(out_unitp,'(a)') ' Qeq         sigma       imp_k       phase'
        DO i=1,size(GWP%tab_GWP1D)
          CALL Write_GWP1D(GWP%tab_GWP1D(i))
        END DO
      
      END SUBROUTINE Write_GWP

         SUBROUTINE Read_tab_GWP(tab_GWP,nb_GWP,nio)
        
           TYPE (GWP_t), allocatable, intent(inout) :: Tab_GWP(:)
           integer,                   intent(in)    :: nb_GWP
           integer,                   intent(in)    :: nio
        
        
        
           !------ initial WP definition -----------------------------
             integer                :: i
        
           !----- for debuging --------------------------------------------------
           character (len=*), parameter :: name_sub='Read_tab_GWP'
           logical, parameter :: debug =.FALSE.
           !logical, parameter :: debug =.TRUE.
           !-----------------------------------------------------------
        
           allocate(Tab_GWP(nb_GWP))
           DO i=1,nb_GWP
             CALL Read_GWP(Tab_GWP(i),nio)
           END DO        
           CALL Write_Tab_GWP(Tab_GWP)
        
        
         END SUBROUTINE Read_tab_GWP


         SUBROUTINE Write_Tab_GWP(Tab_GWP)
           TYPE (GWP_t), intent(in) :: Tab_GWP(:)
        
           integer :: i
           !----- for debuging --------------------------------------------------
           character (len=*), parameter :: name_sub='Write_Tab_GWP'
           logical, parameter :: debug =.FALSE.
           !logical, parameter :: debug =.TRUE.
           !-----------------------------------------------------------
        
           DO i=1,size(Tab_GWP)
             CALL Write_GWP(Tab_GWP(i))
           END DO
        
         END SUBROUTINE Write_Tab_GWP

         SUBROUTINE calc_Tab_GWP(Tab_GWP,Q,psi)
          complex (kind=Rk) ,intent(inout)   :: psi(:)
          real (kind=Rk), intent(in)         :: Q(:)
          TYPE (GWP_t),      intent(in)      :: Tab_GWP(:)
          integer :: iGWP0
          !----- for debuging --------------------------------------------------
          character (len=*), parameter :: name_sub='calc_GWP'
          logical, parameter :: debug =.FALSE.
          !logical, parameter :: debug =.TRUE.
          !-----------------------------------------------------------
          psi(:) = CZERO     

          DO iGWP0 = 1,size(Tab_GWP),1
               psi(Tab_GWP(iGWP0)%Elecindex)  = calc_GWP(Tab_GWP(iGWP0),Q)
          END DO
          
        
        END SUBROUTINE calc_Tab_GWP


      END MODULE param_WP0_m
