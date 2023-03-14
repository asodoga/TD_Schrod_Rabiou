
      MODULE param_WP0_m
      USE UtilLib_m
      USE psi_m
      IMPLICIT NONE

        TYPE GWP1D_t
          ! character(len=:),allocatable    :: name
          ! integer                         :: nb_GWP
          ! complex (kind=Rk)               :: Coef  = CONE

        real (kind=Rk)                  :: sigma = ONETENTH  ! width of WP0
        real (kind=Rk)                  :: Q0    = ZERO      ! position of WP0
        real (kind=Rk)                  :: imp_k = ZERO      ! impultion for WP0
        real (kind=Rk)                  :: phase = ZERO      ! phase for WP0
        TYPE( GWP1D_t) ,allocatable     :: tab(:)     !  for more than one nbGWP


        END TYPE GWP1D_t

        TYPE GWP_t
            integer                   :: ndim           = 1
          complex (kind=Rk)           :: Coef           = CZERO

          TYPE (GWP1D_t), allocatable :: tab_GWP1D(:)

        END TYPE GWP_t

        TYPE param_WP0


        integer             :: nb_WP0             = 1       ! default: 1
        !for each variable Qi : exp[-((Q-Qeq)/sigma)2+i*imp_k*(Q-Qeq)]
        real (kind=Rk), allocatable :: WP0sigma(:)       ! WP0sigma(nb_act1) : sigma for WP0
        real (kind=Rk), allocatable :: WP0Qeq(:)         ! WP0Qeq(nb_act1)   : position of WP0
        real (kind=Rk), allocatable :: WP0imp_k(:)       ! WP0imp_k(nb_act1) : impultion for WP0
        real (kind=Rk), allocatable :: WP0phase(:)       ! WP0imp_k(nb_act1) : phase for WP0
        TYPE (GWP_t), allocatable :: tab_GWP0(:)

        END TYPE param_WP0

        CONTAINS

   RECURSIVE SUBROUTINE Read_GWP1D(GWP1D)
    TYPE (GWP1D_t), intent(inout)          :: GWP1D
    !integer,             intent(in)        :: nio

    !------ initial WP definition -----------------------------
    !     GWP(Q)=exp[-((Q-Qeq)/sigma)2+i*imp_k*(Q-Qeq)+i*phase]
    real (kind=Rk)                         :: sigma,imp_k,Qeq,phase
     ! complex (kind=Rk)                      :: Coef
     ! integer                                :: nb_GWP
     ! character(len=:),allocatable           :: name
    integer                                :: Rerr!,inb

    NAMELIST /defWP0/ sigma,Qeq,imp_k,phase


    !----- for debuging --------------------------------------------------
    character (len=*), parameter :: name_sub='Read_GWP1D'
    logical, parameter :: debug =.FALSE.
    !logical, parameter :: debug =.TRUE.
    !-----------------------------------------------------------
    !name     ='0'
    !nb_GWP   = 1
    !Coef     = CONE

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
       !if(nb_GWP>=2)then
       !     GWP1D%name = 'defWP'
       !    allocate(GWP1D%tab(nb_GWP))
       !    do inb = 1,nb_GWP,1
       !        call  Read_GWP1D(GWP1D%tab(inb))
       !    end do
       !else
          !GWP1D%name = 'defWP0'
          GWP1D = GWP1D_t(sigma,Qeq,imp_k,phase)
       !end if

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
        
        SUBROUTINE Read_GWP(GWP)
          TYPE (GWP_t), intent(inout) :: GWP

          !------ initial WP definition -----------------------------
          integer                      :: i,Rerr
          integer                      ::  ndim

          complex (kind=Rk)            :: Coef

          NAMELIST /defGWP/ ndim,Coef

          !----- for debuging --------------------------------------------------
          character (len=*), parameter :: name_sub='Read_GWP'
          logical, parameter :: debug =.FALSE.
          !logical, parameter :: debug =.TRUE.
          !-----------------------------------------------------------


         ndim            = 1
         Coef            = CONE

          read(in_unitp,defGWP,iostat=Rerr)
          stop 'cc'
          IF (Rerr /= 0) THEN
            write(out_unitp,*) ' ERROR in ',name_sub
            write(out_unitp,*) ' problem while reading the namelist "defGWP"'
            write(out_unitp,defGWP)
            STOP 'ERROR in Read_GWP: problem while reading the namelist "defGWP"'
          END IF

           GWP%ndim          = ndim
           GWP%Coef          = Coef

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
      
        write(out_unitp,*) 'ndim,Coef',GWP%ndim,GWP%Coef
        write(out_unitp,'(a)') ' Qeq         sigma       imp_k       phase'
        DO i=1,size(GWP%tab_GWP1D)
          CALL Write_GWP1D(GWP%tab_GWP1D(i))
        END DO
      
      END SUBROUTINE Write_GWP
        !
        ! SUBROUTINE Read_tab_GWP(tab_GWP,nb_GWP,mole,nb_bi,nb_elec)
        !   USE mod_Coord_KEO
        !
        !   TYPE (GWP_t), allocatable, intent(inout) :: tab_GWP(:)
        !   integer,                   intent(in)    :: nb_GWP,nb_elec,nb_bi
        !
        !   TYPE (CoordType),          intent(in)    :: mole
        !
        !
        !   !------ initial WP definition -----------------------------
        !     integer                :: iGWP,iQdyn,i
        !
        !   !----- for debuging --------------------------------------------------
        !   character (len=*), parameter :: name_sub='Read_tab_GWP'
        !   logical, parameter :: debug =.FALSE.
        !   !logical, parameter :: debug =.TRUE.
        !   !-----------------------------------------------------------
        !
        !   allocate(tab_GWP(nb_GWP))
        !   DO i=1,nb_GWP
        !     CALL Read_GWP(tab_GWP(i),mole%nb_act1,nb_bi,nb_elec)
        !   END DO
        !
        !   ! check if iQact is defined
        !   DO iGWP=1,size(tab_GWP)
        !     IF (count(tab_GWP(iGWP)%tab_GWP1D(:)%iQact == -1) /= 0) THEN ! iQact is not defined => transfert iQdyn to iQact
        !       DO i=1,mole%nb_act1
        !         iQdyn = tab_GWP(iGWP)%tab_GWP1D(i)%iQdyn
        !         tab_GWP(iGWP)%tab_GWP1D(i)%iQact = mole%liste_QdynTOQact(iQdyn)
        !       END DO
        !     END IF
        !   END DO
        !
        !   CALL Write_Tab_GWP(tab_GWP)
        !
        !
        ! END SUBROUTINE Read_tab_GWP
        ! SUBROUTINE Read_tab_GWP_v0(tab_GWP,nb_GWP,ndim,nb_bi,nb_elec)
        !   TYPE (GWP_t), allocatable, intent(inout) :: tab_GWP(:)
        !   integer,                   intent(in)    :: nb_GWP,ndim,nb_elec,nb_bi
        !
        !
        !   !------ initial WP definition -----------------------------
        !   integer                     :: i
        !
        !   !----- for debuging --------------------------------------------------
        !   character (len=*), parameter :: name_sub='Read_tab_GWP_v0'
        !   logical, parameter :: debug =.FALSE.
        !   !logical, parameter :: debug =.TRUE.
        !   !-----------------------------------------------------------
        !
        !   allocate(tab_GWP(nb_GWP))
        !   DO i=1,nb_GWP
        !     CALL Read_GWP(tab_GWP(i),ndim,nb_bi,nb_elec)
        !   END DO
        !
        ! END SUBROUTINE Read_tab_GWP_v0
        ! SUBROUTINE Write_Tab_GWP(tab_GWP)
        !   TYPE (GWP_t), intent(in) :: tab_GWP(:)
        !
        !   integer :: i
        !   !----- for debuging --------------------------------------------------
        !   character (len=*), parameter :: name_sub='Write_Tab_GWP'
        !   logical, parameter :: debug =.FALSE.
        !   !logical, parameter :: debug =.TRUE.
        !   !-----------------------------------------------------------
        !
        !   DO i=1,size(tab_GWP)
        !     CALL Write_GWP(tab_GWP(i))
        !   END DO
        !
        ! END SUBROUTINE Write_Tab_GWP
        !================================================================
        !
        !    alloc / dealloc param_WP0
        !
        !================================================================
        !     !!@description: alloc param_WP0
        !     !!@param: para_poly
        !     SUBROUTINE alloc_param_WP0(para_WP0,WP0Grid_Gaussian,WP0_CleanChannel)
        !     IMPLICIT NONE
        !     TYPE (param_WP0), intent(inout) :: para_WP0
        !     logical, intent(in) :: WP0Grid_Gaussian,WP0_CleanChannel
        !
        !       IF (para_WP0%nb_act1 > 0 .AND. WP0Grid_Gaussian) THEN
        !         IF (.NOT. allocated(para_WP0%WP0sigma)) THEN
        !           CALL alloc_NParray(para_WP0%WP0sigma,[para_WP0%nb_act1],    &
        !                             "para_WP0%WP0sigma","alloc_param_WP0")
        !         END IF
        !
        !         IF (.NOT. allocated(para_WP0%WP0Qeq)) THEN
        !           CALL alloc_NParray(para_WP0%WP0Qeq,[para_WP0%nb_act1],      &
        !                             "para_WP0%WP0Qeq","alloc_param_WP0")
        !         END IF
        !
        !         IF (.NOT. allocated(para_WP0%WP0imp_k)) THEN
        !           CALL alloc_NParray(para_WP0%WP0imp_k,[para_WP0%nb_act1],    &
        !                             "para_WP0%WP0imp_k","alloc_param_WP0")
        !         END IF
        !
        !         IF (.NOT. allocated(para_WP0%WP0phase)) THEN
        !           CALL alloc_NParray(para_WP0%WP0phase,[para_WP0%nb_act1],    &
        !                             "para_WP0%WP0phase","alloc_param_WP0")
        !         END IF
        !
        !       END IF
        !
        !       IF (WP0_CleanChannel .AND.                                      &
        !           .NOT. allocated(para_WP0%WP0_CleanChannellist)) THEN
        !         CALL alloc_NParray(para_WP0%WP0_CleanChannellist,               &
        !                                    [para_WP0%WP0_nb_CleanChannel],  &
        !                         "para_WP0%WP0_CleanChannellist","alloc_param_WP0")
        !       END IF
        !
        !     END SUBROUTINE alloc_param_WP0
        !
        !     !!@description: dealloc param_WP0
        !     !!@param: para_poly
        !     SUBROUTINE dealloc_param_WP0(para_WP0)
        !     IMPLICIT NONE
        !     TYPE (param_WP0), intent(inout) :: para_WP0
        !
        !       IF (allocated(para_WP0%WP0sigma)) THEN
        !         CALL dealloc_NParray(para_WP0%WP0sigma,                         &
        !                             "para_WP0%WP0sigma","dealloc_param_WP0")
        !       END IF
        !
        !       IF (allocated(para_WP0%WP0Qeq)) THEN
        !         CALL dealloc_NParray(para_WP0%WP0Qeq,                           &
        !                             "para_WP0%WP0Qeq","dealloc_param_WP0")
        !       END IF
        !
        !       IF (allocated(para_WP0%WP0imp_k)) THEN
        !         CALL dealloc_NParray(para_WP0%WP0imp_k,                         &
        !                             "para_WP0%WP0imp_k","dealloc_param_WP0")
        !       END IF
        !
        !       IF (allocated(para_WP0%WP0phase)) THEN
        !         CALL dealloc_NParray(para_WP0%WP0phase,                         &
        !                             "para_WP0%WP0phase","dealloc_param_WP0")
        !       END IF
        !
        !       IF (allocated(para_WP0%WP0_CleanChannellist)) THEN
        !         CALL dealloc_NParray(para_WP0%WP0_CleanChannellist,             &
        !                           "para_WP0%WP0_CleanChannellist","dealloc_param_WP0")
        !       END IF
        !     END SUBROUTINE dealloc_param_WP0
        !
      END MODULE param_WP0_m
