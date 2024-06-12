module Propa_m
   USE QDUtil_m
   USE psi_m
   USE Op_m
   USE Ana_psi_m
   Use lanczos_m
   USE Auto_corr_m
   Use Hagedorn_m
   USE sub_propa_m
   USE Basis_m

   implicit none

contains

   SUBROUTINE propagation(psif, psi0, propa)
      USE psi_m
      USE Basis_m

      TYPE(psi_t), intent(inout)          :: psif
      TYPE(psi_t), intent(in)             :: psi0
      TYPE(propa_t), intent(inout)        :: propa
      logical, parameter                  :: debug = .true.

      ! variables locales------------------------------------------------------------------
      TYPE(Basis_t)                       ::  Basis
      TYPE(Op_t)                          :: H
      REAL(kind=Rkind)                    :: t, t_deltat, Norm, E,e0
      REAL(kind=Rkind), allocatable       :: Qt(:), SQt(:), pop(:),Pt(:)
      integer, allocatable                :: Tab_iq(:, :)
      complex(kind=Rkind) ,allocatable    :: At(:)
      REAL(kind=Rkind)                    :: aut_func_arg
      complex(kind=Rkind)                 :: aut_func
      TYPE(REDUCED_DENSIRY_t)             :: Rd
      integer                             :: Ndim
      INTEGER                             :: i, nt, nf,nsurf
      TYPE(psi_t)                         :: psi, psi_dt,psi_t0
      if (debug) then

         write (out_unit, *) 'BEGINNIG propagation', propa%t0, propa%tf, propa%delta_t
         ! write(out_unit,*) ''

         write (out_unit, *) '-------------propagation parameters---------------'
         Call write_propa(propa)
      else
         STOP ' check your data!'
         flush (out_unit)

      end if

      call creat_file_unit(nio=10, name='psi', propa=propa)
      call creat_file_unit(nio=11, name='Qt', propa=propa)
      call creat_file_unit(nio=12, name='E', propa=propa)
      call creat_file_unit(nio=13, name='SQt', propa=propa)
      call creat_file_unit(nio=14, name='Norm', propa=propa)
      call creat_file_unit(nio=18, name='pop', propa=propa)
      call creat_file_unit(nio=19, name='Imp_k', propa=propa)
      call creat_file_unit(nio=20, name='alpha', propa=propa)
      call creat_file_unit(nio=21, name='Rd', propa=propa)
      !call creat_file_unit(nio=22, name='maxcoeff', propa=propa)
      !call creat_file_unit(nio=23, name='psi_int', propa=propa)
      call creat_file_unit(nio=24, name='Norm_13', propa=propa)
      call creat_file_unit(nio=25, name='E_13', propa=propa)
      call creat_file_unit(nio=28, name='file_norm_pics', propa=propa)
       call creat_file_unit(nio=26, name='auto_cor', propa=propa)
       call creat_file_unit(nio=27, name='psi_dt_on_basis0', propa=propa)
      
      call init_Basis1_TO_Basis2(Basis, psi0%Basis)
      call construct_primitive_basis(Basis)
      Ndim = size(psi0%Basis%tab_basis) - 1
      nsurf = psi0%Basis%tab_basis(Ndim + 1)%nb
      allocate (Qt(Ndim), SQt(Ndim),Pt(Ndim),At(Ndim))
      allocate (pop(nsurf))
      nt = int((propa%tf - propa%t0)/propa%delta_t)

      call init_psi(psi, psi0%Basis, cplx=.true., grid=.false.)
      call init_psi(psi_dt,psi0%Basis, cplx=.true., grid=.false.)
      call init_psi(psi_t0,Basis, cplx=.true., grid=.false.)

      psi%CVec(:) = psi0%CVec(:)
      psi_t0%CVec(:) = psi0%CVec(:)
      call Calc_tab_Iq0(Tab_Iq,psi0%Basis)
      call Set_Op(H, psi0%Basis,Tab_Iq)

      ! call Calc_reduced_density(Rd,psi%CVec,psi%Basis)
      ! call Rdensity_Writing(Rd,psi%Basis,nio=21,t=ZERO)
      !  STOP 'cc propa'
      ! ---------------------------------- Beging  propagation----------------------------------------------------------
      
      !call Hagedorn_temp(psi, psi0,propa)
      call Calc_Av_E(e0,psi,H)

      DO i = 0, nt
         t = i*propa%delta_t
         t_deltat = t + propa%delta_t
         write (out_unit, *) propa%propa_name2, i, t, t_deltat
          If (propa%propa_name2 == 'vp') Then
          call Get_Basis_Parameters(psi%Basis,Qt,SQt,At,Pt)
          Else
            call  Calc_Basis_parameters(psi,Qt,SQt,At,Pt)
          End if
         call Calc_Av_E(E,psi,H)
          call Calc_Norm_OF_psi(psi, Norm)
          call Population(psi, pop)
          call  Calc_Auto_corr(psi_t0, psi, aut_func, aut_func_arg, propa%propa_name,propa%renorm,t=t)

          write (11,*) t, Qt
          write (12,*) t, E !FMT= "(F20.10,F20.10)"
          write (13,*) t, SQt
          write (14,FMT= *) t, Norm !"(F20.10,F20.10)"
          write (26,*) t,aut_func
          write (18,*) t, pop
           write (19,*) t, Pt
           !write (20,"(F20.10,F10.5,F10.5,F10.5,F10.5)") t, real(At(:),kind=Rkind), aimag(At(:))
          write (20,*) t, real(At(:),kind=Rkind), aimag(At(:)) !t, sqrt(real(At(:),kind=Rkind)),At(:)
          flush  (11);flush  (18)
          flush  (12);flush  (19)
          flush  (13);flush  (20)
          flush  (14)
          flush  (26)
          call eval_pics(psi,ib=28,t=t_deltat)
         if (mod(i, 60 ) == 0) then
              call write_psi(psi=psi, psi_cplx=.false., print_psi_grid=.false. &
              , print_basis=.false., t=t_deltat, int_print=10, real_part=.false.)
             ! call eval_pics(psi,ib=28,t=t)
              call Calc_reduced_density(Rd,psi%CVec,psi%Basis)
              call Rdensity_Writing(Rd,psi%Basis,nio=21,ib=1,t=t_deltat)
             ! call test_propa_Hagedorn(psi,t)
              write(10,*)
              flush  (10)
              write(21,*)
              flush  (21)

         end if
        call march_Global(psi, psi_dt, t_deltat, propa,H)
        If (propa%propa_name == 'hagedorn') Then
           !> cette partitie reconstruit le potentiel apres reconstruction de la base.
            !deallocate(H%Scal_pot)
            !call Set_Op(H, psi%Basis,Tab_Iq)

        End if
  
      END DO
    psif%CVec(:) = psi_dt%CVec(:)
      CALL Calc_Norm_OF_Psi(psif, Norm)         

      IF(debug) then

         write (out_unit, *) 'END propagation'
         write (out_unit, *) 'norm,psi_dt', Norm
         flush (out_unit)

      END IF

   END SUBROUTINE

   SUBROUTINE test_propa_Hagedorn(psi,t)
    TYPE(psi_t), intent(in)        :: psi
     real(kind=Rkind),intent(in)   :: t
    integer                        :: ib,nb,nsurf,ndim
    real(kind=Rkind),allocatable   :: V(:)

    ndim = size(psi%Basis%tab_basis)
    nsurf = psi%Basis%tab_basis(ndim)%nb
    nb = psi%basis%nb*nsurf
    allocate(V(nb))
    V(:) = abs(psi%CVec(:))
    write(22,*) t,abs(ONE-V(1))**2,abs(maxval(V(2:nb)))**2
   deallocate(V)

   End SUBROUTINE


 End module Propa_m
