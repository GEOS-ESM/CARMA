!! CARMA sulfate module, as in GEOS, but in a box model. To be used for 
!!  running a CARMA box model from a wrapper. Runs one timestep. Should
!!  have the following inputs:
!!
!!   Temperature T
!!   Pressure    P
!!   H2SO4 MMR   H2SO4
!!   H2O MMR     H2O
!!   SU MMR Bins SU001...SU024
!!   28 input variables
!!
!! TODO: 
!!
!! @author Parker Case
!! @version 2023/10/17: Added custom particle properties and variable 
!!                        timesteps
!!          2023/10/10: Updated for f2py setup
!!          2023/04/02: First crack, using files for passing parameters
!!
!! Just have one grid box. Allow for all sulfate processes:
!!   nucleation, condenstation, coagulation, settling.
subroutine carma_box_sulfate(rmrat, rmin, rhop, t_0, p_0, h2so4_0, h2o_0,  mmr_0, dt, nt, constant_h2so4, nbin, mmr_out, t_out, p_out, h2so4_out, h2o_out, rhoa_out, rh_out)
  use carma_precision_mod 
  use carma_constants_mod 
  use carma_enums_mod 
  use carma_types_mod 
  use carmaelement_mod
  use carmagroup_mod
  use carmagas_mod
  use carmastate_mod
  use carma_mod
  use atmosphere_mod

  implicit none

  integer, parameter        :: NX           = 1
  integer, parameter        :: NY           = 1
  integer, parameter        :: NZ           = 1
  integer, parameter        :: NZP1         = NZ+1
  integer, parameter        :: NELEM        = 1
  integer, parameter        :: NGROUP       = 1
  integer, parameter        :: NSOLUTE      = 0
  integer, parameter        :: NGAS         = 2
  integer, parameter        :: NWAVE        = 0
  integer, parameter        :: LUNOPRT      = 6
  
  ! Different sizes for time steps provide different results
  ! because of the satbility issues.
  real(kind=f), parameter   :: deltax = 100._f
  real(kind=f), parameter   :: deltay = 100._f
  real(kind=f), parameter   :: deltaz = 100._f
  real(kind=f), parameter   :: zmin   = 0._f

  integer, parameter        :: I_H2SO4  = 1

  type(carma_type), target            :: carma
  type(carma_type), pointer           :: carma_ptr
  type(carmastate_type)               :: cstate
  integer                             :: rc = 0
  
  real(kind=f), allocatable   :: xc(:,:,:)
  real(kind=f), allocatable   :: dx(:,:,:)
  real(kind=f), allocatable   :: yc(:,:,:)
  real(kind=f), allocatable   :: dy(:,:,:)
  real(kind=f), allocatable   :: zc(:,:,:)
  real(kind=f), allocatable   :: zl(:,:,:)
  real(kind=f), allocatable   :: p(:,:,:)
  real(kind=f), allocatable   :: pl(:,:,:)
  real(kind=f), allocatable   :: t(:,:,:),t_orig(:,:,:)
  real(kind=f), allocatable   :: relhum(:,:,:)
  real(kind=f), allocatable   :: rho(:,:,:)
  
  real(kind=f), allocatable   :: mmr(:,:,:,:,:)
  real(kind=f), allocatable   :: mmr_gas(:,:,:,:)
  real(kind=f), allocatable   :: mmr_gas_old(:,:,:,:)
  real(kind=f), allocatable   :: satliq(:,:,:,:)
  real(kind=f), allocatable   :: satice(:,:,:,:)
  
  real(kind=f), allocatable   :: r(:)
  real(kind=f), allocatable   :: rmass(:)
  real(kind=f), allocatable   :: numberDensity(:,:,:)
  real(kind=f), allocatable   :: r_wet(:,:,:)
  real(kind=f), allocatable   :: rhop_wet(:,:,:)

  real(kind=f), allocatable   :: lat(:,:)
  real(kind=f), allocatable   :: lon(:,:)

  integer      :: nbin, nt, retries
  logical      :: constant_h2so4
  real(kind=f) :: rmrat, rmin, rhop, dt, zc_0, zl_0, dsthresh
  real(kind=f) :: p_0(1), t_0(1), h2o_0(1), h2so4_0(1)
  real(kind=f) :: mmr_0(nbin)
  real(kind=f), intent(out) :: mmr_out(nbin)
  real(kind=f), intent(out) :: p_out(1)
  real(kind=f), intent(out) :: t_out(1)
  real(kind=f), intent(out) :: rhoa_out(1)
  real(kind=f), intent(out) :: rh_out(1)
  real(kind=f), intent(out) :: h2o_out(1)
  real(kind=f), intent(out) :: h2so4_out(1)

  integer               :: outid
  character(len=80)     :: binName(NELEM, nbin)
  character(len=80)     :: gasName(NGAS)
  character(len=80)     :: sname
  character(len=80)     :: groupmassname,grouprname,binndname,binwetrhopname,binwetrname
  character(len=80)     :: satliqName,saticeName,mrname,gname
  character(len=10)      :: pid
  character(len=10)      :: pid_input

  integer               :: i
  integer               :: ix
  integer               :: iy
  integer               :: ixy
  integer               :: istep
  integer               :: igas
  integer               :: igroup
  integer               :: ielem
  integer               :: ibin
  integer               :: nsubsteps
  integer               :: lastsub = 0
  integer               :: date
  integer               :: bins(nbin)

  real(kind=f)          :: nretries
  real(kind=f)          :: lastret = 0._f

  real(kind=f)          :: time
  real(kind=f)          :: rrat(nbin)
  real(kind=f)          :: drh

  ! Allocate the arrays that we need for the model
  allocate(xc(NZ,NY,NX), dx(NZ,NY,NX), yc(NZ,NY,NX), dy(NZ,NY,NX), &
           zc(NZ,NY,NX), zl(NZP1,NY,NX), p(NZ,NY,NX), pl(NZP1,NY,NX), &
           t(NZ,NY,NX),rho(NZ,NY,NX))
  allocate(mmr(NZ,NY,NX,NELEM,nbin))
  allocate(mmr_gas(NZ,NY,NX,NGAS))
  allocate(mmr_gas_old(NZ,NY,NX,NGAS))
  allocate(satliq(NZ,NY,NX,NGAS))
  allocate(satice(NZ,NY,NX,NGAS))
  allocate(r(nbin))
  allocate(rmass(nbin))
  allocate(lat(NY,NX), lon(NY,NX))
  allocate(numberDensity(NZ,NY,NX))
  allocate(r_wet(NZ,NY,NX))
  allocate(rhop_wet(NZ,NY,NX))
  allocate(t_orig(NZ,NY,NX),)

  ! Define the particle-grid extent of the CARMA test
  call CARMA_Create(carma, nbin, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, &
      LUNOPRT=6)
  if (rc /=0) stop "    *** CARMA_Create FAILED ***"
  
  carma_ptr => carma

  call CARMAGROUP_Create(carma, 1, "sulfate", rmin, rmrat, I_SPHERE, 1._f, .false., &
                        rc, irhswell=I_WTPCT_H2SO4, do_drydep=.false., &
                        shortname="SULF", do_vtran=.true., is_sulfate=.true.)
!  call CARMAGROUP_Create(carma, 1, "sulfate", rmin, rmrat, I_SPHERE, 1._f, .false., &
!                        rc, irhswell=I_WTPCT_H2SO4, do_drydep=.false., &
!                        shortname="SULF", is_sulfate=.true.)
  if (rc /=0) stop "    *** CARMAGROUP_Create FAILED ***"
  
  ! Define the elements
  call CARMAELEMENT_Create(carma, 1, 1, "Sulfate", rhop, I_VOLATILE, I_H2SO4, rc, shortname="SULF")    
  if (rc /=0) stop "    *** CARMAELEMENT_Create FAILED ***"
  
  ! Define the gases
  call CARMAGAS_Create(carma, 1, "Water Vapor", WTMOL_H2O, I_VAPRTN_H2O_MURPHY2005, &
    I_GCOMP_H2O, rc, shortname = "H2O", dgc_threshold=0.1_f, ds_threshold=0.1_f)
  if (rc /=0) stop "    *** CARMAGAS_Create FAILED ***"

  call CARMAGAS_Create(carma, 2, "Sulpheric Acid", 98.078479_f, I_VAPRTN_H2SO4_AYERS1980, &
    I_GCOMP_H2SO4, rc, shortname = "H2SO4", dgc_threshold=0.1_f, ds_threshold=0.1_f)
  if (rc /=0) stop "    *** CARMAGAS_Create FAILED ***"

  ! Setup the CARMA processes to exercise
  call CARMA_AddGrowth(carma, 1, 2, rc)   ! set H2SO4 to be the condensing gas
  if (rc /=0) stop "    *** CARMA_AddGrowth FAILED ***"

   call CARMA_AddNucleation(carma, 1, 1, I_HOMNUC, 0._f, rc, igas=2)
  if (rc /=0) stop "    *** CARMA_AddNucleation FAILED ***"

   call CARMA_AddCoagulation(carma, 1, 1, 1, I_COLLEC_FUCHS, rc)
  if (rc /=0) stop "    *** CARMA_AddCoagulation FAILED ***"


  call CARMA_Initialize(carma, rc, do_grow=.true., do_coag=.true., &
          do_substep=.true., do_thermo=.true., maxretries=16, maxsubsteps=32, dt_threshold=1._f)
  if (rc /=0) stop "    *** CARMA_Initialize FAILED ***"
  
  ! For simplicity of setup, do a case with Cartesian coordinates,
  ! which are specified in this interface in meters.
  !
  ! NOTE: For Cartesian coordinates, the first level is the bottom 
  ! of the model (e.g. z = 0), while for sigma and hybrid coordinates
  ! the first level is the top of the model.
  lat = 10.0_f
  lon = 10.0_f

  ! Horizonal centers
  do ix = 1, NX
    do iy = 1, NY
      dx(:,iy,ix) = deltax
      xc(:,iy,ix) = ix*dx(:,iy,ix) / 2._f
      dy(:,iy,ix) = deltay
      yc(:,iy,ix) = iy*dy(:,iy,ix) / 2._f
    end do
  end do

  ! Vertical center
  do i = 1, NZ
    zc(i,:,:) = zmin + (deltaz * (i - 0.5_f))
  end do
  
  call GetStandardAtmosphere(zc, p=p, t=t)

  ! Vertical edge
  do i = 1, NZP1
    zl(i,:,:) = zmin + ((i - 1) * deltaz)
  end do
  call GetStandardAtmosphere(zl, p=pl)

  do ibin = 1,nbin
    bins(ibin) = ibin
  end do

  ! Initial Conditions:
  p(1,:,:)         = p_0(1)
  zc(1,:,:)        = zl_0
  t(1,:,:)         = t_0(1)
  zl(1,:,:)        = zl_0 - deltaz
  zl(2,:,:)        = zl_0 + deltaz
  rho(1,:,:)       = (p_0(1) * 10._f) / (R_AIR * t_0(1)) * (1e-3_f * 1e6_f)
  pl(1,:,:)        = p_0(1) - (zl(1,:,:) - zc(1,:,:)) * rho(1,:,:) * (GRAV / 100._f)
  pl(2,:,:)        = p_0(1) - (zl(2,:,:) - zc(1,:,:)) * rho(1,:,:) * (GRAV / 100._f)

  ! Initial H2O and H2SO4 concentrations
  mmr_gas(:,:,:,1)  = h2o_0(1)     ! H2O
  mmr_gas(:,:,:,2)  = h2so4_0(1)     ! H2SO4
  mmr_gas_old(:,:,:,1) = h2o_0(1)
  mmr_gas_old(:,:,:,2) = h2so4_0(1)

  satliq(:,:,:,:)   = -1._f
  satice(:,:,:,:)   = -1._f
  
  ! Initial sulfate concentration
  do ibin = 1,nbin
    mmr(:,:,:,:,ibin) = mmr_0(ibin)
  end do
  
  t_orig(1,:,:) = t(1,:,:)
  
  ! Iterate the model over a few time steps.
  do istep = 1, nt
    ! Create a CARMASTATE for this column.
    do ixy = 1, NX*NY
      ix = ((ixy-1) / NY) + 1
      iy = ixy - (ix-1)*NY
      call CARMASTATE_Create(cstate, carma_ptr, time, dt, NZ, &
                          I_CART, I_CART, lat(iy,ix), lon(iy,ix), &
                          xc(:,iy,ix), dx(:,iy,ix), &
                          yc(:,iy,ix), dy(:,iy,ix), &
                          zc(:,iy,ix), zl(:,iy,ix), &
                          p(:,iy,ix),  pl(:,iy,ix), &
                          t(:,iy,ix), rc, &
                          told=t(:,iy,ix), &
                          qh2o=mmr_gas(:,iy,ix,1))
      if (rc /=0) stop "    *** CARMASTATE_Create FAILED ***"
      ! Send the bin mmrs to CARMA
      do ielem = 1, NELEM
        do ibin = 1, nbin
          call CARMASTATE_SetBin(cstate, ielem, ibin, mmr(:,iy,ix,ielem,ibin), rc)
          if (rc /=0) stop "    *** CARMASTATE_SetBin FAILED ***"
        end do
      end do

      ! Send the gas mmrs to CARMA
      !
      ! For substepping to do anything, during a step, the old an current
      ! gas mmrs or temperatures need to be different.
      do igas = 1, NGAS
        call CARMASTATE_SetGas(cstate, igas, mmr_gas(:,iy,ix,igas), rc, &
                mmr_old=mmr_gas_old(:,iy,ix,igas),&
                satice_old=satice(:,iy,ix,igas), &
                satliq_old=satliq(:,iy,ix,igas))
        if (rc /=0) stop "    *** CARMASTATE_SetGas FAILED ***"
      end do

      ! Execute the step
      call CARMASTATE_Step(cstate, rc)
      if (rc /=0) stop "    *** CARMASTATE_Step FAILED ***"
       

      ! Get the retry stats and the updated temperature.
      call CARMASTATE_Get(cstate, rc, nsubstep=nsubsteps, nretry=nretries)
      if (rc /=0) stop "    *** CARMASTATE_Get FAILED ***"

      call CARMASTATE_GetState(cstate, rc, rhoa_wet=rhoa_out, t=t(:,iy,ix))
      if (rc /=0) stop "    *** CARMASTATE_GetState FAILED ***"

      ! Get the updated bin mmr.
      do ielem = 1, NELEM
        call CARMAELEMENT_Get(carma, ielem, rc, igroup=igroup, shortname=sname)

        do ibin = 1, nbin
          call CARMASTATE_GetBin(cstate, ielem, ibin, mmr(:,iy,ix,ielem,ibin), rc, numberDensity=numberDensity(:,iy,ix), r_wet=r_wet(:,iy,ix), rhop_wet=rhop_wet(:,iy,ix))
          if (rc /=0) stop "    *** CARMASTATE_GetBin FAILED ***"

        end do
      end do


      ! Get the updated gas mmr.
      do igas = 1, NGAS
        mmr_gas_old(:,iy,ix,igas) = mmr_gas(:,iy,ix,igas)
        call CARMASTATE_GetGas(cstate, igas, &
        mmr_gas(:,iy,ix,igas), rc, &
        satliq=satliq(:,iy,ix,igas), &
        satice=satice(:,iy,ix,igas))
        if (rc /=0) stop "    *** CARMASTATE_GetGas FAILED ***"
      end do

      ! Replace H2SO4 if constant_h2so4
      if (constant_h2so4) then
          mmr_gas(:,:,:,2)  = h2so4_0(1)     ! H2SO4
      end if

      lastsub = nsubsteps
      lastret = nretries

    end do   ! space loop
  end do ! time loop

  h2o_out(1) = mmr_gas(1,1,1,1)
  rh_out(1) = satliq(1,1,1,1)
  h2so4_out(1) = mmr_gas(1,1,1,2)
  t_out(1) = t(1,1,1)
  p_out(1) = p(1,1,1)
  do ibin = 1,nbin
    mmr_out(ibin) = mmr(1,1,1,1,ibin)
  end do

  ! Cleanup the carma state objects
  call CARMASTATE_Destroy(cstate, rc)
  if (rc /=0) stop "    *** CARMASTATE_Destroy FAILED ***"

  call CARMA_Destroy(carma, rc)
  if (rc /=0) stop "    *** CARMA_Destroy FAILED ***"
end subroutine
