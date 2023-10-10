!! CARMA sulfate module, as in GEOS, but in a box model.
!!  To be used for training a neural network as a surrogate model
!!  for CARMA. Should have the following inputs:
!!
!!   Temperature T
!!   Pressure    P
!!   H2SO4 MMR   H2SO4
!!   H2O MMR     H2O
!!   SU MMR Bins SU001...SU024
!!   28 input variables
!!
!! @author Parker Case
!! @version Apr-2023

program carma_box
  implicit none

  write(*,*) "Box Test"
  
  call test_box_simple()  
end program

!! Just have one grid box. In that grid box, but an initial concentration
!! of drops at the smallest size, then allow that to grow using a gas. The
!! total mas of drops + gas should be conserved.
subroutine test_box_simple()
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
  use ncio_mod

  implicit none

  integer, parameter        :: NX           = 1
  integer, parameter        :: NY           = 1
  integer, parameter        :: NZ           = 1
  integer, parameter        :: NZP1         = NZ+1
  integer, parameter        :: NELEM        = 1
  integer, parameter        :: NBIN         = 24
  integer, parameter        :: NGROUP       = 1
  integer, parameter        :: NSOLUTE      = 0
  integer, parameter        :: NGAS         = 2
  integer, parameter        :: NWAVE        = 0
  integer, parameter        :: LUNOPRT      = 6
  


  ! Different sizes for time steps provide different results
  ! because of the satbility issues.
  real(kind=f), parameter   :: dtime  = 900._f ! 15 minute time step
  real(kind=f), parameter   :: deltax = 100._f
  real(kind=f), parameter   :: deltay = 100._f
  real(kind=f), parameter   :: deltaz = 100._f
  real(kind=f), parameter   :: zmin   = 0._f
  real(kind=f), parameter   :: secsperday = 24._f * 3600._f

  integer, parameter        :: nstep  = 1 ! 1 time step
  
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
  real(kind=f), allocatable   :: t(:,:,:),t_orig(:,:,:),deltaT(:,:,:)
  real(kind=f), allocatable   :: relhum(:,:,:)
  real(kind=f), allocatable   :: rho(:,:,:)
  
  real(kind=f), allocatable   :: mmr(:,:,:,:,:),mmrtotal(:,:,:)
  real(kind=f), allocatable   :: mmr_gas(:,:,:,:)
  real(kind=f), allocatable   :: new_gas(:,:,:,:)
  real(kind=f), allocatable   :: satliq(:,:,:,:)
  real(kind=f), allocatable   :: satice(:,:,:,:)
  
  real(kind=f), allocatable   :: r(:)
  real(kind=f), allocatable   :: rmass(:)
  real(kind=f), allocatable   :: numberDensity(:,:,:)
  real(kind=f), allocatable   :: r_wet(:,:,:)
  real(kind=f), allocatable   :: rhop_wet(:,:,:)

  real(kind=f), allocatable   :: lat(:,:)
  real(kind=f), allocatable   :: lon(:,:)

  real(kind=f) :: p_0, zc_0, t_0, zl_0, h2o_0, h2so4_0
  real(kind=f) :: su001, su002, su003, su004, su005, su006, su007, su008, su009, su010, su011, su012, su013, su014, su015, su016, su017, su018, su019, su020, su021, su022, su023, su024

  integer               :: outid
  character(len=80)     :: binName(NELEM, NBIN)
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
  integer               :: bins(NBIN)

  real(kind=f)          :: nretries
  real(kind=f)          :: lastret = 0._f

  real(kind=f)          :: time
  real(kind=f)          :: rmin, rmrat,rrat(NBIN)
  real(kind=f)          :: RHO_SULFATE   
  real(kind=f)          :: drh
  real(kind=f)          :: totalnd(NZ,NGROUP)
  real(kind=f)          :: totalre(NZ,NGROUP)
  real(kind=f)          :: totalad(NZ,NGROUP)
  real(kind=f)          :: totalmr(NZ,NGROUP)
  real(kind=f)          :: re2(NZ,NGROUP)
  real(kind=f)          :: re3(NZ,NGROUP)

  ! Allocate the arrays that we need for the model
  allocate(xc(NZ,NY,NX), dx(NZ,NY,NX), yc(NZ,NY,NX), dy(NZ,NY,NX), &
           zc(NZ,NY,NX), zl(NZP1,NY,NX), p(NZ,NY,NX), pl(NZP1,NY,NX), &
           t(NZ,NY,NX),rho(NZ,NY,NX))
  allocate(mmr(NZ,NY,NX,NELEM,NBIN))
  allocate(mmr_gas(NZ,NY,NX,NGAS))
  allocate(new_gas(NZ,NY,NX,NGAS))
  allocate(satliq(NZ,NY,NX,NGAS))
  allocate(satice(NZ,NY,NX,NGAS))
  allocate(r(NBIN))
  allocate(rmass(NBIN))
  allocate(lat(NY,NX), lon(NY,NX))
  allocate(numberDensity(NZ,NY,NX))
  allocate(r_wet(NZ,NY,NX))
  allocate(rhop_wet(NZ,NY,NX))
  allocate(t_orig(NZ,NY,NX),deltaT(NZ,NY,NX))
  allocate(mmrtotal(NZ,NY,NX))

  ! Get command line argument for name of text file
  call GET_COMMAND_ARGUMENT(1, pid_input)

  ! Define the particle-grid extent of the CARMA test
  call CARMA_Create(carma, NBIN, NELEM, NGROUP, NSOLUTE, NGAS, NWAVE, rc, &
      LUNOPRT=6)
  if (rc /=0) stop "    *** CARMA_Create FAILED ***"
  
  carma_ptr => carma

  ! Define the groups
  rmrat = 3.7515201_f
  rmin  = 2.6686863e-8_f
  RHO_SULFATE = 1.923_f  ! dry density of sulfate particles (g/cm3)
  
  call CARMAGROUP_Create(carma, 1, "sulfate", rmin, rmrat, I_SPHERE, 1._f, .false., &
                        rc, irhswell=I_WTPCT_H2SO4, do_drydep=.false., &
                        shortname="SULF", is_sulfate=.true.)
  if (rc /=0) stop "    *** CARMAGROUP_Create FAILED ***"
  
  ! Define the elements
  call CARMAELEMENT_Create(carma, 1, 1, "Sulfate", RHO_SULFATE, I_VOLATILE, I_H2SO4, rc, shortname="SULF")    
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


  call CARMA_Initialize(carma, rc, do_grow=.true., do_coag=.true., do_substep=.true., &
          do_thermo=.true., maxretries=16, maxsubsteps=32, dt_threshold=1._f)
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

  do ibin = 1,NBIN
    bins(ibin) = ibin
  end do

  ! Initial Conditions:
  !
  ! p, T, z, mmrgas, rmin, particle concentration
  ! 90 hPa, 190 K, 17 km, H2O mmr 6.22e-5 g/g, H2SO4 mmr 8.85e-7 g/g 
  open(unit=2, file='box_input_' // trim(pid_input) // '.txt')
  read(2,*) p_0, zc_0, t_0, zl_0, h2o_0, h2so4_0
  read(2,*) su001, su002, su003, su004, su005, su006, su007, su008, su009, su010, su011, su012, su013, su014, su015, su016, su017, su018, su019, su020, su021, su022, su023, su024
  close(2)
  p(1,:,:)         = p_0
  zc(1,:,:)        = zl_0
  t(1,:,:)         = t_0
  zl(1,:,:)        = zl_0 - deltaz
  zl(2,:,:)        = zl_0 + deltaz
  rho(1,:,:)       = (p_0 * 10._f) / (R_AIR * t_0) * (1e-3_f * 1e6_f)
  pl(1,:,:)        = p_0 - (zl(1,:,:) - zc(1,:,:)) * rho(1,:,:) * (GRAV / 100._f)
  pl(2,:,:)        = p_0 - (zl(2,:,:) - zc(1,:,:)) * rho(1,:,:) * (GRAV / 100._f)
  PRINT *, 'P'
  PRINT *, p_0
  PRINT *, 'ZL'
  PRINT *, zl_0
  PRINT *, 'T'
  PRINT *, t_0
  PRINT *, 'H2O'
  PRINT *, h2o_0
  PRINT *, 'H2SO4'
  PRINT *, h2so4_0

  ! Initial H2O and H2SO4 concentrations
  mmr_gas(:,:,:,1)  = h2o_0     ! H2O
  mmr_gas(:,:,:,2)  = h2so4_0     ! H2SO4

  satliq(:,:,:,:)   = -1._f
  satice(:,:,:,:)   = -1._f
  
  ! Initial sulfate concentration
  mmr(:,:,:,:,1)    = su001
  mmr(:,:,:,:,2)    = su002
  mmr(:,:,:,:,3)    = su003
  mmr(:,:,:,:,4)    = su004
  mmr(:,:,:,:,5)    = su005
  mmr(:,:,:,:,6)    = su006
  mmr(:,:,:,:,7)    = su007
  mmr(:,:,:,:,8)    = su008
  mmr(:,:,:,:,9)    = su009
  mmr(:,:,:,:,10)   = su010
  mmr(:,:,:,:,11)   = su011
  mmr(:,:,:,:,12)   = su012
  mmr(:,:,:,:,13)   = su013
  mmr(:,:,:,:,14)   = su014
  mmr(:,:,:,:,15)   = su015
  mmr(:,:,:,:,16)   = su016
  mmr(:,:,:,:,17)   = su017
  mmr(:,:,:,:,18)   = su018
  mmr(:,:,:,:,19)   = su019
  mmr(:,:,:,:,20)   = su020
  mmr(:,:,:,:,21)   = su021
  mmr(:,:,:,:,22)   = su022
  mmr(:,:,:,:,23)   = su023
  mmr(:,:,:,:,24)   = su024
  
  t_orig(1,:,:) = t(1,:,:)
  
  ! Iterate the model over a few time steps.
  do istep = 1, nstep
    ! Create a CARMASTATE for this column.
    do ixy = 1, NX*NY
      ix = ((ixy-1) / NY) + 1
      iy = ixy - (ix-1)*NY
      call CARMASTATE_Create(cstate, carma_ptr, time, dtime, NZ, &
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
        do ibin = 1, NBIN
          call CARMASTATE_SetBin(cstate, ielem, ibin, mmr(:,iy,ix,ielem,ibin), rc)
          if (rc /=0) stop "    *** CARMASTATE_SetBin FAILED ***"
        end do
      end do
      PRINT *, 'MMR'
      PRINT *, mmr(:,:,:,:,:)

      ! Send the gas mmrs to CARMA
      !
      ! For substepping to do anything, during a step, the old an current
      ! gas mmrs or temperatures need to be different.
      
      ! If you want to add some H2SO4, you can do it here using one or the other
      ! of theses lines.
      new_gas(:,iy,ix,:) = mmr_gas(:,iy,ix,:)
      PRINT *, 'MMR GAS'
      PRINT *, mmr_gas(:,:,:,:)

      do igas = 1, NGAS
        call CARMASTATE_SetGas(cstate, igas, new_gas(:,iy,ix,igas), rc, &
                mmr_old=mmr_gas(:,iy,ix,igas),&
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

      call CARMASTATE_GetState(cstate, rc, t=t(:,iy,ix))
      if (rc /=0) stop "    *** CARMASTATE_GetState FAILED ***"

      ! Get the updated bin mmr.
      ! Get the updated bin mmr.
      totalnd(:,:) = 0._f
      totalad(:,:) = 0._f
      totalre(:,:) = 0._f
      re2(:,:) = 0._f
      re3(:,:) = 0._f

      do ielem = 1, NELEM
        call CARMAELEMENT_Get(carma, ielem, rc, igroup=igroup, shortname=sname)

        do ibin = 1, NBIN
          call CARMASTATE_GetBin(cstate, ielem, ibin, mmr(:,iy,ix,ielem,ibin), rc, numberDensity=numberDensity(:,iy,ix), r_wet=r_wet(:,iy,ix), rhop_wet=rhop_wet(:,iy,ix))
          if (rc /=0) stop "    *** CARMASTATE_GetBin FAILED ***"

          if(numberDensity(1,iy,ix) .ne. CAM_FILL)then
            totalnd(:,igroup) = totalnd(:,igroup) + numberDensity(:,iy,ix)
            re2(:,igroup) = re2(:,igroup) + numberDensity(:,iy,ix) * ((r_wet(:,iy,ix)*rrat(ibin))**2)
            re3(:,igroup) = re3(:,igroup) + numberDensity(:,iy,ix) * ((r_wet(:,iy,ix)*rrat(ibin))**3)
            totalad(:,igroup)  = totalad(:,igroup)  + numberDensity(:,iy,ix) * 4.0_f * PI * (r_wet(:,iy,ix)**2) * 1.0e8_f
          end if
        end do
        mmrtotal(:,iy,ix) = sum(mmr(:,iy,ix,ielem,:),DIM=2)
      end do

      !totalre(:,:) = (re3(:,:) / re2(:,:)) * 1e4_f

      ! Get the updated gas mmr.
      do igas = 1, NGAS
        call CARMASTATE_GetGas(cstate, igas, &
        mmr_gas(:,iy,ix,igas), rc, &
        satliq=satliq(:,iy,ix,igas), &
        satice=satice(:,iy,ix,igas))
        if (rc /=0) stop "    *** CARMASTATE_GetGas FAILED ***"
      end do

      lastsub = nsubsteps
      lastret = nretries

      deltaT(:,iy,ix) = t(:,iy,ix) - t_orig(:,iy,ix)

    end do   ! space loop
  end do ! time loop

  h2o_0 = mmr_gas(1,1,1,1)
  h2so4_0 = mmr_gas(1,1,1,2)
  t_0 = t(1,1,1)
  su001 = mmr(1,1,1,1,1)
  su002 = mmr(1,1,1,1,2)
  su003 = mmr(1,1,1,1,3)
  su004 = mmr(1,1,1,1,4)
  su005 = mmr(1,1,1,1,5)
  su006 = mmr(1,1,1,1,6)
  su007 = mmr(1,1,1,1,7)
  su008 = mmr(1,1,1,1,8)
  su009 = mmr(1,1,1,1,9)
  su010 = mmr(1,1,1,1,10)
  su011 = mmr(1,1,1,1,11)
  su012 = mmr(1,1,1,1,12)
  su013 = mmr(1,1,1,1,13)
  su014 = mmr(1,1,1,1,14)
  su015 = mmr(1,1,1,1,15)
  su016 = mmr(1,1,1,1,16)
  su017 = mmr(1,1,1,1,17)
  su018 = mmr(1,1,1,1,18)
  su019 = mmr(1,1,1,1,19)
  su020 = mmr(1,1,1,1,20)
  su021 = mmr(1,1,1,1,21)
  su022 = mmr(1,1,1,1,22)
  su023 = mmr(1,1,1,1,23)
  su024 = mmr(1,1,1,1,24)

  write(pid, '(I0)') getpid()
  open(unit=3, file='box_output_' // trim(pid) // '.txt')
  write(3,*) p_0, zc_0, t_0, zl_0, h2o_0, h2so4_0
  write(3,*) su001, su002, su003, su004, su005, su006, su007, su008, su009, su010, su011, su012, su013, su014, su015, su016, su017, su018, su019, su020, su021, su022, su023, su024
  close(3)

  ! Cleanup the carma state objects
  call CARMASTATE_Destroy(cstate, rc)
  if (rc /=0) stop "    *** CARMASTATE_Destroy FAILED ***"

  call CARMA_Destroy(carma, rc)
  if (rc /=0) stop "    *** CARMA_Destroy FAILED ***"
end subroutine
