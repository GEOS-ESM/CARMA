#include "MAPL_Generic.h"

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !MODULE: CARMAchem_GridCompMod - The Community Aerosol and Radiation Model
!                                  for Atmospheres
!
! !INTERFACE:
!
   MODULE CARMAchem_GridCompMod
!
! !USES:
!
   USE ESMF
   USE MAPL
   USE Chem_Mod 	                        ! Chemistry Base Class

   USE CARMA_GridCompMod                        ! ESMF parent component
   USE CARMA_UtilMod
   USE Chem_UtilMod
   USE m_inpak90	                        ! Resource file management
   use m_die, only: die

!  CARMA Specific Methods
   use carma_precision_mod 
   use carma_constants_mod 
   use carma_enums_mod 
   use carma_types_mod 
   use carmaelement_mod
   use carmagroup_mod
   use carmastate_mod
   use carma_mod

   IMPLICIT NONE
   PRIVATE

   type(Chem_Mie), dimension(2), save :: carmaMieTable
   integer, parameter :: instanceComputational = 1
   integer, parameter :: instanceData          = 2

!
! !PUBLIC MEMBER FUNCTIONS:

   PUBLIC SetServices
!
! !DESCRIPTION: 
!
!  {\tt CARMAchem_GridComp} is an ESMF gridded component for the Community
!  Aerosol and Radiation Model for Atmospheres aerosol and cloud
!  microphysics packages.
!
!  Developed for GEOS-5 release Eros-beta7p6 and later.
!
! !REVISION HISTORY:
!
!  31Jul2006  da Silva  Created the GMI stub.
!  11Dec2007  Nielsen   Real code for Eros-beta7p17.
!  18May2009  Colarco   Developed based on GMIchem_GridComp.F90
!
!EOP
!-------------------------------------------------------------------------

  TYPE CARMAchem_State
     PRIVATE
     TYPE(Chem_Registry), POINTER      :: chemReg   => null()
     TYPE(CARMA_GridComp), POINTER     :: gcCARMA   => null()
     TYPE(CARMA_Registry), POINTER     :: CARMAReg  => null()
     TYPE(Chem_Array), POINTER         :: qa(:)     => null()
  END TYPE CARMAchem_State

  TYPE CARMAchem_WRAP
     TYPE (CARMAchem_State), pointer :: PTR => null()
  END TYPE CARMAchem_WRAP

CONTAINS


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 910.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE: SetServices --- Sets IRF services for CARMA Grid Component
!
! !INTERFACE:

   SUBROUTINE SetServices ( GC, RC )

! !ARGUMENTS:

    type(ESMF_GridComp)                :: GC  ! gridded component
    integer, intent(OUT)               :: RC  ! return code

! !DESCRIPTION: Sets Initialize, Run and Finalize services. 
!
! !REVISION HISTORY:
!
!  31Jul2006  da Silva  First crack.
!  18May2009  Colarco   Adapted for CARMA
!
!EOP
!-------------------------------------------------------------------------

!   ErrLog Variables
!   ----------------
    character(len=ESMF_MAXSTR)      :: IAm = 'SetServices'
    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: COMP_NAME

!   Local derived type aliases
!   --------------------------
    type (ESMF_Config)                :: CF
    type (CARMAchem_State), pointer   :: state   ! internal, that is
    type (CARMAchem_wrap)             :: wrap
    type(CARMA_Registry), pointer     :: r

    integer                           :: n, i_XX, j_XX, i, j, iq
    character(len=ESMF_MAXSTR)        :: binstr

!                              ------------

!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
    VERIFY_(STATUS)
    Iam = TRIM(COMP_NAME) // '::' // TRIM(Iam)

!   Wrap internal state for storing in GC; rename legacyState
!   -------------------------------------
    allocate ( state, stat=STATUS )
    VERIFY_(STATUS)
    wrap%ptr => state
 
!   Start by loading the CARMA Registry
!   -----------------------------------
    allocate ( state%CARMAReg )
    call registry_ ( state%CARMAReg )
    VERIFY_(STATUS)
    call registry_print_ ( state%CARMAReg )

!   Start by loading the Chem Registry
!   ----------------------------------
    allocate ( state%chemReg )
    state%chemReg = Chem_RegistryCreate ( STATUS )
    VERIFY_(STATUS)


    r => state%CARMAReg   ! short hand

!                       ------------------------
!                       ESMF Functional Services
!                       ------------------------

!   Set the Initialize, Run, Finalize entry points
!   ----------------------------------------------
     if ( r%doing_CARMA ) then

        if (MAPL_AM_I_ROOT()) print *, trim(Iam)//': ACTIVE'

        call MAPL_GridCompSetEntryPoint ( GC, ESMF_METHOD_INITIALIZE,  Initialize_, &
             RC=STATUS)
        VERIFY_(STATUS)
        
        call MAPL_GridCompSetEntryPoint ( GC,  ESMF_METHOD_RUN,  Run_,        &
             RC=STATUS)
        VERIFY_(STATUS)
        
        call MAPL_GridCompSetEntryPoint ( GC,  ESMF_METHOD_FINALIZE,  Finalize_,  &
             RC=STATUS)
        VERIFY_(STATUS)
        
!       Store internal state in GC
!       --------------------------
        call ESMF_UserCompSetInternalState ( GC, 'CARMA_state', wrap, STATUS )
        VERIFY_(STATUS)

     else

        if (MAPL_AM_I_ROOT()) &
        print *, trim(Iam)//': NOT ACTIVE, defaulting to Generic No-op stubs'

     endif
  
!                         ------------------
!                         GEOS Data Services
!                         ------------------
!
! !IMPORT STATE:
#include "CARMAchem_Import___.h"

! !INTERNAL STATE:

!
!  NOTES: 
!    1)  vtitle as it stands is as the CF definition of long name.
!        I may need to create a "standard name" in chemReg and pass
!        this to GEOS Generic
!    2)  Host model MUST provide convective transport as well
!
!    Convention for tracer names is CARMA::groupname::elemname::XXX
!    where XXX is the bin number.
!    Will need to add space for the NGAS array
!
     if ( r%doing_CARMA ) then

!     Add particle tracers
      do j = 1, r%NELEM
       do i = 1, r%NBIN
        write(binstr,'(i3)') i
        binstr = adjustl(binstr)
        if(i .lt. 10)  binstr = '0'//binstr
        if(i .lt. 100) binstr = '0'//binstr

        iq = i + (j-1)*r%NBIN
        r%vname(iq) = trim(r%groupname(r%igroup(j)))// '::' &
                   // trim(r%elemname(j))// '::'            &
                   // trim(binstr)
        call MAPL_AddInternalSpec(GC,                              &
               SHORT_NAME  = trim(COMP_NAME)// '::'                &
                          // trim(r%groupname(r%igroup(j)))// '::' &
                          // trim(r%elemname(j))// '::'            &
                          // trim(binstr),                         &
               LONG_NAME   = trim(COMP_NAME)// '::'                &
                          // trim(r%groupname(r%igroup(j)))// '::' &
                          // trim(r%elemname(j))// '::'            &
                          // trim(binstr),                         &
               UNITS       = 'kg/kg',                              &     ! placeholder
               FRIENDLYTO  = 'DYNAMICS:TURBULENCE',          &
               DIMS        = MAPL_DimsHorzVert,                    &
               VLOCATION   = MAPL_VLocationCenter,     RC=STATUS  )
          VERIFY_(STATUS)

       end do
      end do
       
!     Add gas tracers
      do j = 1, r%NGAS

        iq = j + r%NBIN*r%NELEM
        r%vname(iq) = trim(r%gasname(j))

        call MAPL_AddInternalSpec(GC,                              &
               SHORT_NAME  = trim(COMP_NAME)// '::'                &
                          // trim(r%gasname(j)),                   &
               LONG_NAME   = trim(COMP_NAME)// '::'                &
                          // trim(r%gasname(j)),                   &
               UNITS       = 'kg/kg',                              &     ! placeholder
               FRIENDLYTO  = 'DYNAMICS:TURBULENCE',                &
               DIMS        = MAPL_DimsHorzVert,                    &
               VLOCATION   = MAPL_VLocationCenter,     RC=STATUS  )
          VERIFY_(STATUS)

      end do

!     Non-advected tracers
!     Add the old temperature
      iq = r%NBIN*r%NELEM + r%NGAS + 1
      r%vname(iq) = 't_old'

      call MAPL_AddInternalSpec(GC,                              &
             SHORT_NAME  = trim(COMP_NAME)// '::'                &
                        // trim(r%vname(iq)),                    &
             LONG_NAME   = trim(COMP_NAME)// '::'                &
                        // trim(r%vname(iq)),                    &
             UNITS       = 'K',                                  &     ! placeholder
             ADD2EXPORT      = .TRUE.,                           &
             DIMS        = MAPL_DimsHorzVert,                    &
             VLOCATION   = MAPL_VLocationCenter,     RC=STATUS  )
      VERIFY_(STATUS)


!     Add the prior time step gas tracers
      if(r%NGAS > 0) then
      do j = 1, r%NGAS

        iq = j + r%NBIN*r%NELEM + r%NGAS + 1
        r%vname(iq) = trim(r%gasname(j))//'_old'

        call MAPL_AddInternalSpec(GC,                              &
               SHORT_NAME  = trim(COMP_NAME)// '::'                &
                          // trim(r%vname(iq)),                     &
               LONG_NAME   = trim(COMP_NAME)// '::'                &
                          // trim(r%vname(iq)),                     &
               UNITS       = 'kg/kg',                              &     ! placeholder
               ADD2EXPORT      = .TRUE.,                           &
               DIMS        = MAPL_DimsHorzVert,                    &
               VLOCATION   = MAPL_VLocationCenter,     RC=STATUS  )
          VERIFY_(STATUS)

      end do
       
!     Add the prior time step saturation wrt liquid
      do j = 1, r%NGAS

        iq = j + r%NBIN*r%NELEM + r%NGAS + r%NGAS + 1
        r%vname(iq) = 'satliq_'//trim(r%gasname(j))//'_old'
        call MAPL_AddInternalSpec(GC,                              &
               SHORT_NAME  = trim(COMP_NAME)// '::'                &
                          // trim(r%vname(iq)),                     &
               LONG_NAME   = trim(COMP_NAME)// '::'                &
                          // trim(r%vname(iq)),                     &
               UNITS       = '1',                                  &     ! placeholder
               ADD2EXPORT      = .TRUE.,                           &
               DIMS        = MAPL_DimsHorzVert,                    &
               VLOCATION   = MAPL_VLocationCenter,     RC=STATUS  )
          VERIFY_(STATUS)

      end do
       
!     Add the prior time step saturation wrt ice
      do j = 1, r%NGAS

        iq = j + r%NBIN*r%NELEM + r%NGAS + r%NGAS + r%NGAS + 1
        r%vname(iq) = 'satice_'//trim(r%gasname(j))//'_old'

        call MAPL_AddInternalSpec(GC,                              &
               SHORT_NAME  = trim(COMP_NAME)// '::'                &
                          // trim(r%vname(iq)),                     &
               LONG_NAME   = trim(COMP_NAME)// '::'                &
                          // trim(r%vname(iq)),                     &
               UNITS       = '1',                                  &     ! placeholder
               ADD2EXPORT      = .TRUE.,                           &
               DIMS        = MAPL_DimsHorzVert,                    &
               VLOCATION   = MAPL_VLocationCenter,     RC=STATUS  )
          VERIFY_(STATUS)

      end do
      endif  ! NGAS > 0
       
     endif

! !EXPORT STATE:
!   This bundle is needed by radiation - It will contain the 
!   basically the same as the internal state for aerosols
!   and aerosol optics
!   --------------------------------------------------------
    call MAPL_AddExportSpec(GC,                                 &
        SHORT_NAME         = 'AERO',                              &
        LONG_NAME          = 'aerosol_mass_mixing_ratios',        &
        UNITS              = 'kg kg-1',                           &
        DIMS               = MAPL_DimsHorzVert,                   &
        VLOCATION          = MAPL_VLocationCenter,                &
        DATATYPE           = MAPL_StateItem,                      &
                                                       RC=STATUS  )
    VERIFY_(STATUS)

!   This state is needed by MOIST - It will contain aerosols
!   This bundle is not currently filled in by CARMA, just a 
!   place holder for symmetry with GOCART
!   --------------------------------------------------------
    call MAPL_AddExportSpec(GC,                    &
       SHORT_NAME = 'AERO_ACI',                    &
       LONG_NAME  = 'aerosol_cloud_interaction',   &
       UNITS      = 'kg kg-1',                     &
       DIMS       = MAPL_DimsHorzVert,             &
       VLOCATION  = MAPL_VLocationCenter,          &
       DATATYPE   = MAPL_StateItem, __RC__)

!   This bundle is needed by surface for snow albedo 
!   modification by aerosol settling and deposition
!   This bundle is not currently filled in by CARMA, just a 
!   place holder for symmetry with GOCART
!   --------------------------------------------------------
    call MAPL_AddExportSpec(GC,                                   &
        SHORT_NAME         = 'AERO_DP',                           &
        LONG_NAME          = 'aerosol_deposition',                &
        UNITS              = 'kg m-2 s-1',                        &
        DIMS               = MAPL_DimsHorzOnly,                   &
        DATATYPE           = MAPL_BundleItem,                     &
                                                       RC=STATUS  )
   VERIFY_(STATUS)

!
#include "CARMAchem_Export___.h"


!   Set the profiling timers
!   ------------------------
    CALL MAPL_TimerAdd(GC, NAME="INITIALIZE", RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_TimerAdd(GC, NAME="RUN", RC=STATUS)
    VERIFY_(STATUS)
    CALL MAPL_TimerAdd(GC, NAME="FINALIZE", RC=STATUS)
    VERIFY_(STATUS)

!   Generic Set Services
!   --------------------
    call MAPL_GenericSetServices ( GC, RC=STATUS )
    VERIFY_(STATUS)

!   All done
!   --------

    RETURN_(ESMF_SUCCESS)

  END SUBROUTINE SetServices


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Initialize_ --- Initialize CARMA
!
! !INTERFACE:
!

   SUBROUTINE Initialize_ ( gc, import, export, clock, rc )

! !USES:

   implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: clock      ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: gc      ! Grid Component
   type(ESMF_State), intent(inout) :: import     ! Import State
   type(ESMF_State), intent(inout) :: export     ! Export State
   integer, intent(out) ::  rc                    ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  27Feb2005 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------


!  ErrLog Variables
!  ----------------
   character(len=ESMF_MAXSTR)      :: IAm = 'Initialize_'
   integer                         :: STATUS
   character(len=ESMF_MAXSTR)      :: COMP_NAME

   type(CARMA_GridComp), pointer   :: gcCARMA     ! Grid Component
   integer                         :: nymd, nhms  ! time of day
   real                            :: cdt         ! chemistry timestep (secs)

   type(ESMF_Grid)                 :: grid        
   type(ESMF_Config)               :: CF, carmaCF
 
   integer                         :: i1=1, i2, ig=0, im  ! dist grid indices
   integer                         :: j1=1, j2, jg=0, jm  ! dist grid indices
   integer                         :: km, nq              ! dist grid indices
   integer                         :: n, dims(3), l

   type(Chem_Array), pointer       :: qa(:)       ! array of pointers
   type(MAPL_MetaComp), pointer    :: ggState      ! GEOS Generic State
   type(ESMF_State)                :: internal
   type(ESMF_Field)                :: field
   type(ESMF_Field)                :: fld
   type(ESMF_FieldBundle)          :: bundle
   type(ESMF_State)                :: aero
   type(ESMF_FieldBundle)          :: aero_state_aerosols
   character(len=ESMF_MAXSTR)      :: fld_name
   integer                         :: n_aerosols
   type(MAPL_VarSpec), pointer     :: InternalSpec(:)
   integer                         :: instance

   character(len=ESMF_MAXSTR)      :: short_name
   type(CARMA_Registry), pointer   :: reg => null()
   type(Chem_Registry), pointer    :: mieReg => null()  ! pointer to registry for CARMA aero_provider
   integer :: nCARMAbegin, nCARMAend, ibin, ielem, igroup, igas
   real    :: fscav, maxclean, CCNtuning
   character(LEN=ESMF_MAXSTR)      :: CLDMICRO
   
   integer :: i, j, k, iq, istart, iend
   real, parameter :: rad2deg = 180. / MAPL_PI

   character(len=ESMF_MAXSTR), allocatable :: aero_aci_modes(:)
   real(MAPL_R4), allocatable :: aero_rmin(:), aero_rmrat(:)

!  Declare pointers to IMPORT/EXPORT/INTERNAL states 
!  -------------------------------------------------
#  include "CARMAchem_DeclarePointer___.h"

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, CONFIG=CF, __RC__ )

   Iam = trim(COMP_NAME) // '::Initialize_'

   if (MAPL_AM_I_ROOT()) then
      PRINT *, TRIM(Iam)//': Starting...'
      PRINT *,' '
   end if

!  Get my internal MAPL_Generic state
!  -----------------------------------
   call MAPL_GetObjectFromGC ( GC, ggState, __RC__)

   call MAPL_TimerOn(ggState, 'TOTAL')
   call MAPL_TimerOn(ggState, 'INITIALIZE')

!  Initialize GEOS Generic
!  ------------------------
   call MAPL_GenericInitialize ( gc, import, export, clock,  __RC__ )

!  Get parameters from gc and clock
!  --------------------------------
   call extract_ ( gc, clock, gcCARMA, qa, nymd, nhms, cdt, STATUS )
   VERIFY_(STATUS)

!  Get the grid
!  ------------
   call ESMF_GridCompGet ( GC, grid=grid, __RC__ )

   call MAPL_GridGet ( grid, globalCellCountPerDim=DIMS, __RC__ )

   gcCARMA%im = dims(1)
   gcCARMA%jm = dims(2)

   call ESMF_GridGet(GRID, localDE=0, &
        staggerloc=ESMF_STAGGERLOC_CENTER, &
        computationalCount=DIMS, __RC__ )

   gcCARMA%grid = grid
   reg => gcCARMA%carmaReg

!  Get pointers to IMPORT/EXPORT/INTERNAL states 
!  ---------------------------------------------
#  include "CARMAchem_GetPointer___.h"

!  Associate the Internal State fields with our legacy state 
!  ---------------------------------------------------------
   call MAPL_Get ( ggSTATE, INTERNALSPEC=InternalSpec, &
                            INTERNAL_ESMF_STATE=internal, &
                            LONS=gcCARMA%LONS, &
                            LATS=gcCARMA%LATS, __RC__ )

! Convert radians to degrees
  gcCARMA%lons = gcCARMA%lons * rad2deg
  gcCARMA%lats = gcCARMA%lats * rad2deg

! Local sizes of three dimensions
!--------------------------------
   gcCARMA%i2 = dims(1)
   gcCARMA%j2 = dims(2)
   gcCARMA%km = dims(3)

!   Initialize the tracer array
!   ---------------------------
    _ASSERT( size(InternalSpec) == reg%nq, 'needs informative message' )

    do L = 1, size(InternalSpec)

      call MAPL_VarSpecGet(InternalSpec(L), SHORT_NAME=short_name, __RC__)

      call MAPL_GetPointer(internal,NAME=short_name,ptr=qa(L)%data3d, __RC__ )

    end do

!   Bootstrapping option
!   --------------------
!   For sub-stepping of gases we need to provide the prior time-step
!   temperature, gas mmr, and saturation ratios for liquid and ice.
!   If these are not set from the carma_internal_rst we need to 
!   initialize.

!   Prior time-step temperature
    iq = reg%NBIN*reg%NELEM + reg%NGAS + 1
    if(qa(iq)%data3d(gcCARMA%i2,gcCARMA%j2,gcCARMA%km) < 0.) qa(iq)%data3d = -1.

!   Prior time-step gases
    istart = reg%NBIN*reg%NELEM + reg%NGAS + 1 + 1
    iend   = istart + reg%NGAS - 1
    do iq = istart, iend
     if(qa(iq)%data3d(gcCARMA%i2,gcCARMA%j2,gcCARMA%km) < 0.) qa(iq)%data3d = qa(iq-reg%NGAS-1)%data3d
    enddo

!   Prior time-step saturation wrt liquid
    istart = reg%NBIN*reg%NELEM + reg%NGAS + 1 + reg%NGAS + 1
    iend   = istart + reg%NGAS - 1
    do iq = istart, iend
     if(qa(iq)%data3d(gcCARMA%i2,gcCARMA%j2,gcCARMA%km) < 0.) qa(iq)%data3d = -1.
    enddo

!   Prior time-step saturation wrt ice
    istart = reg%NBIN*reg%NELEM + reg%NGAS + 1 + reg%NGAS + reg%NGAS + 1
    iend   = istart + reg%NGAS - 1
    do iq = istart, iend
     if(qa(iq)%data3d(gcCARMA%i2,gcCARMA%j2,gcCARMA%km) < 0.) qa(iq)%data3d = -1.
    enddo

!   Call initialize
!   ---------------
    call CARMA_GridCompInitialize ( gcCARMA, import, export, nymd, nhms, cdt, &
                                    STATUS )
    VERIFY_(STATUS)

#ifdef PRINT_STATES

   if (MAPL_AM_I_ROOT()) then
       print *, trim(Iam)//': INTERNAL State during Initialize():'
       call ESMF_StatePrint ( internal )
       print *, trim(Iam)//': IMPORT   State during Initialize():'
       call ESMF_StatePrint ( import  )
       print *, trim(Iam)//': EXPORT   State during Initialize():' 
       call ESMF_StatePrint ( export  )
    end if

#endif

!   Get the chemistry coupling information from the configuration
!   -------------------------------------------------------------
    call ESMF_ConfigGetAttribute(CF, reg%sulfuric_acid_source, &
         Label="SULFURIC_ACID_SOURCE:" , DEFAULT='tendency', __RC__)
  

!   Fill in the scavenging attribute
!   --------------------------------
    nCARMABegin =  1
    nCARMAEnd   =  reg%nq
    do ielem = 1, reg%NELEM
     igroup = reg%igroup(ielem)
     do ibin = 1, reg%NBIN
      n = nCARMAbegin + (ielem-1)*reg%NBIN + ibin - 1

      call MAPL_VarSpecGet(InternalSpec(n), SHORT_NAME=short_name, __RC__ )
      call ESMF_StateGet(internal, short_name, field, __RC__ )
      fscav = reg%fscav(igroup)
      call ESMF_AttributeSet(field,NAME="ScavengingFractionPerKm",VALUE=fscav, __RC__ )

     end do
    end do

!   Get the Mie tables (GEOS-5 like Mie tables)
!   -------------------------------------------
    call CARMA_GetMieTables(gcCARMA, rc)
    if(rc /= 0) then
     if(MAPL_AM_I_ROOT()) print *, 'CARMA: Failed reading Mie tables'
     RETURN_(ESMF_FAILURE)
    endif

!   Fill the AERO bundle - For now we add all concentration elements
!   --------------------
    call ESMF_StateGet(export, 'AERO', aero, __RC__ )

    aero_state_aerosols = ESMF_FieldBundleCreate(name='AEROSOLS', __RC__)
    call MAPL_StateAdd(aero, aero_state_aerosols, __RC__)

    do ielem = 1, reg%NELEM
     igroup = reg%igroup(ielem)
     if(ielem /= gcCARMA%carma%f_group(igroup)%f_ienconc ) cycle
     do ibin = 1, reg%NBIN
      n = nCARMAbegin + (ielem-1)*reg%NBIN + ibin - 1
      call ESMF_StateGet ( INTERNAL,                     &
                           trim(COMP_NAME) // '::'//     &
                           trim(reg%vname(n)),           &
                           FIELD, __RC__ )
      fld = MAPL_FieldCreate(FIELD, name=reg%vname(n), __RC__)
      call MAPL_FieldBundleAdd(aero_state_aerosols, fld, __RC__)
     end do
    end do

    call ESMF_FieldBundleGet(aero_state_aerosols, fieldCount=n_aerosols, __RC__)

    if (n_aerosols > 0) then
        
!        This is a placeholder code in case sometime I want CARMA to be "data_driven"
!        if (myState%data_driven) then
!            instance = instanceData
!        else
            instance = instanceComputational
!        end if

        carmaCF = ESMF_ConfigCreate(__RC__)
        call ESMF_ConfigLoadFile(carmaCF,'CARMAchem_Registry.rc', __RC__)
        allocate(mieReg, stat=STATUS)
        VERIFY_(STATUS)
        mieReg = Chem_RegistryCreate(rc,rcfile='CARMAchem_MieRegistry.rc')
        if ( rc /= 0 ) call die('CARMA', 'Cannot read CARMAchem_MieRegistry.rc' )
        carmaMieTable(instance) = Chem_MieCreate(carmaCF, chemReg=mieReg, __RC__)
        deallocate(mieReg,stat=STATUS)
        VERIFY_(STATUS)
        call ESMF_ConfigDestroy(carmaCF, __RC__)

        ! Mie Table instance/index
        call ESMF_AttributeSet(aero, name='mie_table_instance', value=instance, __RC__)

        ! state of the atmosphere
        call ESMF_AttributeSet(aero, name='air_pressure_for_aerosol_optics',             value='PLE', __RC__)
        call ESMF_AttributeSet(aero, name='relative_humidity_for_aerosol_optics',        value='RH',  __RC__)
        call ESMF_AttributeSet(aero, name='cloud_area_fraction_for_aerosol_optics',      value='',    __RC__) ! 'cloud_area_fraction_in_atmosphere_layer_for_aerosol_optics'

        ! aerosol optics
        call ESMF_AttributeSet(aero, name='implements_aerosol_optics_method', value=.true., __RC__)

        call ESMF_AttributeSet(aero, name='band_for_aerosol_optics',                     value=0,     __RC__)
        call ESMF_AttributeSet(aero, name='wavelength_for_aerosol_optics',               value=0,     __RC__)
        call ESMF_AttributeSet(aero, name='extinction_in_air_due_to_ambient_aerosol',    value='EXT', __RC__)
        call ESMF_AttributeSet(aero, name='single_scattering_albedo_of_ambient_aerosol', value='SSA', __RC__)
        call ESMF_AttributeSet(aero, name='asymmetry_parameter_of_ambient_aerosol',      value='ASY', __RC__)

        call ESMF_AttributeSet(aero, name='aerosolName', value='',      __RC__)
        call ESMF_AttributeSet(aero, name='im',          value=dims(1), __RC__)
        call ESMF_AttributeSet(aero, name='jm',          value=dims(2), __RC__)
        call ESMF_AttributeSet(aero, name='km',          value=dims(3), __RC__)

        call add_aero (aero, label='air_pressure_for_aerosol_optics', label2='PLE', grid=grid, typekind=MAPL_R4, __RC__)
        call add_aero (aero, label='relative_humidity_for_aerosol_optics', label2='RH', grid=grid, typekind=MAPL_R4, __RC__)
        call add_aero (aero, label='extinction_in_air_due_to_ambient_aerosol', label2='EXT', grid=grid, typekind=MAPL_R4, __RC__)
        call add_aero (aero, label='single_scattering_albedo_of_ambient_aerosol', label2='SSA', grid=grid, typekind=MAPL_R4, __RC__)
        call add_aero (aero, label='asymmetry_parameter_of_ambient_aerosol', label2='ASY', grid=grid, typekind=MAPL_R4, __RC__)
        call add_aero (aero, label='monochromatic_extinction_in_air_due_to_ambient_aerosol', label2='monochromatic_EXT', grid=grid, typekind=MAPL_R4, __RC__)

        call ESMF_MethodAdd(aero, label='run_aerosol_optics', userRoutine=run_aerosol_optics, __RC__)

        ! aerosol cloud interaction
        ! PAC: Not the prettiest way to build a list of aerosol modes
        allocate(aero_aci_modes(n_aerosols))
        allocate(aero_rmin(n_aerosols))
        allocate(aero_rmrat(n_aerosols))
        n = 1
        do ielem = 1, reg%NELEM
         igroup = reg%igroup(ielem)
         if(ielem /= gcCARMA%carma%f_group(igroup)%f_ienconc ) cycle
         do ibin = 1, reg%NBIN
          aero_aci_modes(n) = reg%vname((ielem-1)*reg%NBIN + ibin)
          aero_rmin(n) = reg%rmin(igroup)
          aero_rmrat(n) = reg%rmrat(igroup)
          n = n+1
         end do
        end do
        call ESMF_AttributeSet(aero, name='number_of_aerosol_modes', value=n_aerosols, __RC__)
        call ESMF_AttributeSet(aero, name='aerosol_modes', itemcount=n_aerosols, valuelist=aero_aci_modes, __RC__)
        call ESMF_AttributeSet(aero, name='aerosol_rmin', itemcount=n_aerosols, valuelist=aero_rmin, __RC__)
        call ESMF_AttributeSet(aero, name='aerosol_rmrat', itemcount=n_aerosols, valuelist=aero_rmrat, __RC__)
        deallocate(aero_aci_modes)
        deallocate(aero_rmin)
        deallocate(aero_rmrat)
        call ESMF_ConfigGetAttribute(CF, maxclean, default=1.0e-9, label='MAXCLEAN:', __RC__)
        call ESMF_AttributeSet(aero, name='max_q_clean', value=maxclean, __RC__)
        call ESMF_ConfigGetAttribute(CF, CCNtuning, default=1.8, label='CCNTUNING:', __RC__)
        call ESMF_AttributeSet(aero, name='ccn_tuning', value=CCNtuning, __RC__)
        call ESMF_ConfigGetAttribute( CF, CLDMICRO, Label='CLDMICR_OPTION:',  default="BACM_1M", RC=STATUS)
        call ESMF_AttributeSet(aero, name='cldmicro', value=CLDMICRO, __RC__)


        call add_aero (aero, label='air_temperature', label2='T', grid=grid, typekind=MAPL_R4, __RC__)
        call add_aero (aero, label='fraction_of_land_type', label2='FRLAND', grid=grid, typekind=MAPL_R4, __RC__)
        call add_aero (aero, label='width_of_aerosol_mode', label2='SIGMA', grid=grid, typekind=MAPL_R4, __RC__)
        call add_aero (aero, label='aerosol_number_concentration', label2='NUM', grid=grid, typekind=MAPL_R4, __RC__)
        call add_aero (aero, label='aerosol_dry_size', label2='DGN', grid=grid, typekind=MAPL_R4, __RC__)
        call add_aero (aero, label='aerosol_density', label2='density', grid=grid, typekind=MAPL_R4, __RC__)
        call add_aero (aero, label='aerosol_hygroscopicity', label2='KAPPA', grid=grid, typekind=MAPL_R4, __RC__)
        call add_aero (aero, label='fraction_of_dust_aerosol', label2='FDUST', grid=grid, typekind=MAPL_R4, __RC__)
        call add_aero (aero, label='fraction_of_soot_aerosol', label2='FSOOT', grid=grid, typekind=MAPL_R4, __RC__)
        call add_aero (aero, label='fraction_of_organic_aerosol', label2='FORGANIC', grid=grid, typekind=MAPL_R4, __RC__)

        call ESMF_MethodAdd(aero, label='aerosol_activation_properties', userRoutine=aerosol_activation_properties, __RC__)
       

    end if

#ifdef PRINT_STATES

   if (MAPL_AM_I_ROOT()) then
       print *, trim(Iam)//': AERO Bundle during Initialize():' 
       call ESMF_FieldBundlePrint ( bundle )
   end if

#endif

!  Stop timers
!  -----------
   CALL MAPL_TimerOff(ggState, "INITIALIZE")
   CALL MAPL_TimerOff(ggState, "TOTAL")


    RETURN_(ESMF_SUCCESS)

   END SUBROUTINE Initialize_

!====================================================================================
   subroutine add_aero (state, label, label2, grid, typekind, ptr, rc)

!    Description: Adds fields to aero state for aerosol optics calcualtions. 
!                    Original version in GOCART Shared:
!                    @GOCART/ESMF/Shared/Chem_AeroGeneric.F90 

     implicit none

     type (ESMF_State),                          intent(inout)     :: state
     character (len=*),                          intent(in   )     :: label
     character (len=*),                          intent(in   )     :: label2
     type (ESMF_Grid),                           intent(inout)     :: grid
     integer,                                    intent(in   )     :: typekind
     real, pointer, dimension(:,:,:), optional,  intent(in   )     :: ptr
     integer,                                    intent(  out)     :: rc

     ! locals
     type (ESMF_Field)                                             :: field
     character (len=ESMF_MAXSTR)                                   :: field_name

     __Iam__('add_aero')

!----------------------------------------------------------------------------------
!    Begin...

     call ESMF_AttributeSet (state, name=trim(label), value=trim(label2),  __RC__)

     call ESMF_AttributeGet (state, name=trim(label), value=field_name, __RC__)
     if (field_name /= '') then
        field = MAPL_FieldCreateEmpty(trim(field_name), grid, __RC__)
        if (trim(field_name) == 'PLE') then
           call MAPL_FieldAllocCommit (field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationEdge, typekind=typekind, hw=0, __RC__)
        else if ((trim(field_name) == 'FRLAND') .or. (trim(field_name) == 'monochromatic_EXT')) then
            call MAPL_FieldAllocCommit(field, dims=MAPL_DimsHorzOnly, location=MAPL_VLocationCenter, typekind=typekind, hw=0, __RC__)
        else
           call MAPL_FieldAllocCommit (field, dims=MAPL_DimsHorzVert, location=MAPL_VLocationCenter, typekind=typekind, hw=0, __RC__)
        end if
        call MAPL_StateAdd (state, field, __RC__)
     end if

!   if (field_name /= '') then
!       field = ptr
!       call MAPL_StateAdd (state, field, __RC__)
!   end if

     RETURN_(ESMF_SUCCESS)

  end subroutine add_aero
!=====================================================================================


!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Run_ --- Runs CARMA
!
! !INTERFACE:
!

   SUBROUTINE Run_ ( gc, import, export, clock, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: clock      ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: gc      ! Grid Component
   type(ESMF_State), intent(inout) :: import     ! Import State
   type(ESMF_State), intent(inout) :: export     ! Export State
   integer, intent(out) ::  rc                    ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  27Feb2005 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------

!  ErrLog Variables
!  ----------------
   character(len=ESMF_MAXSTR)      :: IAm = 'Run_'
   integer                         :: STATUS
   character(len=ESMF_MAXSTR)      :: COMP_NAME

   type(CARMA_GridComp), pointer   :: gcCARMA     ! Grid Component
   integer                         :: nymd, nhms  ! time
   real                            :: cdt         ! chemistry timestep (secs)
   type(Chem_Array), pointer       :: qa(:)
   integer                         :: i1=1, i2, ig=0, im  ! dist grid indices
   integer                         :: j1=1, j2, jg=0, jm  ! dist grid indices
   integer                         :: km, nq              ! dist grid indices
   integer                         :: k, n, dims(3), l, ijl, iq
   real                            :: qmin, qmax

   type(ESMF_Config)               :: CF
   type(ESMF_Grid)                 :: grid
   type(ESMF_Time)                 :: TIME

   type(ESMF_State)                :: internal
   type(MAPL_VarSpec), pointer     :: InternalSpec(:)
   type(MAPL_MetaComp), pointer    :: ggState      ! GEOS Generic State


   real, pointer, dimension(:,:)   :: LATS
   real, pointer, dimension(:,:)   :: LONS
   type (MAPL_SunOrbit)            :: ORBIT


!  Declare pointers to IMPORT/EXPORT/INTERNAL states 
!  -------------------------------------------------
#  include "CARMAchem_DeclarePointer___.h"

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, __RC__)
   Iam = trim(COMP_NAME) // '::Run_'

!  Get my internal MAPL_Generic state
!  -----------------------------------
   CALL MAPL_GetObjectFromGC(GC, ggState, __RC__)

!  Start a comprehensive timer
!  ---------------------------
   CALL MAPL_TimerOn(ggState, "TOTAL")
   CALL MAPL_TimerOn(ggState, "RUN")

!  Get pointers to IMPORT/EXPORT/INTERNAL states 
!  ---------------------------------------------
#  include "CARMAchem_GetPointer___.h"

!  Get parameters from gc and clock
!  --------------------------------
   call extract_ ( gc, clock, gcCARMA, qa, nymd, nhms, cdt, STATUS )
   VERIFY_(STATUS)

!  Run
!  ---
   call CARMA_Emissions  ( gcCARMA, qa, import, export, nymd, nhms, &
                           cdt, STATUS )
   VERIFY_(STATUS)

   call CARMA_GridCompRun ( gcCARMA, qa, import, export, nymd, nhms, &
                            cdt, STATUS )
   VERIFY_(STATUS)

   call CARMA_DryDeposition  ( gcCARMA, qa, import, export, nymd, nhms, &
                               cdt, STATUS )
   VERIFY_(STATUS)

   call CARMA_WetRemoval  ( gcCARMA, qa, import, export, nymd, nhms, &
                            cdt, STATUS )
   VERIFY_(STATUS)

   call CARMA_Convection  ( gcCARMA, qa, import, export, nymd, nhms, &
                            cdt, STATUS )
   VERIFY_(STATUS)

   call CARMA_ComputeDiags  ( gcCARMA, qa, import, export, nymd, nhms, &
                              cdt, STATUS )
   VERIFY_(STATUS)

   CALL MAPL_TimerOff(ggState, "RUN")
   CALL MAPL_TimerOff(ggState, "TOTAL")

   RETURN_(ESMF_SUCCESS)

   END SUBROUTINE Run_

!-------------------------------------------------------------------------
!     NASA/GSFC, Global Modeling and Assimilation Office, Code 610.1     !
!-------------------------------------------------------------------------
!BOP
!
! !IROUTINE:  Finalize_ --- Finalize CARMA
!
! !INTERFACE:
!

   SUBROUTINE Finalize_ ( gc, import, export, clock, rc )

! !USES:

  implicit NONE

! !INPUT PARAMETERS:

   type(ESMF_Clock),  intent(inout) :: clock      ! The clock

! !OUTPUT PARAMETERS:

   type(ESMF_GridComp), intent(inout)  :: gc      ! Grid Component
   type(Chem_Array), pointer           :: qa(:)
   type(ESMF_State), intent(inout) :: import     ! Import State
   type(ESMF_State), intent(inout) :: export     ! Export State
   integer, intent(out) ::  rc                    ! Error return code:
                                                  !  0 - all is well
                                                  !  1 - 

! !DESCRIPTION: This is a simple ESMF wrapper.
!
! !REVISION HISTORY:
!
!  27Feb2005 da Silva  First crack.
!
!EOP
!-------------------------------------------------------------------------


!  ErrLog Variables
!  ----------------
   character(len=ESMF_MAXSTR)      :: IAm = 'Finalize_'
   integer                         :: STATUS
   character(len=ESMF_MAXSTR)      :: COMP_NAME

   type(CARMA_GridComp), pointer   :: gcCARMA     ! Grid Component
   integer                         :: nymd, nhms  ! time
   real                            :: cdt         ! chemistry timestep (secs)

   type(CARMAchem_state), pointer     :: state

   type(MAPL_MetaComp), pointer    :: ggState      ! GEOS Generic State

!  Get my name and set-up traceback handle
!  ---------------------------------------
   call ESMF_GridCompGet( GC, NAME=COMP_NAME, RC=STATUS )
   VERIFY_(STATUS)
   Iam = trim(COMP_NAME) // 'Finalize_'

!  Get my internal MAPL_Generic state
!  -----------------------------------
   CALL MAPL_GetObjectFromGC(GC, ggState, RC=STATUS)

!  Start timers
!  ------------
   CALL MAPL_TimerOn(ggState, "TOTAL")
   CALL MAPL_TimerOn(ggState, "FINALIZE")

!  Get ESMF parameters from gc and clock
!  -------------------------------------
   call extract_ ( gc, clock, gcCARMA, qa, nymd, nhms, cdt, STATUS, &
                   state = state )
   VERIFY_(STATUS)

!  Call ESMF version
!  -----------------
   call CARMA_GridCompFinalize ( gcCARMA, import, export, &
                                 nymd, nhms, cdt, STATUS )
   VERIFY_(STATUS)

!  Destroy Mie Tables
!  ------------------
   call CARMA_DestroyMieTables(gcCARMA, rc)

!  Destroy emissions
!  -----------------
   if(associated(gcCARMA%vLat))              deallocate( gcCARMA%vLat, __STAT__)
   if(associated(gcCARMA%vLon))              deallocate( gcCARMA%vLon, __STAT__)
   if(associated(gcCARMA%vSO2))              deallocate( gcCARMA%vSO2, __STAT__)
   if(associated(gcCARMA%vElev))             deallocate( gcCARMA%vElev, __STAT__)
   if(associated(gcCARMA%vCloud))            deallocate( gcCARMA%vCloud, __STAT__)

!  Destroy Legacy state
!  --------------------
   call registry_destroy_ (state%CARMAreg)
   deallocate ( state%CARMAreg, state%qa, state%gcCARMA, state%chemReg, __STAT__)
   VERIFY_(STATUS)

!  Stop timers
!  -----------
   CALL MAPL_TimerOff(ggState, "FINALIZE")
   CALL MAPL_TimerOff(ggState, "TOTAL")

!  Finalize MAPL Generic.  Atanas says, "Do not deallocate foreign objects."
!  -------------------------------------------------------------------------
   call MAPL_GenericFinalize ( gc, import, export, clock,  RC=STATUS )
   VERIFY_(STATUS)

   RETURN_(ESMF_SUCCESS)

   END SUBROUTINE Finalize_

!.......................................................................
    subroutine extract_ ( gc, clock, gcCARMA, qa, nymd, nhms, cdt, &
                          rc, state )

    type(ESMF_GridComp), intent(INout)  :: gc
    type(ESMF_Clock), intent(in)     :: clock
    type(CARMA_GridComp), pointer    :: gcCARMA
    type(Chem_Array), pointer        :: qa(:)
    integer, intent(out)             :: nymd, nhms
    real, intent(out)                :: cdt
    integer, intent(out)             :: rc
    type(MAPL_MetaComp), pointer               :: ggState
    type(CARMAchem_state), pointer, optional   :: state


    type(CARMAchem_state), pointer    :: myState

!   ErrLog Variables
!   ----------------
    character(len=ESMF_MAXSTR)      :: IAm
    integer                         :: STATUS
    character(len=ESMF_MAXSTR)      :: COMP_NAME

    type(ESMF_Alarm)                :: ALARM
    type(ESMF_TimeInterval)         :: RingInterval

    type(ESMF_Time)      :: TIME
    type(ESMF_Config)    :: CF
    type(CARMAchem_Wrap) :: wrap
    integer              :: IYR, IMM, IDD, IHR, IMN, ISC
    real(ESMF_KIND_R8)   :: dt_r8


!   Get my name and set-up traceback handle
!   ---------------------------------------
    call ESMF_GridCompGet( GC, NAME=COMP_NAME, __RC__ )
    Iam = trim(COMP_NAME) // '::' // 'extract_'

    rc = 0

!   Get my internal MAPL_Generic state
!   -----------------------------------
    call MAPL_GetObjectFromGC ( GC, ggState, __RC__ )


!   Get my internal state
!   ---------------------
    call ESMF_UserCompGetInternalState(gc, 'CARMA_state', WRAP, STATUS)
    VERIFY_(STATUS)
    myState => wrap%ptr
    if ( present(state) ) then
         state => wrap%ptr
    end if

    if ( .not. associated(myState%gcCARMA) ) then
         allocate ( myState%gcCARMA, stat=STATUS )
         VERIFY_(STATUS)
    end if

    if ( .not. associated(myState%CARMAreg) ) then
         allocate ( myState%CARMAreg, stat=STATUS )
         VERIFY_(STATUS)
    end if

    if ( .not. associated(myState%qa) ) then
         allocate ( myState%qa(myState%CARMAreg%nq), stat=STATUS )
         VERIFY_(STATUS)
    end if

    gcCARMA => myState%gcCARMA
    gcCARMA%CARMAreg => myState%CARMAreg
    qa => myState%qa

!   Get the configuration
!   ---------------------
    call ESMF_GridCompGet ( GC, CONFIG = CF, __RC__ )

!   Get time step
!   -------------
    call MAPL_Get(ggState, RUNALARM=ALARM, __RC__ )
    call ESMF_AlarmGet(ALARM, ringInterval=RingInterval, __RC__)

    call ESMF_TimeIntervalGet(RingInterval, s_r8=dt_r8, __RC__)
    cdt = real(dt_r8)

    call ESMF_ClockGet(CLOCK,currTIME=TIME,rc=STATUS)
    VERIFY_(STATUS)

!   Need code to extract nymd(20050205), nhms(120000) from clock
!   ------------------------------------------

    call ESMF_ClockGet(CLOCK,currTIME=TIME, __RC__ )
    call ESMF_TimeGet(TIME ,YY=IYR, MM=IMM, DD=IDD, H=IHR, M=IMN, S=ISC, __RC__ )
    call MAPL_PackTime(NYMD,IYR,IMM,IDD)
    call MAPL_PackTime(NHMS,IHR,IMN,ISC)

    RETURN_(ESMF_SUCCESS)

   end subroutine extract_

!.......................................................................
   subroutine registry_ (r)

   type(CARMA_Registry), pointer   :: r
   CHARACTER(LEN=255) :: string
   integer :: ios, ier(20), i, j, n, rc
   character(len=1)   :: elem_str

!  Load resource file
!  ------------------
   CALL I90_loadf ( TRIM(r%rcfilen), ier(1) )
   IF ( ier(1) .NE. 0 ) THEN
    CALL final_(10)
    RETURN
   END IF
   ier(1)=0

!  Particle/Gas/Radiation structure
!  --------------------------------
!  NBIN
   call i90_label ( 'NBIN:', ier(1) )
   r%NBIN = i90_gint ( ier(2) )
   if ( any(ier(1:2) /= 0) ) then
      call final_(20)
      return
   end if

!  NGROUP
   call i90_label ( 'NGROUP:', ier(1) )
   r%NGROUP = i90_gint ( ier(2) )
   if ( any(ier(1:2) /= 0) ) then
      call final_(20)
      return
   end if

!  NELEM
   call i90_label ( 'NELEM:', ier(1) )
   r%NELEM = i90_gint ( ier(2) )
   if ( any(ier(1:2) /= 0) ) then
      call final_(20)
      return
   end if

!  NGAS
   call i90_label ( 'NGAS:', ier(1) )
   r%NGAS = i90_gint ( ier(2) )
   if ( any(ier(1:2) /= 0) ) then
      call final_(20)
      return
   end if

!  NSOLUTE
   call i90_label ( 'NSOLUTE:', ier(1) )
   r%NSOLUTE = i90_gint ( ier(2) )
   if ( any(ier(1:2) /= 0) ) then
      call final_(20)
      return
   end if

!  NWAVE
   call i90_label ( 'NWAVE:', ier(1) )
   r%NWAVE = i90_gint ( ier(2) )
   if ( any(ier(1:2) /= 0) ) then
      call final_(20)
      return
   end if

!  Group Characteristics
!  ---------------------
   allocate ( r%rmrat(r%NGROUP),  &
              r%rmin(r%NGROUP),   &
              r%ishape(r%NGROUP), &
              r%eshape(r%NGROUP), &
              r%fscav(r%NGROUP),  &
              r%irhswell(r%NGROUP), &
              r%irhswcomp(r%NGROUP), &
              r%groupname(r%NGROUP), stat=ios )
   if ( ios /= 0) then
    call final_(100)
    return
   endif

   call i90_label ( 'RMRAT:', ier(1) )
   do j = 1, r%NGROUP
    r%rmrat(j) = i90_gfloat(ier(j+1))
   end do
   if(any(ier(1:r%NGROUP+1) /= 0)) then
    call final_(101)
    return
   endif
    
   call i90_label ( 'RMIN:', ier(1) )
   do j = 1, r%NGROUP
    r%rmin(j) = i90_gfloat(ier(j+1))
   end do
   if(any(ier(1:r%NGROUP+1) /= 0)) then
    call final_(101)
    return
   endif
    
   call i90_label ( 'ISHAPE:', ier(1) )
   do j = 1, r%NGROUP
    r%ishape(j) = i90_gint(ier(j+1))
   end do
   if(any(ier(1:r%NGROUP+1) /= 0)) then
    call final_(101)
    return
   endif
    
   call i90_label ( 'ESHAPE:', ier(1) )
   do j = 1, r%NGROUP
    r%eshape(j) = i90_gfloat(ier(j+1))
   end do
   if(any(ier(1:r%NGROUP+1) /= 0)) then
    call final_(101)
    return
   endif
    
   call i90_label ( 'FSCAV:', ier(1) )
   do j = 1, r%NGROUP
    r%fscav(j) = i90_gfloat(ier(j+1))
   end do
   if(any(ier(1:r%NGROUP+1) /= 0)) then
    call final_(101)
    return
   endif
    
   call i90_label ( 'IRHSWELL:', ier(1) )
   do j = 1, r%NGROUP
    r%irhswell(j) = i90_gint(ier(j+1))
   end do
   if(any(ier(1:r%NGROUP+1) /= 0)) then
    call final_(101)
    return
   endif
    
   call i90_label ( 'IRHSWCOMP:', ier(1) )
   do j = 1, r%NGROUP
    r%irhswcomp(j) = i90_gint(ier(j+1))
   end do
   if(any(ier(1:r%NGROUP+1) /= 0)) then
    call final_(101)
    return
   endif
    
   call i90_label ( 'GROUPNAME:', ier(1) )
   do j = 1, r%NGROUP
    call i90_gtoken(string, ier(j+1))
    r%groupname(j) = trim(string)
   end do
   if(any(ier(1:r%NGROUP+1) /= 0)) then
    call final_(101)
    return
   endif
    

!  Element Characteristics
!  -----------------------
   allocate ( r%rhop(r%NELEM), r%igroup(r%NELEM), r%itype(r%NELEM), &
              r%elemname(r%NELEM), r%icomposition(r%NELEM), &
              r%distribution(r%NELEM, r%NBIN), stat=ios )
   if ( ios /= 0) then
    call final_(100)
    return
   endif

   call i90_label ( 'IGROUP:', ier(1) )
   do i = 1, r%NELEM
    r%igroup(i) = i90_gint(ier(i+1))
   end do
   if(any(ier(1:r%NELEM+1) /= 0)) then
    call final_(101)
    return
   endif
    
   call i90_label ( 'ICOMPOSITION:', ier(1) )
   do i = 1, r%NELEM
    r%icomposition(i) = i90_gint(ier(i+1))
   end do
   if(any(ier(1:r%NELEM+1) /= 0)) then
    call final_(101)
    return
   endif
    
   call i90_label ( 'ITYPE:', ier(1) )
   do i = 1, r%NELEM
    r%itype(i) = i90_gint(ier(i+1))
   end do
   if(any(ier(1:r%NELEM+1) /= 0)) then
    call final_(101)
    return
   endif
    
   call i90_label ( 'RHOP:', ier(1) )
   do i = 1, r%NELEM
    r%rhop(i) = i90_gfloat(ier(i+1))
   end do
   if(any(ier(1:r%NELEM+1) /= 0)) then
    call final_(101)
    return
   endif
    
   call i90_label ( 'ELEMNAME:', ier(1) )
   do i = 1, r%NELEM
    call i90_gtoken(string, ier(i+1))
    r%elemname(i) = trim(string)
   end do
   if(any(ier(1:r%NELEM+1) /= 0)) then
    call final_(101)
    return
   endif

   do i = 1, r%NELEM
    write(elem_str, '(i0)') i
    call i90_label ( 'DISTRIBUTION'//elem_str//':', ier(1) )
    do j = 1, r%NBIN
     r%distribution(i,j) = i90_gfloat(ier(j+1))
    end do
    if(any(ier(1:r%NBIN+1) /= 0)) then
     call final_(101)
     return
    endif
   end do

!  Gas Characteristics
!  -------------------
   if(r%NGAS .gt. 0) then
    allocate ( r%gasname(r%NGAS), r%igcomp(r%NGAS), &
               r%igvapreq(r%NGAS), stat = ios )
    if ( ios /= 0) then
     call final_(200)
     return
    endif

    call i90_label ( 'GASNAME:', ier(1) )
    do i = 1, r%NGAS
     call i90_gtoken(string, ier(i+1))
     r%gasname(i) = trim(string)
    end do
    if(any(ier(1:r%NGAS+1) /= 0)) then
     call final_(201)
     return
    endif

    call i90_label ( 'IGCOMP:', ier(1) )
    do i = 1, r%NGAS
     r%igcomp(i) = i90_gint(ier(i+1))
    end do
    if(any(ier(1:r%NGAS+1) /= 0)) then
     call final_(202)
     return
    endif

    call i90_label ( 'IGVAPREQ:', ier(1) )
    do i = 1, r%NGAS
     r%igvapreq(i) = i90_gint(ier(i+1))
    end do
    if(any(ier(1:r%NGAS+1) /= 0)) then
     call final_(203)
     return
    endif

  endif



!  Microphysical process flags and timesteps
!  -----------------------------------------
!  Note that the labels are not required to be present
   call i90_label ( 'DO_COAG:', ios )
   if(ios .eq. 0) r%do_coag = i90_gint( ier(1) )
   call i90_label ( 'DO_GROW:', ios )
   if(ios .eq. 0) r%do_grow = i90_gint( ier(2) )
   call i90_label ( 'DO_SUBSTEP:', ios )
   if(ios .eq. 0) r%do_substep = i90_gint( ier(3) )
   call i90_label ( 'DO_THERMO:', ios )
   if(ios .eq. 0) r%do_thermo = i90_gint( ier(4) )
   call i90_label ( 'DO_VDIFF:', ios )
   if(ios .eq. 0) r%do_vdiff = i90_gint( ier(5) )
   call i90_label ( 'DO_VTRAN:', ios )
   if(ios .eq. 0) r%do_vtran = i90_gint( ier(6) )
   call i90_label ( 'DO_FIXEDINIT:', ios )
   if(ios .eq. 0) r%do_fixedinit = i90_gint( ier(6) )
   call i90_label ( 'VF_CONST:', ios )
   if(ios .eq. 0) r%vf_const = i90_gfloat( ier(7) )
   call i90_label ( 'MINSUBSTEPS:', ios )
   if(ios .eq. 0) r%minsubsteps = i90_gint( ier(8) )
   call i90_label ( 'MAXSUBSTEPS:', ios )
   if(ios .eq. 0) r%maxsubsteps = i90_gint( ier(9) )
   call i90_label ( 'MAXRETRIES:', ios )
   if(ios .eq. 0) r%maxretries = i90_gint( ier(10) )
   call i90_label ( 'CONMAX:', ios )
   if(ios .eq. 0) r%conmax = i90_gfloat( ier(11) )
   if(any(ier(1:11) /= 0)) then
    call final_(102)
    return
   endif

!  Species specific
!  ----------------
!  Dust
   call i90_label ( 'dust_emissions_fudgefactor:', ier(1) )
   r%dust_emissions_fudgefactor = i90_gfloat ( ier(2))
   if ( any(ier(1:2) /= 0 )) then
         call final_(40)
         return
   end if
   allocate(r%dmass_dust(r%NBIN), stat = ios)

!  Sea Salt
   call i90_label ( 'seasalt_emissions_fudgefactor:', ier(1) )
   r%seasalt_emissions_fudgefactor = i90_gfloat ( ier(2))
   if ( any(ier(1:2) /= 0 )) then
         call final_(40)
         return
   end if

!  Black Carbon

!  Organic Carbon 
   call i90_label ( 'organic_matter_to_organic_carbon_ratio:', ier(1) )
   r%organic_matter_to_organic_carbon_ratio = i90_gfloat ( ier(2))
   if ( any(ier(1:2) /= 0 )) then
         call final_(40)
         return
   end if
   call i90_label ( 'fraction_terpene_to_organic_carbon:', ier(1) )
   r%fraction_terpene_to_organic_carbon = i90_gfloat ( ier(2))
   if ( any(ier(1:2) /= 0 )) then
         call final_(40)
         return
   end if

!  Get any requested point emissions
!  ---------------------------------
!  Sulfate
   ier(:) = 0
   call i90_label  ( 'point_emissions_srcfilen_sulfate:', ier(1) )
   call i90_gtoken ( r%point_emissions_srcfilen_sulfate,  ier(2) )
   if ( ier(1) /= 0 ) then
        r%doing_point_emissions_sulfate = .FALSE. ! if rc is missing, don't fuss
   else if ( any(ier(2:2) /= 0) ) then
         call final_(42) ! this means point emissions info is messed up, abort
         return
   else
         if ( (index(r%point_emissions_srcfilen_sulfate,'/dev/null')>0) ) then
                     r%doing_point_emissions_sulfate = .FALSE. ! disable it if no file specified
         else
                     r%doing_point_emissions_sulfate = .TRUE.  ! we are good to go
         end if
   end if

!  Ash
   ier(:) = 0
   call i90_label  ( 'point_emissions_srcfilen_ash:', ier(1) )
   call i90_gtoken ( r%point_emissions_srcfilen_ash,  ier(2) )
   if ( ier(1) /= 0 ) then
        r%doing_point_emissions_ash = .FALSE. ! if rc is missing, don't fuss
   else if ( any(ier(2:2) /= 0) ) then
         call final_(42) ! this means point emissions info is messed up, abort
         return
   else
         if ( (index(r%point_emissions_srcfilen_ash,'/dev/null')>0) ) then
                     r%doing_point_emissions_ash = .FALSE. ! disable it if no file specified
         else
                     r%doing_point_emissions_ash = .TRUE.  ! we are good to go
         end if
   end if
!  Dust
   ier(:) = 0
   call i90_label  ( 'point_emissions_srcfilen_dust:', ier(1) )
   call i90_gtoken ( r%point_emissions_srcfilen_dust,  ier(2) )
   if ( ier(1) /= 0 ) then
        r%doing_point_emissions_dust = .FALSE. ! if rc is missing, don't fuss
   else if ( any(ier(2:2) /= 0) ) then
         call final_(42) ! this means point emissions info is messed up, abort
         return
   else
         if ( (index(r%point_emissions_srcfilen_dust,'/dev/null')>0) ) then
                     r%doing_point_emissions_dust = .FALSE. ! disable it if no file specified
         else
                     r%doing_point_emissions_dust = .TRUE.  ! we are good to go
         end if
   end if

!  Mie Tables
!  ----------
!   Set the number of channels to calculate over
!   --------------------------------------------
    call i90_label ( 'n_channels:', ios )
    if ( ios /= 0 ) then
     call final_(60)
    else
     r%nchannels = i90_gint ( ios )
     if ( ios /= 0 ) call final_(61)
    end if

!   Set the number of moments
!   -------------------------
    call i90_label ( 'n_moments:', ios )
    if ( ios /= 0 ) then
     r%nmoments = 0
    else
     r%nmoments = i90_gint ( ios )
     if ( ios /= 0 ) call final_(62)
    end if

!   Set the channels to calculate over
!   ----------------------------------
    allocate( r%channels(r%nchannels), stat = ios )
    call i90_label ( 'r_channels:', ios )
    if ( ios /= 0 ) then
     call final_(63)
    else
     do n = 1, r%nchannels
      r%channels(n) = i90_gfloat ( ios )
      if ( ios /= 0 ) call final_(64)
     enddo
    end if

    call i90_label ( 'filename_optical_properties_DU:', ios )
     if ( ios /= 0 ) then
      call final_(65)
     else
      call i90_gtoken ( r%du_optics_file, ios )
      if ( ios /= 0 ) call final_(66)
     end if

    call i90_label ( 'filename_optical_properties_SS:', ios )
     if ( ios /= 0 ) then
      call final_(67)
     else
      call i90_gtoken ( r%ss_optics_file, ios )
      if ( ios /= 0 ) call final_(68)
     end if

    call i90_label ( 'filename_optical_properties_BC:', ios )
     if ( ios /= 0 ) then
      call final_(67)
     else
      call i90_gtoken ( r%bc_optics_file, ios )
      if ( ios /= 0 ) call final_(68)
     end if

    call i90_label ( 'filename_optical_properties_OC:', ios )
     if ( ios /= 0 ) then
      call final_(69)
     else
      call i90_gtoken ( r%oc_optics_file, ios )
      if ( ios /= 0 ) call final_(70)
     end if

    call i90_label ( 'filename_optical_properties_SU:', ios )
     if ( ios /= 0 ) then
      call final_(69)
     else
      call i90_gtoken ( r%su_optics_file, ios )
      if ( ios /= 0 ) call final_(70)
     end if


   r%nq = r%NBIN*r%NELEM + r%NGAS + r%NGAS + r%NGAS + r%NGAS + 1

   if(r%nq .gt. 0) r%doing_CARMA = .true.

   if(r%doing_CARMA) allocate(r%vname(r%nq), stat=ios)
   if(ios /= 0) call final_(100)

   CALL I90_release()

   end subroutine registry_





   subroutine registry_destroy_ (r)

   type(CARMA_Registry), pointer   :: r
   integer :: ios

!  Group Characteristics
!  ---------------------
   deallocate ( r%rmrat, &
                r%rmin, &
                r%ishape, &
                r%eshape, &
                r%fscav, &
                r%irhswell, &
                r%irhswcomp, &
                r%groupname, stat=ios )
   if ( ios /= 0) then
    call final_(100)
    return
   endif

!  Element Characteristics
!  -----------------------
   deallocate ( r%rhop, r%igroup, r%itype, &
                r%elemname, r%icomposition, &
                r%distribution, stat=ios )
   if ( ios /= 0) then
    call final_(100)
    return
   endif

!  Gas Characteristics
!  -------------------
   if(r%NGAS .gt. 0) then
    deallocate ( r%gasname, r%igcomp, &
                 r%igvapreq, stat = ios )
    if ( ios /= 0) then
     call final_(200)
     return
    endif
   endif

!   Mie Tables
!   ----------
    deallocate( r%channels, stat = ios )
    if ( ios /= 0) then
     call final_(300)
     return
    endif

    if(r%doing_CARMA) deallocate(r%vname, stat=ios)
    if(ios /= 0) then
     call final_(400)
     return
    endif

!   Other
!   -----
    deallocate( r%dmass_dust, stat=ios)

   end subroutine registry_destroy_





!  print the CARMA registry (for checking)
   subroutine registry_print_ (r)

   type(CARMA_Registry), pointer   :: r
   CHARACTER(LEN=255) :: string
   integer :: ios, ier(20), i, j, rc

   if(MAPL_AM_I_ROOT()) then
    print *, 'CARMAchem_GridCompMod: registry_print_'
    print *, 'NBIN:    ', r%NBIN
    print *, 'NGROUP:  ', r%NGROUP
    print *, 'NELEM:   ', r%NELEM
    print *, 'NGAS:    ', r%NGAS
    print *, 'NSOLUTE: ', r%NSOLUTE
    print *, 'NWAVE:   ', r%NWAVE

!   Requested point emissions
!   -------------------------
    print *, 'Requested point emissions (sulfate) : ', r%doing_point_emissions_sulfate
    if(r%doing_point_emissions_sulfate) then
     print *, ' Point emissions template (sulfate): ', r%point_emissions_srcfilen_sulfate
    endif
    print *, 'Requested point emissions (dust) : ', r%doing_point_emissions_dust
    if(r%doing_point_emissions_dust) then
     print *, ' Point emissions template (dust): ', r%point_emissions_srcfilen_dust
    endif
    print *, 'Requested point emissions (ash) : ', r%doing_point_emissions_ash
    if(r%doing_point_emissions_ash) then
     print *, ' Point emissions template (ash): ', r%point_emissions_srcfilen_ash
    endif

!   Group Characteristics
!   -----------------------
    if(r%NGROUP .gt. 0) then
     do j = 1, r%NGROUP
      print *, 'GROUP (',j,'): GROUPNAME = ',trim(r%groupname(j)), &
               ', RMIN = ',r%rmin(j),', RMRAT = ',r%rmrat(j), &
               ', ISHAPE = ',r%ishape(j),', ESHAPE = ',r%eshape(j), &
               ', FSCAV = ',r%fscav(j),', IRHSWELL = ', r%irhswell(j), &
               ', IRHSWCOMP = ',r%irhswcomp(j)
     enddo
    endif 

!   Element Characteristics
!   -----------------------
    if(r%NELEM .gt. 0) then
     do j = 1, r%NELEM
      print *, 'ELEMENT (',j,'): ELEMNAME = ',trim(r%elemname(j)), &
               ', IGROUP = ',r%igroup(j),', RHOP = ',r%rhop(j)
     enddo
    endif 

!   Gas Characteristics
!   -------------------
    if(r%NGAS .gt. 0) then
     do j = 1, r%NGAS
      print *, 'GAS (',j,'): GASNAME = ',trim(r%gasname(j)), &
               ', IGCOMP = ',r%igcomp(j),', IGVAPREQ = ',r%igvapreq(j)
     enddo
    endif 

!   Microphysical process flags and timesteps
!   -----------------------------------------
    print *, 'DO_COAG:       ', r%do_coag
    print *, 'DO_GROW:       ', r%do_grow
    print *, 'DO_SUBSTEP:    ', r%do_substep
    print *, 'DO_THERMO:     ', r%do_thermo
    print *, 'DO_VDIFF:      ', r%do_vdiff
    print *, 'DO_VTRAN:      ', r%do_vtran
    print *, 'VF_CONST:      ', r%vf_const
    print *, 'DO_FIXEDINIT:  ', r%do_fixedinit
    print *, 'MINSUBSTEPS:   ', r%minsubsteps
    print *, 'MAXSUBSTEPS:   ', r%maxsubsteps
    print *, 'MAXRETRIES:    ', r%maxretries
    print *, 'CONMAX:        ', r%conmax

   endif

   end subroutine registry_print_





   SUBROUTINE final_(ierr)
   INTEGER :: ios, ierr, rc
   CALL I90_release()
   rc = ierr
   END SUBROUTINE final_


subroutine run_aerosol_optics(state, rc)

  implicit none

! Arguments
! ---------
  type(ESMF_State)     :: state
  integer, intent(out) :: rc


! Local
! ---------
  integer                                 :: n_aerosols
  character(len=ESMF_MAXSTR), allocatable :: aerosol_names(:)
  type(ESMF_FieldBundle)                  :: aerosols

  real, dimension(:,:,:), pointer         :: ple
  real, dimension(:,:,:), pointer         :: rh
  real, dimension(:,:,:), pointer         :: var
  real, dimension(:,:,:), pointer         :: q
  real, dimension(:,:,:,:), pointer       :: q_4d

  real, dimension(:,:,:), allocatable     :: dp, f_p

  character(len=ESMF_MAXSTR)              :: fld_name
  type(ESMF_Field)                        :: fld

  real, dimension(:,:,:,:), allocatable   :: ext, ssa, asy  ! (lon:,lat:,lev:,band:)

  integer                                 :: n
  integer                                 :: i1, j1, i2, j2, km

  integer                                 :: band, offset

  integer                                 :: instance

  integer                                 :: STATUS
  character(len=ESMF_MAXSTR)              :: Iam

  integer, parameter                      :: n_bands = 1

  real    :: x
  integer :: i, j, k

  Iam = 'CARMA::run_aerosol_optics()'


! Mie Table instance/index
! ------------------------
  call ESMF_AttributeGet(state, name='mie_table_instance', value=instance, __RC__)

! Radiation band
! --------------
  band = 0
  call ESMF_AttributeGet(state, name='band_for_aerosol_optics', value=band, __RC__)
  offset = band - n_bands

! Pressure at layer edges 
! ------------------------
  call ESMF_AttributeGet(state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)
  call MAPL_GetPointer(state, ple, trim(fld_name), __RC__)

  i1 = lbound(ple, 1); i2 = ubound(ple, 1)
  j1 = lbound(ple, 2); j2 = ubound(ple, 2)
                       km = ubound(ple, 3)

! Relative humidity
! -----------------
  call ESMF_AttributeGet(state, name='relative_humidity_for_aerosol_optics', value=fld_name, __RC__)
  call MAPL_GetPointer(state, rh, trim(fld_name), __RC__)

  i1 = lbound(rh, 1); i2 = ubound(rh, 1)
  j1 = lbound(rh, 2); j2 = ubound(rh, 2)
                      km = ubound(rh, 3)
  
  call ESMF_StateGet(state, 'AEROSOLS', aerosols, __RC__)
  call ESMF_FieldBundleGet(aerosols, fieldCount=n_aerosols, __RC__)

  allocate(aerosol_names(n_aerosols), __STAT__)
 
  call ESMF_FieldBundleGet(aerosols, FieldNameList=aerosol_names, __RC__)
 
  allocate(ext(i1:i2,j1:j2,km,n_bands), &
           ssa(i1:i2,j1:j2,km,n_bands), &
           asy(i1:i2,j1:j2,km,n_bands), __STAT__)

  allocate(q_4d(i1:i2,j1:j2,km,n_aerosols), __STAT__)

#if (0)
  allocate(dp(i1:i2,j1:j2,km), f_p(i1:i2,j1:j2,km), __STAT__)

  dp  = ple(:,:,1:km) - ple(:,:,0:km-1)
  f_p = dp / MAPL_GRAV

  do n = 1, n_aerosols
      call ESMF_FieldBundleGet(aerosols, trim(aerosol_names(n)), field=fld, __RC__)
      call ESMF_FieldGet(fld, farrayPtr=q, __RC__)

      q_4d(:,:,:,n) = f_p * q
  end do

  call ESMF_AttributeGet(state, name='mie_table_instance', value=instance, __RC__)
  call mie_(carmaMieTable(instance),aerosol_names, n_bands, offset, q_4d, rh, ext, ssa, asy, __RC__)

  deallocate(dp, f_p, __STAT__)
#else
  do n = 1, n_aerosols
      call ESMF_FieldBundleGet(aerosols, trim(aerosol_names(n)), field=fld, __RC__)
      call ESMF_FieldGet(fld, farrayPtr=q, __RC__)

      do k = 1, km
          do j = j1, j2
              do i = i1, i2
                  x = ((PLE(i,j,k) - PLE(i,j,k-1))*0.01)*(100./MAPL_GRAV)
                  q_4d(i,j,k,n) = x * q(i,j,k)
              end do
          end do
      end do
  end do

  call mie_(carmaMieTable(instance), aerosol_names, n_bands, offset, q_4d, rh, ext, ssa, asy, __RC__)
#endif
  
  call ESMF_AttributeGet(state, name='extinction_in_air_due_to_ambient_aerosol', value=fld_name, __RC__)
  if (fld_name /= '') then 
      call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
      var = ext(:,:,:,1)
  end if

  call ESMF_AttributeGet(state, name='single_scattering_albedo_of_ambient_aerosol', value=fld_name, __RC__)
  if (fld_name /= '') then 
      call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
      var = ssa(:,:,:,1)
  end if

  call ESMF_AttributeGet(state, name='asymmetry_parameter_of_ambient_aerosol', value=fld_name, __RC__)
  if (fld_name /= '') then 
      call MAPL_GetPointer(state, var, trim(fld_name), __RC__)
      var = asy(:,:,:,1)
  end if

  deallocate(aerosol_names, ext, ssa, asy, q_4d, __STAT__)

  RETURN_(ESMF_SUCCESS)

contains 

    subroutine mie_(mie_table, aerosol, nb, offset, q, rh, ext, ssa, asy, rc)
     
     implicit none

     type(Chem_Mie),    intent(inout):: mie_table    ! mie table
     character(len=*),  intent(in )  :: aerosol(:)   ! list of aerosols
     integer,           intent(in )  :: nb           ! number of bands
     integer,           intent(in )  :: offset       ! bands offset 
     real,              intent(in )  :: q(:,:,:,:)   ! aerosol mass mixing ratio, kg kg-1
     real,              intent(in )  :: rh(:,:,:)    ! relative humidity

     real,              intent(out)  :: ext(:,:,:,:) ! extinction
     real,              intent(out)  :: ssa(:,:,:,:) ! SSA
     real,              intent(out)  :: asy(:,:,:,:) ! asymmetry parameter

     integer,           intent(out)  :: rc

     ! local
     integer :: STATUS
     character(len=ESMF_MAXSTR) :: Iam='aerosol_optics::mie_' 

     integer :: l, idx, na

     real(kind=8) :: ext_(size(ext,1),size(ext,2),size(ext,3),size(ext,4))
     real(kind=8) :: ssa_(size(ext,1),size(ext,2),size(ext,3),size(ext,4))
     real(kind=8) :: asy_(size(ext,1),size(ext,2),size(ext,3),size(ext,4))

     na = size(aerosol)

     _ASSERT(na == size(q,4), 'needs informative message')

     ext_ = 0.0d0
     ssa_ = 0.0d0
     asy_ = 0.0d0

     do l = 1, na
        idx = Chem_MieQueryIdx(mie_table, 'CARMA::'//trim(aerosol(l)), __RC__)

        call Chem_MieQueryAllBand4D(mie_table, idx, nb, offset, q(:,:,:,l), rh, ext, ssa, asy, __RC__)

        ext_ = ext_ +          ext     ! total extinction
        ssa_ = ssa_ +     (ssa*ext)    ! total scattering
        asy_ = asy_ + asy*(ssa*ext)    ! sum of (asy * sca)
     end do

     ext = ext_
     ssa = ssa_
     asy = asy_

     RETURN_(ESMF_SUCCESS)

    end subroutine mie_

 end subroutine run_aerosol_optics

!=====================================================================================================

  subroutine aerosol_activation_properties(state, rc)

    implicit none

!   Arguments
!   ---------
    type(ESMF_State)     :: state
    integer, intent(out) :: rc

!   Local
!   ---------
    character(len=ESMF_MAXSTR)      :: mode              ! mode name
    character(len=ESMF_MAXSTR)      :: mode_             ! lowercase mode name

    type(ESMF_State)                :: child_state

    real(MAPL_R4), dimension(:,:,:), pointer :: ple, temperature, num, diameter, sigma, density, hygroscopicity
    real(MAPL_R4), dimension(:,:,:), pointer     :: f_dust, f_soot, f_organic
    !real(ESMF_KIND_R8), dimension(:,:,:), pointer :: ptr_3d_r8
    real(MAPL_R4), dimension(:,:,:), pointer :: ptr_3d
    type(ESMF_TypeKind_Flag) :: tk
    !real(ESMF_KIND_R8), dimension(:,:,:), allocatable :: q
    real(MAPL_R4), dimension(:,:,:), allocatable :: q

    real                            :: max_clean          ! max mixing ratio before considered polluted
    real                            :: ccn_tuning         ! tunes conversion factors for sulfate
    character(LEN=ESMF_MAXSTR)      :: cld_micro

    character(len=ESMF_MAXSTR)      :: fld_name
    type(ESMF_Field)                :: fld

    integer                         :: i2, j2, k2, n_aerosols
    integer                         :: b, i, j, n, aerosol_bin
    integer                         :: varNameLen, c_idx

!    real(ESMF_KIND_R8), allocatable           :: aeroRmin(:), aeroRmrat(:)
    real(MAPL_R4), allocatable           :: aeroRmin(:), aeroRmrat(:)
    character (len=ESMF_MAXSTR), allocatable  :: aeroList(:)
    type (ESMF_FieldBundle)                   :: aerosols

!   auxilliary parameters
!   ---------------------
    real, parameter :: densSO4 = 1700.0
    real, parameter :: densORG = 1600.0
    real, parameter :: densSS  = 2200.0
    real, parameter :: densDU  = 1700.0
    real, parameter :: densBC  = 1600.0
    real, parameter :: densOC  =  900.0
    real, parameter :: densBR  =  900.0

    real, parameter :: k_SO4   = 0.65
    real, parameter :: k_ORG   = 0.20
    real, parameter :: k_SS    = 1.28
    real, parameter :: k_DU    = 0.0001
    real, parameter :: k_BC    = 0.0001
    real, parameter :: k_OC    = 0.0001
    real, parameter :: k_BR    = 0.0001

    integer, parameter :: UNKNOWN_AEROSOL_MODE = 2015

    __Iam__('CARMA::aerosol_activation_properties')

!   Begin...

!   Get list of child states within state and add to aeroList
!   ---------------------------------------------------------
!    call ESMF_StateGet (state, itemCount=n, __RC__)
!    allocate (itemList(n), __STAT__)
!    allocate (itemTypes(n), __STAT__)
!    call ESMF_StateGet (state, itemNameList=itemList, itemTypeList=itemTypes, __RC__)

!
!    b=0
!    do i = 1, n
!       if ((itemTypes(i) == ESMF_StateItem_State) .and. (trim(itemList(i)(1:2)) /= 'NI')) then
!          b = b + 1
!       end if
!    end do
!
!    allocate (aeroList(b), __STAT__)
!
!    j = 1
!    do i = 1, n
!       if ((itemTypes(i) == ESMF_StateItem_State) .and. (trim(itemList(i)(1:2)) /= 'NI')) then
!          aeroList(j) = trim(itemList(i))
!          j = j + 1
!       end if
!    end do


! The way I did this above...
! Why don't we just get the list of aerosols from the "aerosol_modes" attribute in the AERO bundle?
!        allocate(aero_aci_modes(n_aerosols))
!        do ielem = 1, reg%NELEM
!         igroup = reg%igroup(ielem)
!         if(ielem /= gcCARMA%carma%f_group(igroup)%f_ienconc ) cycle
!         do ibin = 1, reg%NBIN
!          aero_aci_modes((ielem-1)*reg%NBIN + ibin - 1) = reg%vname(n)
!         end do
!        end do
!        call ESMF_AttributeSet(aero, name='number_of_aerosol_modes', value=n_aerosols, __RC__)
!        call ESMF_AttributeSet(aero, name='aerosol_modes', itemcount=n_aerosols, valuelist=aero_aci_modes, __RC__)
!        deallocate(aer_aci_modes)

!   CARMA maintains list of modes
!   -----------------------------
    call ESMF_StateGet(state, 'AEROSOLS', aerosols, __RC__)
    call ESMF_FieldBundleGet(aerosols, fieldCount=n_aerosols, __RC__)
    allocate(aeroList(n_aerosols), __STAT__)
    allocate(aeroRmin(n_aerosols), __STAT__)
    allocate(aeroRmrat(n_aerosols), __STAT__)
    call ESMF_AttributeGet(state, name='aerosol_modes', valuelist=aeroList, __RC__)
    call ESMF_AttributeGet(state, name='aerosol_rmin', valuelist=aeroRmin, __RC__)
    call ESMF_AttributeGet(state, name='aerosol_rmrat', valuelist=aeroRmrat, __RC__)

!   Aerosol mode (set by moist before this function is called)
!   ----------------------------------------------------------
    call ESMF_AttributeGet(state, name='aerosol_mode', value=mode, __RC__)

!   Land fraction
!   -------------
!    call ESMF_AttributeGet(state, name='fraction_of_land_type', value=fld_name, __RC__)
!    call MAPL_GetPointer(state, f_land, trim(fld_name), __RC__)

!   Pressure at layer edges
!   ------------------------
    call ESMF_AttributeGet(state, name='air_pressure_for_aerosol_optics', value=fld_name, __RC__)
    call MAPL_GetPointer(state, ple, trim(fld_name), __RC__)

!   Temperature
!   -----------
    call ESMF_AttributeGet(state, name='air_temperature', value=fld_name, __RC__)
    call MAPL_GetPointer(state, temperature, trim(fld_name), __RC__)

    i2 = ubound(temperature, 1)
    j2 = ubound(temperature, 2)
    k2 = ubound(temperature, 3)

!   Activation activation properties
!   --------------------------------
    call ESMF_AttributeGet(state, name='aerosol_number_concentration', value=fld_name, __RC__)
    call MAPL_GetPointer(state, num, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='aerosol_dry_size', value=fld_name, __RC__)
    call MAPL_GetPointer(state, diameter, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='width_of_aerosol_mode', value=fld_name, __RC__)
    call MAPL_GetPointer(state, sigma, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='aerosol_density', value=fld_name, __RC__)
    call MAPL_GetPointer(state, density, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='aerosol_hygroscopicity', value=fld_name, __RC__)
    call MAPL_GetPointer(state, hygroscopicity, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='fraction_of_dust_aerosol', value=fld_name, __RC__)
    call MAPL_GetPointer(state, f_dust, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='fraction_of_soot_aerosol', value=fld_name, __RC__)
    call MAPL_GetPointer(state, f_soot, trim(fld_name), __RC__)

    call ESMF_AttributeGet(state, name='fraction_of_organic_aerosol', value=fld_name, __RC__)
    call MAPL_GetPointer(state, f_organic, trim(fld_name), __RC__)

!   Sea salt scaling fctor
!   ----------------------
    call ESMF_AttributeGet(state, name='max_q_clean', value=max_clean, __RC__)
    call ESMF_AttributeGet(state, name='cldmicro', value=cld_micro, __RC__)
    call ESMF_AttributeGet(state, name='ccn_tuning', value=ccn_tuning, __RC__)

!   Aerosol mass mixing ratios
!   --------------------------
    mode_ = trim(mode)
    mode_ = ESMF_UtilStringLowerCase(mode_, __RC__)
    c_idx = index(mode_, '::', back=.true.)
    read (mode_(c_idx+2:c_idx+4),*) aerosol_bin

    allocate(q(i2,j2,k2),  __STAT__)
       
    q = 0.0

! ============================================================================
! CARMA version...
! Since only concentration elements are added to the ACI list, here we are
!    determining speciation solely based on group name.
! ============================================================================
    if (index(mode_, 'dust') > 0) then ! Pure dust
        do i = 1, size(aeroList)
           if (index(aeroList(i), 'dust') > 0 .and. index(aeroList(i), '::', back=.true.) == c_idx) then

              call ESMF_FieldBundleGet(aerosols, trim(aeroList(i)), field=fld, __RC__)
              call ESMF_FieldGet(fld, typekind=tk, __RC__)

              call ESMF_FieldGet(fld, farrayPtr=ptr_3d, __RC__)
              q = q + ptr_3d

              hygroscopicity = k_DU
              density = densDU
              ! Diameter of bin center
              diameter = 2*(((4./3.*MAPL_PI*density*aeroRmin(i)**3.) * &
                  aeroRmrat(i)**(aerosol_bin - 1))/(4./3.*MAPL_PI*density))**(1./3.)
              ! Geometric dispersion log(sigma_g)
              sigma(:,:,:) = log(aeroRmrat(i))/3./4.
              ! Number concentration to match mass in bin
              num = q / ((MAPL_PI/6.0) * density * diameter**3 * exp(4.5*sigma**2))
              ! Fraction of dust
              f_dust = 1.0
              ! Fraction of soot
              f_soot = 0.0
              ! Fraction of organic
              f_organic = 0.0
          end if
        end do
    else if (index(mode_, 'sulfate') > 0) then ! Pure sulfate
        do i = 1, size(aeroList)
           if (index(aeroList(i), 'sulfate') > 0 .and. index(aeroList(i), '::', back=.true.) == c_idx) then

              call ESMF_FieldBundleGet(aerosols, trim(aeroList(i)), field=fld, __RC__)
              call ESMF_FieldGet(fld, typekind=tk, __RC__)

              call ESMF_FieldGet(fld, farrayPtr=ptr_3d, __RC__)
              q = q + ptr_3d

              hygroscopicity = k_SO4
              density = densSO4 ! PAC: should this be the wet density?
              ! Diameter of bin center
              diameter = 2*(((4./3.*MAPL_PI*density*aeroRmin(i)**3.) * &
                  aeroRmrat(i)**(aerosol_bin - 1))/(4./3.*MAPL_PI*density))**(1./3.)
              ! Geometric dispersion log(sigma_g)
              sigma(:,:,:) = log(aeroRmrat(i))/3./4.
              ! Number concentration to match mass in bin
              num = q / ((MAPL_PI/6.0) * density * diameter**3 * exp(4.5*sigma**2))
              ! Fraction of dust
              f_dust = 0.0
              ! Fraction of soot
              f_soot = 0.0
              ! Fraction of organic
              f_organic = 0.0
          end if
        end do
    else if (index(mode_, 'seasalt') > 0) then ! Pure sea salt
        do i = 1, size(aeroList)
           if (index(aeroList(i), 'seasalt') > 0 .and. index(aeroList(i), '::', back=.true.) == c_idx) then

              call ESMF_FieldBundleGet(aerosols, trim(aeroList(i)), field=fld, __RC__)
              call ESMF_FieldGet(fld, typekind=tk, __RC__)

              call ESMF_FieldGet(fld, farrayPtr=ptr_3d, __RC__)
              q = q + ptr_3d

              hygroscopicity = k_SS
              density = densSS ! PAC: should this be the wet density?
              ! Diameter of bin center
              diameter = 2*(((4./3.*MAPL_PI*density*aeroRmin(i)**3.) * &
                  aeroRmrat(i)**(aerosol_bin - 1))/(4./3.*MAPL_PI*density))**(1./3.)
              ! Geometric dispersion log(sigma_g)
              sigma(:,:,:) = log(aeroRmrat(i))/3./4.
              ! Number concentration to match mass in bin
              num = q / ((MAPL_PI/6.0) * density * diameter**3 * exp(4.5*sigma**2))
              ! Fraction of dust
              f_dust = 0.0
              ! Fraction of soot
              f_soot = 0.0
              ! Fraction of organic
              f_organic = 0.0
          end if
        end do
    else if (index(mode_, 'organiccarbon') > 0) then ! Organic carbon
        do i = 1, size(aeroList)
           if (index(aeroList(i), 'organiccarbon') > 0 .and. index(aeroList(i), '::', back=.true.) == c_idx) then

              call ESMF_FieldBundleGet(aerosols, trim(aeroList(i)), field=fld, __RC__)
              call ESMF_FieldGet(fld, typekind=tk, __RC__)

              call ESMF_FieldGet(fld, farrayPtr=ptr_3d, __RC__)
              q = q + ptr_3d

              hygroscopicity = k_OC
              density = densOC ! PAC: should this be the wet density?
              ! Diameter of bin center
              diameter = 2*(((4./3.*MAPL_PI*density*aeroRmin(i)**3.) * &
                  aeroRmrat(i)**(aerosol_bin - 1))/(4./3.*MAPL_PI*density))**(1./3.)
              ! Geometric dispersion log(sigma_g)
              sigma(:,:,:) = log(aeroRmrat(i))/3./4.
              ! Number concentration to match mass in bin
              num = q / ((MAPL_PI/6.0) * density * diameter**3 * exp(4.5*sigma**2))
              ! Fraction of dust
              f_dust = 0.0
              ! Fraction of soot
              f_soot = 0.0
              ! Fraction of organic
              f_organic = 1.0
          end if
        end do
    else if (index(mode_, 'blackcarbon') > 0) then ! Black carbon
        do i = 1, size(aeroList)
           if (index(aeroList(i), 'blackcarbon') > 0 .and. index(aeroList(i), '::', back=.true.) == c_idx) then

              call ESMF_FieldBundleGet(aerosols, trim(aeroList(i)), field=fld, __RC__)
              call ESMF_FieldGet(fld, typekind=tk, __RC__)

              call ESMF_FieldGet(fld, farrayPtr=ptr_3d, __RC__)
              q = q + ptr_3d

              hygroscopicity = k_BC
              density = densBC ! PAC: should this be the wet density?
              ! Diameter of bin center
              diameter = 2*(((4./3.*MAPL_PI*density*aeroRmin(i)**3.) * &
                  aeroRmrat(i)**(aerosol_bin - 1))/(4./3.*MAPL_PI*density))**(1./3.)
              ! Geometric dispersion log(sigma_g)
              sigma(:,:,:) = log(aeroRmrat(i))/3./4.
              ! Number concentration to match mass in bin
              num = q / ((MAPL_PI/6.0) * density * diameter**3 * exp(4.5*sigma**2))
              ! Fraction of dust
              f_dust = 0.0
              ! Fraction of soot
              f_soot = 1.0
              ! Fraction of organic
              f_organic = 0.0
          end if
        end do
    else if (index(mode_, 'mixedp') > 0) then ! mixed group
        do i = 1, size(aeroList)
           if (index(aeroList(i), 'mixedp') > 0 .and. index(aeroList(i), '::', back=.true.) == c_idx) then

              call ESMF_FieldBundleGet(aerosols, trim(aeroList(i)), field=fld, __RC__)
              call ESMF_FieldGet(fld, typekind=tk, __RC__)

              call ESMF_FieldGet(fld, farrayPtr=ptr_3d, __RC__)
              q = q + ptr_3d

              hygroscopicity = k_SO4 ! PAC: We need to have some type of mixing here
              density = densSO4
              ! Diameter of bin center
              diameter = 2*(((4./3.*MAPL_PI*density*aeroRmin(i)**3.) * &
                  aeroRmrat(i)**(aerosol_bin - 1))/(4./3.*MAPL_PI*density))**(1./3.)
              ! Geometric dispersion log(sigma_g)
              sigma(:,:,:) = log(aeroRmrat(i))/3./4.
              ! Number concentration to match mass in bin
              num = q / ((MAPL_PI/6.0) * density * diameter**3 * exp(4.5*sigma**2))
              ! Fraction of dust PAC: need to speciate mixed group
              f_dust = 0.0
              ! Fraction of soot
              f_soot = 0.0
              ! Fraction of organic
              f_organic = 0.0
          end if
        end do
    end if ! Speciation
    ! ============================================================================

    deallocate(q, __STAT__)
    deallocate(aeroList, __STAT__)
    deallocate(aeroRmin, __STAT__)
    deallocate(aeroRmrat, __STAT__)

    RETURN_(ESMF_SUCCESS)

  end subroutine aerosol_activation_properties



 END MODULE CARMAchem_GridCompMod
