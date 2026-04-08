esma_set_this ()

set (src_directories
  .
  CARMA
  CARMA/source/base
  )

set (srcs)
foreach (dir ${src_directories})
  file (GLOB tmpsrcs CONFIGURE_DEPENDS ${dir}/*.[fF] ${dir}/*.[fF]90 ${dir}/*.c)
  list (APPEND srcs ${tmpsrcs})
endforeach()

set (dependencies Chem_Shared Chem_Base GMAO_mpeu ESMF::ESMF)
esma_add_library (${this} SRCS ${srcs} DEPENDENCIES ${dependencies})
target_include_directories (${this} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/CARMA/source/base>)

#new_esma_generate_automatic_code (
#  ${this} CARMAchem_StateSpecs.rc
#  "CARMA_ExportSpec___.h;CARMA_ImportSpec___.h;CARMA_DeclarePointer___.h;CARMA_GetPointer___.h"
#  CARMA_History___.rc
#  ${include_GEOSchem_GridComp} ${esma_etc}
#   "-v;-F"
#  )
mapl_acg (${this}  CARMAchem_StateSpecs.rc
          IMPORT_SPECS EXPORT_SPECS INTERNAL_SPECS
          GET_POINTERS DECLARE_POINTERS)

set (resource_files
   CARMAchem_GridComp_ExtData.rc
   CARMAchem_GridComp_ExtData.yaml
   CARMAchem_MieRegistry.rc
   CARMAchem_Registry.rc
   CARMAchem_StateSpecs.rc
   )

install(
   FILES ${resource_files}
   DESTINATION etc
   )
