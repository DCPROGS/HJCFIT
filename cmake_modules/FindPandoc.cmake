find_program(PANDOC_EXECUTABLE NAMES pandoc
  HINTS
  $ENV{PANDOC_DIR}
  PATH_SUFFIXES bin
  DOC "Pandoc text converter"
)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Pandoc DEFAULT_MSG
  PANDOC_EXECUTABLE
)

mark_as_advanced(
  PANDOC_EXECUTABLE
)
