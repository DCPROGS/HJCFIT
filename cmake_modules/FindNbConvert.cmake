find_program(NBCONVERT_EXECUTABLE NAMES jupyter-nbconvert
  HINTS
  $ENV{NBCONVERT_DIR}
  PATH_SUFFIXES bin
  DOC "Jupyter nbconvert"
)

include(FindPackageHandleStandardArgs)

find_package_handle_standard_args(Nbconvert DEFAULT_MSG
  NBCONVERT_EXECUTABLE
)

mark_as_advanced(
  NBCONVERT_EXECUTABLE
)
