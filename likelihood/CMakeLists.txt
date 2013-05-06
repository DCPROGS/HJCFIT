add_library(likelihood SHARED errors.h idealG.h idealG.cc equilibrium.h equilibrium.cc)
if(tests) 
  add_subdirectory(tests)
endif(tests)
