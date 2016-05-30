##############################################################
#                                                            #
#        SYSTEM DEFINITION DEPENDENT OPTIONS                 #
#                                                            #
# These options are used to build a Debug (G), Release (O)   #
# versions of the code:                                      #
#                                                            #
#    * Fortran compiler and its options (FC,FFG,FFO)         #
#    * Linker and its options (LINK,LFG,LFO)                 #
#                                                            #
##############################################################
 
# Intel Fortran compiler and linker
ifeq ($(COMPILER),INTEL)

   OMP=-openmp

   FFG=-O0 -traceback -g -debug all -check all -ftrapuv 
   FFO=-O3 -opt-report-file=$(NAME).opt.log $(IPO) # -xHost

   FC=mpif90 -c -r8 -fpp -W1 -warn nointerfaces

   I=I
   J=module 

   LINK=mpif90 
   LFG=-g # -tcheck
   LFO=$(IPO)

   #-----------------------------------------------------------------------
   # Configure MKL_LIBS link line accordingly to 
   # user-prescribed environment variables
   # Link lines were extracted from the following
   # URL: http://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/
   #
   # - MKLROOT: absolute path, Intel MKL library installation directory.
   #
   # - MKLVERSION: version of the MKL library. Only two possible values: 
   #               "10.0-10.2" or "10.3".

   ifdef MKLROOT
     ifdef MKLVERSION
         ifeq ($(MKLVERSION),10.0-10.2)
            MKL_LIB_DIR = -L$(MKLROOT)/lib/em64t

            #32-bit integer, MKL multithreaded
            MKL_LIB_THREAD=-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

            #32-bit integer, MKL sequential
            MKL_LIB_SEQ=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
         else	
            ifeq ($(MKLVERSION),10.3)
              MKL_LIB_DIR = -L$(MKLROOT)/lib/intel64

              #32-bit integer, MKL multithreaded
              MKL_LIB_THREAD=-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -openmp -lpthread

              #32-bit integer, MKL sequential
              MKL_LIB_SEQ=-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
            endif
         endif 
       endif 
     endif
endif

#############################################################
# GNU Fortran compiler and linker
ifeq ($(COMPILER),GNU)
   OMP = -fopenmp

   # -Warray-temporaries -Wimplicit-interface -Waliasing
   FFG=-g -fbacktrace -fbounds-check
   FFO=-O3
   I=I
   J=J

   # -Wimplicit-interface 
   FC=mpif90 -c -fdefault-real-8 -ffree-line-length-0 -cpp -Wimplicit-interface 

   LINK=mpif90 -Warray-temporaries
   LFG=-g
   LFO=$(IPO)


   #-----------------------------------------------------------------------
   # Configure MKL_LIBS link line accordingly to 
   # user-prescribed environment variables
   # Link lines were extracted from the following
   # URL: http://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/
   #
   # - MKLROOT: absolute path, Intel MKL library installation directory.
   #
   # - MKLVERSION: version of the MKL library. Only two possible values: 
   #               "10.0-10.2" or "10.3". 
   ifdef MKLROOT
     ifdef MKLVERSION
            ifeq ($(MKLVERSION),10.0-10.2)
               MKL_LIB_DIR = $(MKLROOT)/lib/em64t
               #32-bit integer, MKL multithreaded, libiomp5 
               MKL_LIB_THREAD= -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -fopenmp -lpthread
               #32-bit integer, MKL sequential
               MKL_LIB_SEQ= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
            else
               ifeq ($(MKLVERSION),10.3)
                  MKL_LIB_DIR = $(MKLROOT)/lib/intel64
                  #32-bit integer, MKL multithreaded, libiomp5
                  MKL_LIB_THREAD= -lmkl_gf_lp64 -lmkl_gnu_thread -lmkl_core -fopenmp -lpthread
                  #32-bit integer, MKL sequential
                  MKL_LIB_SEQ= -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lpthread
               endif  
            endif
     endif
   endif

endif
#############################################################
