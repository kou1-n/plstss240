#.SUFFIXES:
TARGET  =  $(HOME)/bin/plstss2393_nagasaku
FC = ifort
#
F90 = $(FC)
#
##### MKL for EM64T #####
MKLPATH =$(MKLROOT)/lib/em64t
MKLINCLUDE = $(MKLROOT)/include
MKL= -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp -lpthread
MKL= -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp
MKL= -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp -lpthread
#

##### for Debugging #####
#FFLAGS = -O0 -C
FFLAGS = -O3 -C
#FFLAGS = -O3
#
#
#( for Core i7 (AVX2) )
#FFLAGS = -O3 -xcore-avx512 -qopenmp -parallel
#FFLAGS = -O3 -xcore-avx512 -qopenmp
##### for General settings #####
#( for Core i7 (AVX) ) for HONBAN
#FFLAGS = -O3 -xavx -static
#
#( for Core i7 (SSE4.2) )
#FFLAGS = -O3 -xsse4.2 -parallel -static
#FFLAGS = -O3 -xsse4.2 -parallel
#FFLAGS = -O3 -xsse4.2 -static
#
#FFLAGS = -O3 -xhost
#
#( for Windows Compilers )
#FFLAGS = -O3 /Qvec-report3 /C
#FFLAGS = -O3 /Qvec-report3
#
#( for Core 2 Duo/Quad (SSE4.1) )
#FFLAGS = -O3 -xsse4.1 -parallel
#FFLAGS = -O3 -xsse4.1
#FFLAGS = -O3 -msse4.1
#
#
F90FLAGS = $(FFLAGS)
LINKER = $(FC)
# --------------------------------------
OBJMAIN  = main.o flopen.o analys.o redcml.o inv_22.o inv_33.o zerodt.o gaussj.o dgaussj.o
OBJPRE   = chkary.o estima.o presky.o addres.o inform.o initia.o initq4.o initt3.o inith8.o initq8.o initt6.o initt4.o initp1.o 
OBJPARDIS= prepar.o rowupr.o parsol.o pars00.o pars99.o pcgsol.o 
OBJASMB  = forces.o assemb.o elastc.o plastc.o mapcrs.o mapsky.o 
OBJELM   = quad4a.o pquad4.o hexa8a.o phexa8.o tria3a.o ptria3.o pres2d.o pres3d.o triap1.o ptrip1.o tria6a.o ptria6.o quad8a.o pquad8.o
#OBJBASE  = trans1.o trans2.o trans3.o trans4.o
OBJSLV   = constr.o skylin.o
OBJPOST  = update.o postpr.o output.o stored.o restor.o stress.o st_gtn.o
#
#OBJRCM   = genrcm.o
#
# --------------------------------------
OBJALL = $(OBJMAIN) $(OBJPRE) $(OBJPARDIS) $(OBJASMB) $(OBJELM) $(OBJSLV) $(OBJPOST)
# $(OBJBASE)
# --------------------------------------
$(TARGET):$(OBJALL)
	$(FC) $(FFLAGS) $(MKL) $(OBJALL) -o $(TARGET) $(MKL)
	@echo make PLSTss Version 2.2 done.
#
clean:
	rm -f $(TARGET) $(OBJALL)
#
#.SUFFIXES: $(SUFFIXES) .f90
#.f90.o:
#	$(F90) $(F90FLAGS) -c $<
#
.SUFFIXES: $(SUFFIXES) .f
.f.o:
#	$(FC) $(FFLAGS) $(MKL) -c $<
#	$(FC) $(FFLAGS) -c $<
#	$(F90) $(F90LAGS) /compile_only /object:test.o  $<
#	$(F90) $(F90LAGS) /c /object:$*.o  $<
#

