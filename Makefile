#.SUFFIXES:
TARGET  =  $(HOME)/bin/plstss2393_nagasaku
FC = ifort
# Directory where object files are placed
OBJDIR = obj

F90 = $(FC)
##### MKL for EM64T #####
MKLPATH =$(MKLROOT)/lib/em64t
MKLINCLUDE = $(MKLROOT)/include
#MKL= -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp -lpthread
#MKL= -Wl,--start-group  $(MKLROOT)/lib/intel64/libmkl_intel_lp64.a $(MKLROOT)/lib/intel64/libmkl_intel_thread.a $(MKLROOT)/lib/intel64/libmkl_core.a -Wl,--end-group -qopenmp
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
OBJMAIN  = $(OBJDIR)/main.o $(OBJDIR)/flopen.o $(OBJDIR)/analys.o $(OBJDIR)/redcml.o $(OBJDIR)/inv_22.o \
           $(OBJDIR)/inv_33.o $(OBJDIR)/zerodt.o $(OBJDIR)/gaussj.o $(OBJDIR)/dgaussj.o
OBJPRE   = $(OBJDIR)/chkary.o $(OBJDIR)/estima.o $(OBJDIR)/presky.o $(OBJDIR)/addres.o $(OBJDIR)/inform.o \
           $(OBJDIR)/initia.o $(OBJDIR)/initq4.o $(OBJDIR)/initt3.o $(OBJDIR)/inith8.o $(OBJDIR)/initq8.o \
           $(OBJDIR)/initt6.o $(OBJDIR)/initt4.o $(OBJDIR)/initp1.o
OBJPARDIS= $(OBJDIR)/prepar.o $(OBJDIR)/rowupr.o $(OBJDIR)/parsol.o $(OBJDIR)/pars00.o $(OBJDIR)/pars99.o $(OBJDIR)/pcgsol.o
OBJASMB  = $(OBJDIR)/forces.o $(OBJDIR)/assemb.o $(OBJDIR)/elastc.o $(OBJDIR)/plastc.o $(OBJDIR)/mapcrs.o $(OBJDIR)/mapsky.o
OBJELM   = $(OBJDIR)/quad4a.o $(OBJDIR)/pquad4.o $(OBJDIR)/hexa8a.o $(OBJDIR)/phexa8.o $(OBJDIR)/tria3a.o \
           $(OBJDIR)/ptria3.o $(OBJDIR)/pres2d.o $(OBJDIR)/pres3d.o $(OBJDIR)/triap1.o $(OBJDIR)/ptrip1.o \
           $(OBJDIR)/tria6a.o $(OBJDIR)/ptria6.o $(OBJDIR)/quad8a.o $(OBJDIR)/pquad8.o
#OBJBASE  = $(OBJDIR)/trans1.o $(OBJDIR)/trans2.o $(OBJDIR)/trans3.o $(OBJDIR)/trans4.o
OBJSLV   = $(OBJDIR)/constr.o $(OBJDIR)/skylin.o
OBJPOST  = $(OBJDIR)/update.o $(OBJDIR)/postpr.o $(OBJDIR)/output.o $(OBJDIR)/stored.o $(OBJDIR)/restor.o \
           $(OBJDIR)/stress.o $(OBJDIR)/stress_vm.o $(OBJDIR)/st_gtn.o \
           $(OBJDIR)/stress_dp.o $(OBJDIR)/stress_dp_rm.o $(OBJDIR)/hardfunc.o
#
#OBJRCM   = genrcm.o
#
# --------------------------------------
OBJALL = $(OBJMAIN) $(OBJPRE) $(OBJPARDIS) $(OBJASMB) $(OBJELM) $(OBJSLV) $(OBJPOST)
# $(OBJBASE)
# --------------------------------------
$(TARGET): $(OBJALL)
	$(FC) $(FFLAGS) $(MKL) $(OBJALL) -o $(TARGET) $(MKL)
	@echo make PLSTss Version 2.2 done.

# ensure object directory exists
$(OBJDIR):
	mkdir -p $(OBJDIR)

clean:
	rm -f $(TARGET) $(OBJALL)
	rm -rf $(OBJDIR)

# Pattern rule to build object files in $(OBJDIR)
.SUFFIXES: $(SUFFIXES) .f
$(OBJDIR)/%.o: %.f | $(OBJDIR)
	$(FC) $(FFLAGS) -c $< -o $@

