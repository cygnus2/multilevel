
 .SUFFIXES: .f90 .o .F90


 .f90.o:
	$(FC) -c -I.  $(FFLAGS) $<

 .F90.o:
	$(FC) -c $(FFLAGS) $<

FOBJ= headers.o ranlux.o init.o unitar.o su2-update.o matops.o storegauge.o addr.o staples.o bndstaples.o plaq.o plaqfull.o matops1.o su2-overrelax.o subrelax.o boundrelax.o boundsp-relax.o makePOL.o makecorr.o makecorr1.o b-update.o bs-update.o sub-update.o hit1.o ahit.o pol.o 

FC = mpif90
LD = mpif90
FFLAGS = -O
LDFLAGS = -L /opt/intel/mkl/lib/intel64/ 
#/opt/intel/mkl/lib/intel64/libmkl_rt.so

su3: $(FOBJ)
	$(LD)  $(FOBJ) -o $@  $(LDFLAGS)

ranlux.o : ranlux.f90
headers.o : headers.f90
addr.o : headers.o
storegauge.o : headers.o storegauge.f90
unitar.o : headers.o unitar.f90 init.o
init.o : headers.o init.f90  
matops.o : matops.f90 init.o
staples.o : staples.f90 matops.o headers.o addr.o
bndstaples.o : bndstaples.f90 matops.o headers.o addr.o
su2-overrelax.o : su2-overrelax.f90 headers.o 
subrelax.o : subrelax.f90 headers.o matops.o su2-overrelax.o addr.o
boundrelax.o : boundrelax.f90 headers.o matops.o su2-overrelax.o addr.o
boundsp-relax.o : boundsp-relax.f90 headers.o matops.o su2-overrelax.o addr.o
su2-update.o : su2-update.f90 headers.o ranlux.o 
b-update.o : b-update.f90 headers.o matops.o su2-update.o addr.o 
bs-update.o : bs-update.f90 headers.o matops.o su2-update.o addr.o
sub-update.o : sub-update.f90 headers.o matops.o su2-update.o addr.o
ahit.o : ahit.f90 headers.o  
hit1.o : hit1.f90 ahit.o headers.o 
makePOL.o : makePOL.f90 headers.o addr.o matops.o hit1.o init.o
plaq.o : plaq.f90 headers.o addr.o matops.o init.o
matops1.o : matops1.f90 init.o
plaqfull.o : plaqfull.f90 headers.o addr.o matops1.o init.o
makecorr.o : makecorr.f90 headers.o addr.o matops.o hit1.o init.o
makecorr1.o : makecorr1.f90 headers.o addr.o matops.o hit1.o init.o
pol.o :  pol.f90 init.o unitar.o ranlux.o cab-mar.o su2-update.o matops.o storegauge.o addr.o headers.o bndstaples.o su2-overrelax.o subrelax.o boundrelax.o boundsp-relax.o makePOL.o makecorr.o makecorr1.o b-update.o bs-update.o sub-update.o hit1.o ahit.o staples.o plaq.o
clean : 
	rm *.mod *.o su3
