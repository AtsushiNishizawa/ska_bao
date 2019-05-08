DIR = 
F77 = gfortran -O3
#F77 = /opt/intel/composer_xe_2015.1.108/bin/intel64/ifort -O3

B_term:			B_term.f Gamma1_1loop.f			
	$(F77)		B_term.f Gamma1_1loop.f \
			-o B_term.exe

A_term_I:		A_term_I.f \
			bispec_RegPT_I.f Gamma2d.f Gamma2v.f \
			Gamma1_1loop.f
	$(F77)		A_term_I.f \
			bispec_RegPT_I.f Gamma2d.f Gamma2v.f \
			Gamma1_1loop.f \
			-o A_term_I.exe

A_term_II:		A_term_II.f bispec_RegPT_II_IRsafe.f libcuba.a
	$(F77)		A_term_II.f bispec_RegPT_II_IRsafe.f -L./Cuba-4.2 -lcuba \
			-o A_term_II.exe

clean           :
		        rm -f *.*~ *.o *.exe
