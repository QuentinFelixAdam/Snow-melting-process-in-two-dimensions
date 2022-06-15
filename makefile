# This is an commentary line in a makefile
# Start of the makefile
objects = m_oalloc.o m_glob.o temp.o sub_init.o sub_output.o sub_NTO.o sub_readpara.o sub_annex.o
exe: $(objects)
	f90 -o exe $(objects)
m_oalloc.mod: m_oalloc.o m_oalloc.f
	f90 -free -c m_oalloc.f
m_oalloc.o: m_oalloc.f
	f90 -free -c m_oalloc.f
m_glob.mod: m_glob.o m_glob.f
	f90 -free -c m_glob.f
m_glob.o: m_glob.f
	f90 -free -c m_glob.f
temp.o: m_glob.mod temp.f
	f90 -free -c -xopenmp  temp.f
sub_init.o: sub_init.f
	f90 -free -c sub_init.f
sub_output.o: sub_output.f
	f90 -free -c sub_output.f
sub_NTO.o: sub_NTO.f
	f90 -free -c sub_NTO.f
sub_readpara.o: sub_readpara.f
	f90 -free -c sub_readpara.f
sub_annex.o: sub_annex.f
	f90 -free -c sub_annex.f
clean:
	rm $(objects) exe m_glob.mod m_oalloc.mod 
# End of the makefile
