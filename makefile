# makefile
# Tabs *must* be used for the indentations below;
# spaces cause syntax errors.

CC=g++
LIBS=-lm
GSLLIBS=-lgsl -lgslcblas 



mlesearch:		
		$(CC) -O3 -c -o readdata.o readdata.cpp $(LIBS) $(GSLLIBS)
		$(CC) -O3 -c -o likelihood.o likelihood.cpp $(LIBS) $(GSLLIBS)
		$(CC) -O3 -c -o prms.o prms.cpp $(LIBS) $(GSLLIBS)
		$(CC) -O3 -c -o main_mlesearch_and_profile.o main_mlesearch_and_profile.cpp $(LIBS) $(GSLLIBS)
		$(CC) $(CFLAGS) -o mlesearch_and_profile main_mlesearch_and_profile.o readdata.o prms.o likelihood.o $(LIBS) $(GSLLIBS) 

mcmc:		
		$(CC) -O3 -c -o readdata.o readdata.cpp $(LIBS) $(GSLLIBS)
		$(CC) -O3 -c -o likelihood.o likelihood.cpp $(LIBS) $(GSLLIBS)
		$(CC) -O3 -c -o prms.o prms.cpp $(LIBS) $(GSLLIBS)
		$(CC) -O3 -c -o chain.o chain.cpp $(LIBS) $(GSLLIBS)
		$(CC) -O3 -c -o main_Senterica_mcmc.o main_Senterica_mcmc.cpp $(LIBS) $(GSLLIBS)
		$(CC) $(CFLAGS) -o mcmc_topt main_Senterica_mcmc.o chain.o readdata.o prms.o likelihood.o $(LIBS) $(GSLLIBS) 
		

clean:
		rm -f *.o core mlesearch_and_profile mcmc_topt




