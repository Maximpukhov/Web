CXXFLAGS = -g -O2 -mavx512f
#CXXFLAGS = -g -msse2 -msse3 -msse4 -mmmx

all: mc_growth.exe

mc_growth.exe :  mc_growth.o geometry.o deform.o deform2.o deform3.o jump.o cycle.o output.o 
	g++ $(CXXFLAGS)  -Wall -lm -o mc_growth.exe  mc_growth.o geometry.o deform.o deform2.o deform3.o jump.o cycle.o output.o


clean :
	rm -f mc_growth.exe *.o
	rm -f A[0-9]\.[0-9]\.xyz
	rm -f A[0-9][0-9]\.[0-9]\.xyz
	rm -f initial.xyz
	rm -f M*\.txt
	rm -f core*
	rm -f _A[0-9]\.[0-9][0-9]\.xyz
	rm -f A[0-9]\.[0-9][0-9]\.xyz
	rm -f log.txt
	rm -f _initial.xyz

clean_dos :
	del mc_growth.exe *.o
