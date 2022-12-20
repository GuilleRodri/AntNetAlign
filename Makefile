TARGET = AntNetAlign_Neg AntNetAlign
CXXFLAGS = -ansi -O3 -fpermissive -std=c++17 

SYSTEM     = x86-64_linux
LIBFORMAT  = static_pic
CPLEXDIR      = /home/guillem/Desktop/CPLEX/CPLEX_Studio201/cplex
CONCERTDIR    = /home/guillem/Desktop/CPLEX/CPLEX_Studio201/concert
GCC = gcc
CCC = g++
CCOPT = -m64 -O -fPIC -fexceptions -DNDEBUG -DIL_STD -std=c++17 -fpermissive -w
CPLEXBINDIR   = $(CPLEXDIR)/bin/$(BINDIST)
CPLEXLIBDIR   = $(CPLEXDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CONCERTLIBDIR = $(CONCERTDIR)/lib/$(SYSTEM)/$(LIBFORMAT)
CCLNFLAGS = -L$(CPLEXLIBDIR) -lilocplex -lcplex -L$(CONCERTLIBDIR) -lconcert -lm -pthread -lpthread -ldl
CLNFLAGS  = -L$(CPLEXLIBDIR) -lcplex -lm -pthread -lpthread -ldl
CONCERTINCDIR = $(CONCERTDIR)/include
CPLEXINCDIR   = $(CPLEXDIR)/include

CCFLAGS = $(CCOPT) -I$(CPLEXINCDIR) -I$(CONCERTINCDIR)

all: ${TARGET}
	
AntNetAlign_Neg: AntNetAlign_Neg.o
	$(CCC) $(CCFLAGS) AntNetAlign_Neg.o -o AntNetAlign_Neg $(CCLNFLAGS)

AntNetAlign_Neg.o: AntNetAlign_Neg.cpp
	$(CCC) -c $(CCFLAGS) AntNetAlign_Neg.cpp -o AntNetAlign_Neg.o
	
AntNetAlign: AntNetAlign.o
	$(CCC) $(CCOPT) AntNetAlign.o -o AntNetAlign

AntNetAlign.o: AntNetAlign.cpp
	$(CCC) -c $(CCOPT) AntNetAlign.cpp -o AntNetAlign.o 
	
clean:
	@rm -f *~ *.o ${TARGET} core


