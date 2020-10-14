# Brendon Phillips
# PhD Candidate
# Bauch computational epidemiology research group
# Department of Applied Mathematics
# Faculty of Mathematics
# Universiity of Waterloo

include Makefile.config


## Main application file
DEMOS = \
	Simulation_hesitance \
	## simple
	## node_test_2
	\

all: $(DEMOS)

# COMPILE

%: %.cpp $(CSNAP)/Snap.o
	$(CC) $(CXXFLAGS) -g -o Skynet_hesitance $@.cpp -std=c++17 -lstdc++fs $(CSNAP)/Snap.o -I$(CSNAP) -I$(CGLIB) $(LDFLAGS) $(LIBS) ## 2> compilation.output

$(CSNAP)/Snap.o:
	$(MAKE) -C $(CSNAP)

clean:
	rm -f *.o $(DEMOS) *.exe
	rm -rf Debug Release
	rm -rf *.Err demo*.dat
