program = monte_carlo

DEBUG = yes

headers = $(source:.cpp=.h)
source = main.cpp
source += Tally.cpp
source += Neutron.cpp
source += MCSolver.cpp
source += Plotter.cpp
source += Fission.cpp
source += Enumerations.cpp
source += ../../OpenMOC/src/Point.cpp
source += ../../OpenMOC/src/Universe.cpp
source += ../../OpenMOC/src/Surface.cpp
source += ../../OpenMOC/src/LocalCoords.cpp
source += ../../OpenMOC/src/Cell.cpp
source += ../../OpenMOC/src/log.cpp
source += ../../OpenMOC/src/Material.cpp
source += ../../OpenMOC/src/Geometry.cpp
source += ../../OpenMOC/src/linalg.cpp
source += ../../OpenMOC/src/Vector.cpp
source += ../../OpenMOC/src/Matrix.cpp
source += ../../OpenMOC/src/Cmfd.cpp
source += ../../OpenMOC/src/Track.cpp
source += ../../OpenMOC/src/PolarQuad.cpp
source += ../../OpenMOC/src/ExpEvaluator.cpp
source += ../../OpenMOC/src/TrackGenerator.cpp
source += ../../OpenMOC/src/Solver.cpp
source += ../../OpenMOC/src/Timer.cpp

obj = $(source:.cpp=.o)

CC = g++

CFLAGS := -DFP_PRECISION=double
CFLAGS += -DVEC_LENGTH=8
CFLAGS += -fopenmp
CFLAGS += -std=c++11
CFLAGS += -DOPENMP
ifeq ($(DEBUG),yes)
	CFLAGS += -g
endif

$(program): $(obj) $(headers)
	$(CC) $(CFLAGS) $(obj) -o $@ -lm

%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(program) $(obj)

edit:
	vim -p $(source) $(headers)

run:
	./$(program)
