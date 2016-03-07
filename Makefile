program = monte_carlo

obj = $(source:.cpp=.o)

headers = $(source:.cpp=.h)
source = Boundaries.cpp
source += Surface.cpp
source += main.cpp
source += Material.cpp
source += Tally.cpp
source += Neutron.cpp
source += Mesh.cpp
source += Monte_carlo.cpp
source += Plotter.cpp
source += Fission.cpp
source += ../../OpenMOC/src/Point.cpp
source += ../../OpenMOC/src/Universe.cpp
source += ../../OpenMOC/src/Python.cpp
source += ../../OpenMOC/src/LocalCoords.cpp
source += ../../OpenMOC/src/boundary_type.cpp
source += ../../OpenMOC/src/Material.cpp
source += ../../OpenMOC/src/linalg.cpp
source += ../../OpenMOC/src/Vector.cpp
source += ../../OpenMOC/src/Matrix.cpp

CFLAGS += -DFP_PRECISION=double

CC = g++

$(program): $(obj) $(headers)
	$(CC) $(obj) -o $@ -lm

%.o: %.cpp
	$(CC) -c $< -o $@

clean:
	rm -rf $(program) $(obj)

edit:
	vim -p $(source) $(headers)

run:
	./$(program)
