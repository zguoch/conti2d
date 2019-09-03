CC = g++-9 -fopenmp 
PROM = conti2d
binDIR=bin
CFLAGS = -I./fftw-3.3.7/build/include -I./netcdf/include
DEPS = $(shell find src -name "*.h")
SRC = $(shell find src -name "*.cpp")
OBJ = $(SRC:%.cpp=%.o)

$(PROM): $(OBJ)
	$(CC) -o $(binDIR)/$(PROM) $(OBJ) $(CFLAGS)  -L${PWD}/fftw-3.3.7/build/lib -lfftw3 -L./netcdf/lib -lnetcdf

%.o:%.cpp $(DEPS)
	$(CC) -c $< -o $@ $(CFLAGS) 

clean:
	rm -rf $(OBJ)