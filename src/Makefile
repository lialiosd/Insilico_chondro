COMPILER_FLAGS = -qopenmp -qmkl=parallel

app : main.cc
	icx -L$(LD_LIBRARY_PATH) -I$(MKLROOT)/include main.cc -o app $(COMPILER_FLAGS)

.PHONY : clean
	rm app
