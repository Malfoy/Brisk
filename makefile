CC=g++

CFLAGS+=-std=c++14 -pipe -fopenmp -mavx2 ${WARNS}
CFLAGS+=-Os -fno-fat-lto-objects -flto=jobserver -march=native -mtune=native -mcmodel=large -msse3 -mavx

LDFLAGS+=-fuse-linker-plugin
LDFLAGS+=-lpthread -fopenmp -lz 

WARNS=-Wfatal-errors -Wall
        
EXEC=counter



all: $(EXEC) 


counter: counter.o buckets.o SuperKmerLight.o DenseMenuYo.o Kmers.o
	@echo "[LD] $@"
	+@$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

counter.o: counter.cpp 
	@echo "[CC] $<"
	@$(CC) $(CFLAGS) -MMD -o $@ -c $<


%.o: %.cpp 
	@echo "[CC] $<"
	@$(CC) $(CFLAGS) -MMD -o $@ -c $<



clean:
	@echo "[clean]"
	@rm -rf *.o $(EXEC)


rebuild: clean
	@$(MAKE) -s all
