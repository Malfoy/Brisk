CC=g++


CFLAGS+=-std=c++14 -pipe -fopenmp -mavx2 ${WARNS}
CFLAGS+=-Os -fno-fat-lto-objects -flto=jobserver -march=native -mtune=native -mcmodel=large -msse3 -mavx


LDFLAGS+=-fuse-linker-plugin
LDFLAGS+=-lpthread -fopenmp -lz 



WARNS=-Wfatal-errors -Wall


        
EXEC=main



all: $(EXEC) 


main: main.o
	@echo "[LD] $@"
	+@$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

main.o: main.cpp 
	@echo "[CC] $<"
	@$(CC) $(CFLAGS) -MMD -o $@ -c $<


Brisk: Brisk.o
	@echo "[LD] $@"
	+@$(CC) -o $@ $^ $(CFLAGS) $(LDFLAGS)

Brisk.o: Brisk.cpp 
	@echo "[CC] $<"
	@$(CC) $(CFLAGS) -MMD -o $@ -c $<


%.o: %.cpp 
	@echo "[CC] $<"
	@$(CC) $(CFLAGS) -MMD -o $@ -c $<



clean:
	@echo "[clean]"
	@rm -rf *.o Brisk


rebuild: clean
	@$(MAKE) -s all
