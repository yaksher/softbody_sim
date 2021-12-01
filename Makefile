CC = clang
CFLAGS = -g -pthread -Wall -Wextra -O3 -march=native# -fsanitize=address -fno-omit-frame-pointer
DEPS = solve_ivp.h kdtree.h vec_funcs.h constructors.h phys_types.h data.h initialize.h
OBJ = sim.o solve_ivp.o kdtree.o vec_funcs.o

clean:
	rm *.o bin/*

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

sim.out: $(OBJ)
	$(CC) $(CFLAGS) -o bin/$@ $^

run: sim.out
	bin/sim.out && python3 parse_spring_sim.py