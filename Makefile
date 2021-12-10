CC = clang
CFLAGS = -g -pthread -Wall -Wextra -Ofast -march=native #-fsanitize=thread# -fsanitize=address
DEPS = solve_ivp.h kdtree.h vec_funcs.h constructors.h phys_types.h data.h initialize.h
OBJ = sim.o solve_ivp.o kdtree.o vec_funcs.o

%.o: %.c $(DEPS)
	$(CC) $(CFLAGS) -c -o $@ $<

sim.out: $(OBJ)
	$(CC) $(CFLAGS) -o bin/$@ $^

clean:
	rm *.o bin/*

run: sim.out
	bin/sim.out && python3 parse_spring_sim.py