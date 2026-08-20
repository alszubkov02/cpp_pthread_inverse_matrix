OPTS = -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
a.out: main.o func.o jordan.o
	g++ $(OPTS) main.o func.o jordan.o -o a.out
main.o: main.cpp func.h
	g++ -c $(OPTS) main.cpp -lm
func.o: func.cpp func.h
	g++ -c $(OPTS) func.cpp
jordan.o: jordan.cpp func.h
	g++ -c $(OPTS) jordan.cpp
clean:
	rm -f main.o func.o jordan.o a.out
