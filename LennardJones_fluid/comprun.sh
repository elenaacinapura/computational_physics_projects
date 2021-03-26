gcc fluid.c util.c FCC.c -lm -o fluid.out -O3 -march=native -ffast-math
./fluid.out > output.txt
