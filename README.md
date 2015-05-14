# combinational

A combinational circuit simulator using optimized ternary analysis. Incorporates an elite genetic algorithm for evolving circuits towards a Boolean goal.

Compilation:
`gcc -Wall -DSFMT_MEXP=19937 -O3 -o combinational combinational.c mt/SFMT.c`

Usage:
`./combinational`

This program currently runs 50 experiments of modular Boolean goal evolution, evolving 4-input/1-output circuits towards two Boolean goal functions that share overlapping subproblems. Given input variables *a*, *b*, *c*, and *d*, these functions are `(a xor b) and (c xor d)` and `(a xor b) or (c xor d)`. This evolutionary process is described in depth in Nadav Kashtan and Uri Alon's paper *Spontaneous evolution of modularity and network motifs* (http://www.pnas.org/content/102/39/13773). This program attains comparable results to Kashtan and Alon but uses a rigorously correct circuit model.

To change the number of inputs, outputs, and gates, alter the macros near the top of combinational.c — these modifications will automatically update the genome encoding. The evolutionary parameters can additionally be manipulated by changing other defined macros at the top. To use a different evolutionary goal, modify the main function and create a goal function following the guide of functions `goal1` and `goal2` (number of inputs and outputs dervied from constants `INPUTS` and `OUTPUTS`).
