# combinational

A combinational circuit simulator using optimized ternary analysis. Incorporates an elite genetic algorithm for evolving circuits towards a Boolean goal.

Compilation:
`gcc -Wall -DSFMT_MEXP=19937 -g -o combinational combinational.c mt/SFMT.c`

Usage:
`./combinational`

To change the number of inputs, outputs, and gates, change the macros near the top of combinational.c -- changing these will automatically update the genome encoding. The evolutionary parameters can additionally be manipulated by changing defined macros. To change the evolutionary goal, modify the eval_network_fitness_vector argument in the time_to_perfect function (currently set to the rivest function).
