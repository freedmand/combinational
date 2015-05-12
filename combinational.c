#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <string.h>

#include "mt/SFMT.h"

typedef struct gate gate;

#define INDETERMINATE -1
#define MAX_FAN_IN 2

#define SEED 18

#define INPUTS 3
#define GATES 20
#define INPUTS_PER_GATE 2
#define OUTPUTS 6
#define DNA_LENGTH (GATES * INPUTS_PER_GATE + OUTPUTS)

#define MUTATION 0.7
#define TRIALS 10000
#define EPOCH 20
#define CIRCUITS 2000

#define ELITE 800

#define EXPERIMENTS 50

#define CACHE_UNDEFINED -2
#define MAX_CACHE 8

#define DEGREE 11
#define DEGREE_PENALTY 0.2

struct gate {
  int (*fn)(int*);
  int fan_in;
  gate** inputs;
  int num_inputs;
  int input_array_size;
  gate** outputs;
  int num_outputs;
  int output_array_size;
  int value;
  int* cache;
  int cache_size;
};

typedef struct {
  gate** gates;
  int num_gates;
  gate** inputs;
  int num_inputs;
  gate** output;
  int num_outputs;
} network;

int subset(int i, int j) {
  if (i == INDETERMINATE) {
    return 1;
  } else {
    return i == j;
  }
}

int subset_v(int* i, int* j, int length) {
  int x;
  for (x = 0; x < length; x++) {
    if (!subset(i[x], j[x])) {
      return 0;
    }
  }
  return 1;
}

int subset_v_gates(gate** i, int* j, int length) {
  int x;
  for (x = 0; x < length; x++) {
    if (!subset(i[x]->value, j[x])) {
      return 0;
    }
  }
  return 1;
}

int partial_join(int i, int j) {
  if (i == j && (i == 0 || i == 1)) {
    return i;
  } else if ((i == 0 || i == 1) && j == INDETERMINATE) {
    return i;
  } else if ((j == 0 || j == 1) && i == INDETERMINATE) {
    return j;
  } else {
    return INDETERMINATE;
  }
}

void partial_join_v(int* result, int* l1, int* l2, int length) {
  int i;
  for (i = 0; i < length; i++) {
    result[i] = partial_join(l1[i], l2[i]);
  }
}

int ternary_ext(int (*gate)(int*), int fan_in, int* inputs) {
  int check0 = 1;
  int check1 = 1;
  int i, j;

  static int bin_input[MAX_FAN_IN];
  int max_val = 1 << fan_in;
  for (i = 0; i < max_val; i++) {
    for (j = 0; j < fan_in; j++) {
      bin_input[j] = (i >> j) & 1;
    }
    if (subset_v(inputs, bin_input, fan_in)) {
      int output = gate(bin_input);
      if (output == 0) {
        check1 = 0;
      } else if (output == 1) {
        check0 = 0;
      }
    }
  }
  if (check0) {
    return 0;
  } else if (check1) {
    return 1;
  } else {
    return INDETERMINATE;
  }
}

int ternary_ext_gates(int (*gate)(int*), int fan_in, struct gate** inputs) {
  int check0 = 1;
  int check1 = 1;
  int i, j;

  static int bin_input[MAX_FAN_IN];
  int max_val = 1 << fan_in;
  for (i = 0; i < max_val; i++) {
    for (j = 0; j < fan_in; j++) {
      bin_input[j] = (i >> j) & 1;
    }
    if (subset_v_gates(inputs, bin_input, fan_in)) {
      int output = gate(bin_input);
      if (output == 0) {
        check1 = 0;
      } else if (output == 1) {
        check0 = 0;
      }
    }
  }
  if (check0) {
    return 0;
  } else if (check1) {
    return 1;
  } else {
    return INDETERMINATE;
  }
}

void gen_truth_table(int* result, int (*gate)(int*), int fan_in, int ternary) {
  int i, j;

  static int bin_input[MAX_FAN_IN];
  int max_val = 1;
  for (j = 0; j < fan_in; j++) {
    max_val *= (ternary ? 3 : 2);
  }
  for (i = 0; i < max_val; i++) {
    int div_val = i;
    for (j = fan_in - 1; j >= 0; j--) {
      bin_input[j] = div_val % (ternary ? 3 : 2);
      bin_input[j] = (bin_input[j] == 2 ? INDETERMINATE : bin_input[j]);
      div_val /= (ternary ? 3 : 2);
    }
    if (ternary) {
      result[i] = ternary_ext(gate, fan_in, bin_input);
    } else {
      result[i] = gate(bin_input);
    }
  }
}

void connect(gate* g1, gate* g2) {
  int i;
  int not_in = 1;
  // for (i = 0; i < g2->num_inputs; i++) {
  //   if (g1 == g2->inputs[i]) {
  //     not_in = 0;
  //   }
  // }
  // if (not_in) {
    if (g2->num_inputs + 1 > g2->input_array_size) {
      g2->input_array_size *= 2;
      g2->inputs = (gate**)realloc(g2->inputs, g2->input_array_size * sizeof(gate*));
    }
    g2->inputs[g2->num_inputs++] = g1;
  // }

  // not_in = 1;
  for (i = 0; i < g1->num_outputs; i++) {
    if (g2 == g1->outputs[i]) {
      not_in = 0;
    }
  }
  if (not_in) {
    if (g1->num_outputs + 1 > g1->output_array_size) {
      g1->output_array_size *= 2;
      g1->outputs = (gate**)realloc(g1->outputs, g1->output_array_size * sizeof(gate*));
    }
    g1->outputs[g1->num_outputs++] = g2;
  }
}

#define CACHE_BASE 3

static inline int cache_val(int fan_in, int* vals) {
  int v = 0;
  int multiplier = 1;
  int i;
  for (i = 0; i < fan_in; i++) {
    v += (vals[i] == INDETERMINATE ? 2 : vals[i]) * multiplier;
    multiplier *= CACHE_BASE;
  }
  return v;
}

static inline int cache_val_gate(int fan_in, gate** vals) {
  int v = 0;
  int multiplier = 1;
  int i;
  for (i = 0; i < fan_in; i++) {
    v += (vals[i]->value == INDETERMINATE ? 2 : vals[i]->value) * multiplier;
    multiplier *= CACHE_BASE;
  }
  return v;
}

int eval_gate(gate* g) {
  int cache_idx = cache_val_gate(g->num_inputs, g->inputs);
  if (g->cache[cache_idx] != CACHE_UNDEFINED) {
    g->value = g->cache[cache_idx];
    return g->cache[cache_idx];
  }

  int result = ternary_ext_gates(g->fn, g->fan_in, g->inputs);
  g->value = result;
  g->cache[cache_idx] = result;
  return result;
}

int eval_gate_inp(gate* g, int* vals) {
  int cache_idx = cache_val(g->fan_in, vals);
  if (g->cache[cache_idx] != CACHE_UNDEFINED) {
    g->value = g->cache[cache_idx];
    return g->cache[cache_idx];
  }

  int result = ternary_ext(g->fn, g->fan_in, vals);
  g->value = result;
  g->cache[cache_idx] = result;
  return result;
}

int and_g(int* vals) {
  return vals[0] && vals[1];
}

int or_g(int* vals) {
  return vals[0] || vals[1];
}

int xor_g(int* vals) {
  return vals[0] ^ vals[1];
}

int not_g(int* vals) {
  return !vals[0];
}

int nand_g(int* vals) {
  return !(vals[0] && vals[1]);
}

int input_g(int* vals) {
  return vals[0];
}

int goal1(int* inputs) {
  return (inputs[0] ^ inputs[1]) & (inputs[2] ^ inputs[3]);
}

int goal2(int* inputs) {
  return (inputs[0] ^ inputs[1]) | (inputs[2] ^ inputs[3]);
}

void rivest(int* output, int* inputs) {
  int total = 0;
  int multiplier = 1;
  int i;
  for (i = 0; i < 3; i++) {
    total += inputs[i] * multiplier;
    multiplier *= 2;
  }

  if (total == 0) {
    output[0] = 0;
    output[1] = 0;
    output[2] = 0;
    output[3] = 0;
    output[4] = 0;
    output[5] = 0;
  } else if (total == 1) {
    output[0] = 0;
    output[1] = 0;
    output[2] = 0;
    output[3] = 0;
    output[4] = 0;
    output[5] = 1;
  } else if (total == 2) {
    output[0] = 0;
    output[1] = 1;
    output[2] = 0;
    output[3] = 0;
    output[4] = 0;
    output[5] = 0;
  } else if (total == 3) {
    output[0] = 0;
    output[1] = 1;
    output[2] = 1;
    output[3] = 1;
    output[4] = 1;
    output[5] = 1;
  } else if (total == 4) {
    output[0] = 0;
    output[1] = 0;
    output[2] = 0;
    output[3] = 1;
    output[4] = 0;
    output[5] = 0;
  } else if (total == 5) {
    output[0] = 1;
    output[1] = 1;
    output[2] = 1;
    output[3] = 1;
    output[4] = 0;
    output[5] = 1;
  } else if (total == 6) {
    output[0] = 1;
    output[1] = 1;
    output[2] = 0;
    output[3] = 1;
    output[4] = 1;
    output[5] = 1;
  } else if (total == 7) {
    output[0] = 1;
    output[1] = 1;
    output[2] = 1;
    output[3] = 1;
    output[4] = 1;
    output[5] = 1;
  } else {
    assert(0);
  }
}

void digital_display(int* output, int* inputs) {
  int total = 0;
  int multiplier = 1;
  int i;
  for (i = 0; i < 3; i++) {
    total += inputs[i] * multiplier;
    multiplier *= 2;
  }

  if (total == 0) {
    output[0] = 1;
    output[1] = 1;
    output[2] = 1;
    output[3] = 0;
    output[4] = 1;
    output[5] = 1;
    output[6] = 1;
  } else if (total == 1) {
    output[0] = 0;
    output[1] = 0;
    output[2] = 0;
    output[3] = 0;
    output[4] = 0;
    output[5] = 1;
    output[6] = 1;
  } else if (total == 2) {
    output[0] = 0;
    output[1] = 1;
    output[2] = 1;
    output[3] = 1;
    output[4] = 1;
    output[5] = 1;
    output[6] = 0;
  } else if (total == 3) {
    output[0] = 0;
    output[1] = 0;
    output[2] = 1;
    output[3] = 1;
    output[4] = 1;
    output[5] = 1;
    output[6] = 1;
  } else if (total == 4) {
    output[0] = 1;
    output[1] = 0;
    output[2] = 0;
    output[3] = 1;
    output[4] = 0;
    output[5] = 1;
    output[6] = 1;
  } else if (total == 5) {
    output[0] = 1;
    output[1] = 0;
    output[2] = 1;
    output[3] = 1;
    output[4] = 1;
    output[5] = 0;
    output[6] = 1;
  } else if (total == 6) {
    output[0] = 1;
    output[1] = 1;
    output[2] = 0;
    output[3] = 1;
    output[4] = 1;
    output[5] = 0;
    output[6] = 1;
  } else if (total == 7) {
    output[0] = 0;
    output[1] = 0;
    output[2] = 1;
    output[3] = 0;
    output[4] = 0;
    output[5] = 1;
    output[6] = 1;
  } else {
    assert(0);
  }
}

void make_gate(gate* g, int (*fn)(int*), int fan_in) {
  gate** inputs = (gate**)malloc(sizeof(gate*) * fan_in);
  gate** outputs = (gate**)malloc(sizeof(gate*));

  g->cache_size = 1;
  int i;
  for (i = 0; i < fan_in; i++) {
    g->cache_size *= 3;
  }
  int* cache = (int*)malloc(sizeof(int) * g->cache_size);
  for (i = 0; i < g->cache_size; i++) {
    cache[i] = CACHE_UNDEFINED;
  }
  g->cache = cache;

  // g = (gate*)malloc(sizeof(gate));
  g->fn = fn;
  g->fan_in = fan_in;
  g->inputs = inputs;
  g->num_inputs = 0;
  g->input_array_size = fan_in;
  g->outputs = outputs;
  g->num_outputs = 0;
  g->output_array_size = 1;
  g->value = INDETERMINATE;
  // g->value = 0;
}

void reset_gate(gate* g) {
  g->num_inputs = 0;
  g->num_outputs = 0;
  g->value = INDETERMINATE;
  // g->value = 0;
}

int degree(int* cyclic, network* n) {
  int i, j;
  for (i = 0; i < n->num_gates; i++) {
    n->gates[i]->value = INDETERMINATE;
  }

  for (i = 0; i < n->num_inputs; i++) {
    n->inputs[i]->value = INDETERMINATE;
  }

  static gate* gates_encountered[GATES];
  static gate* next_layer[GATES];

  int gates_encountered_size = n->num_outputs;

  int effective = 0;
  for (i = 0; i < n->num_outputs; i++) {
    gates_encountered[i] = n->output[i];
    n->output[i]->value = 1;
    effective++;
  }

  printf("Hi\n");
  *cyclic = 0;

  while (1) {
    int next_layer_size = 0;
    for (i = 0; i < gates_encountered_size; i++) {
      for (j = 0; j < gates_encountered[i]->num_inputs; j++) {
        if ((gates_encountered[i]->inputs[j]->value == INDETERMINATE || gates_encountered[i]->inputs[j]->value == 0) && gates_encountered[i]->inputs[j] != gates_encountered[i]) {
          next_layer[next_layer_size++] = gates_encountered[i]->inputs[j];
          gates_encountered[i]->inputs[j]->value = 0;
          effective++;
        } else {
          printf("%d\n", gates_encountered[i]->inputs[j]->value);
          printf("Wabam!\n");
          *cyclic = 1;
        }
      }
      for (j = 0; j < gates_encountered[i]->num_inputs; j++) {
        if (gates_encountered[i]->inputs[j]->value == 0) {
          gates_encountered[i]->inputs[j]->value = 1;
        }
      }
    }
    for (i = 0; i < next_layer_size; i++) {
      gates_encountered[i] = next_layer[i];
    }
    gates_encountered_size = next_layer_size;

    if (gates_encountered_size == 0) {
      return effective;
    }
  }
}

int eval_network(int* output, network* n, int* vals) {
  int i;
  for (i = 0; i < n->num_gates; i++) {
    n->gates[i]->value = INDETERMINATE;
  }
  for (i = 0; i < n->num_inputs; i++) {
    n->inputs[i]->value = vals[i];
  }

  // int* gates_to_eval = (int*)malloc(sizeof(int) * n->num_gates);
  static int gates_to_eval[GATES];
  for (i = 0; i < n->num_gates; i++) {
    gates_to_eval[i] = 1;
  }

  while (1) {
    int all_indeterminate = 1;
    int new_gates = 0;
    for (i = 0; i < n->num_gates; i++) {
      if (gates_to_eval[i]) {
        if (eval_gate(n->gates[i]) == INDETERMINATE) {
          new_gates++;
        } else {
          gates_to_eval[i] = 0;
          all_indeterminate = 0;
        }
      }
    }
    if (all_indeterminate) {
      output[0] = INDETERMINATE;
      // output[0] = 0;
      return 0;
    }
    if (new_gates == 0) {
      break;
    }
  }

  for (i = 0; i < n->num_outputs; i++) {
    output[i] = n->output[i]->value;
  }
  return 1;
}

void eval_network_all(int** outputs, network* n) {
  int i,j;
  int max_val = 1 << n->num_inputs;

  int* bin_input = (int*)malloc(sizeof(int) * n->num_inputs);
  int* output = (int*)malloc(sizeof(int) * n->num_outputs);

  for (i = 0; i < max_val; i++) {
    for (j = 0; j < n->num_inputs; j++) {
      bin_input[n->num_inputs - j - 1] = (i >> j) & 1;
    }
    eval_network(output, n, bin_input);
    outputs[i] = output;
    output = (int*)malloc(sizeof(int) * n->num_outputs);
  }
  free(bin_input);
}

double eval_network_fitness(network* n, int (*fn)(int*)) {
  assert(n->num_outputs == 1);
  int i,j;
  int max_val = 1 << n->num_inputs;

  // int* bin_input = (int*)malloc(sizeof(int) * n->num_inputs);
  static int bin_input[INPUTS];
  int output[1];

  int total = 0;
  int correct = 0;
  for (i = 0; i < max_val; i++) {
    for (j = 0; j < n->num_inputs; j++) {
      bin_input[n->num_inputs - j - 1] = (i >> j) & 1;
    }
    if (eval_network(output, n, bin_input)) {
      if (output[0] == fn(bin_input)) {
        correct++;
      }
    }
    total++;
  }
  // free(bin_input);
  // free(output);

  return (double)correct / total;
}

double eval_network_fitness_vector(network* n, void (*fn)(int*, int*)) {
  int i,j;
  int max_val = 1 << n->num_inputs;

  // int* bin_input = (int*)malloc(sizeof(int) * n->num_inputs);
  static int bin_input[INPUTS];
  static int output[OUTPUTS];
  static int test_output[OUTPUTS];

  int total = 0;
  int correct = 0;
  for (i = 0; i < max_val; i++) {
    for (j = 0; j < n->num_inputs; j++) {
      bin_input[n->num_inputs - j - 1] = (i >> j) & 1;
    }
    if (eval_network(output, n, bin_input)) {
      fn(test_output, bin_input);
      for (j = 0; j < n->num_outputs; j++) {
        if (output[j] == test_output[j]) {
          correct++;
        }
        total++;
      }
    } else {
      total += n->num_outputs;
    }
  }
  // free(bin_input);
  // free(output);

  return (double)correct / total;
}

void assertTrue(const char* test, int expr) {
  if (expr) {
    printf("%s: PASSED\n", test);
  } else {
    printf("%s: FAILED\n", test);
  }
}

void assertListEq(const char* test, int* l1, int* l2, int length) {
  int i;
  for (i = 0; i < length; i++) {
    if (l1[i] != l2[i]) {
      printf("%s: FAILED\n", test);
      return;
    }
  }
  printf("%s: PASSED\n", test);
}

void assertMatrixEq(const char* test, int** l1, int** l2, int length, int sublength) {
  int i, j;
  for (i = 0; i < length; i++) {
    for (j = 0; j < sublength; j++) {
      if (l1[i][j] == l2[i][j] && l1[i][j] == INDETERMINATE) {
        break;
      }
      if (l1[i][j] != l2[i][j]) {
        printf("%s: FAILED\n", test);
        return;
      }
    }
  }
  printf("%s: PASSED\n", test);
}

void RunTests() {
  int test1_1[4] = {INDETERMINATE, 1, INDETERMINATE, 0};
  int test1_2[4] = {1,1,1,0};
  assertTrue("Subset 1", subset_v(test1_1, test1_2, 4));

  int test2_1[4] = {INDETERMINATE, 1, INDETERMINATE, 0};
  int test2_2[4] = {1,1,1,INDETERMINATE};
  assertTrue("Subset 2", !subset_v(test2_1, test2_2, 4));

  int test3_1[4] = {INDETERMINATE, 1, INDETERMINATE, 0};
  int test3_2[4] = {1,1,1,INDETERMINATE};
  int test3_3[4] = {1,1,1,0};
  int test3_r[4];
  partial_join_v(test3_r, test3_1, test3_2, 4);
  assertListEq("Partial 1", test3_r, test3_3, 4);

  int test4_1[4] = {1,1,1,INDETERMINATE};
  int test4_2[4] = {INDETERMINATE, 1, INDETERMINATE, 0};
  int test4_3[4] = {1,1,1,0};
  int test4_r[4];
  partial_join_v(test4_r, test4_1, test4_2, 4);
  assertListEq("Partial 2", test4_r, test4_3, 4);

  int test5[9] = {0,0,0,0,1,INDETERMINATE,0,INDETERMINATE,INDETERMINATE};
  int test5_r[9];
  gen_truth_table(test5_r, and_g, 2, 1);
  assertListEq("Truth 1", test5_r, test5, 9);

  int test6[9] = {0,1,INDETERMINATE,1,1,1,INDETERMINATE,1,INDETERMINATE};
  int test6_r[9];
  gen_truth_table(test6_r, or_g, 2, 1);
  assertListEq("Truth 2", test6_r, test6, 9);

  int test7[9] = {0,1,INDETERMINATE,1,0,INDETERMINATE,INDETERMINATE,INDETERMINATE,INDETERMINATE};
  int test7_r[9];
  gen_truth_table(test7_r, xor_g, 2, 1);
  assertListEq("Truth 3", test7_r, test7, 9);

  int test8[3] = {1,0,INDETERMINATE};
  int test8_r[3];
  gen_truth_table(test8_r, not_g, 1, 1);
  assertListEq("Truth 4", test8_r, test8, 3);


  int i;

  gate* n1[1];
  n1[0] = (gate*)malloc(sizeof(gate));
  make_gate(n1[0], and_g, 2);
  gate* n_i1[2];
  n_i1[0] = (gate*)malloc(sizeof(gate));
  n_i1[1] = (gate*)malloc(sizeof(gate));
  make_gate(n_i1[0], input_g, 1);
  make_gate(n_i1[1], input_g, 1);
  connect(n_i1[0], n1[0]);
  connect(n_i1[1], n1[0]);
  network n_1 = {n1, 1, n_i1, 2, n1, 1};
  int* output_1[4];
  int* t1[4];
  int t_1_1[1] = {0},
      t_1_2[1] = {0},
      t_1_3[1] = {0},
      t_1_4[1] = {1};
  t1[0] = t_1_1;
  t1[1] = t_1_2;
  t1[2] = t_1_3;
  t1[3] = t_1_4;
  eval_network_all(output_1, &n_1);
  assertMatrixEq("Network 1", output_1, t1, 4, 1);
  int cyclic;
  int deg = degree(&cyclic, &n_1);
  printf("%d\n", cyclic);
  assert(!cyclic);
  free(n1[0]);
  free(n_i1[0]);
  free(n_i1[1]);
  for (i = 0; i < 4; i++) {
    free(output_1[i]);
  }

  gate* n2[2];
  n2[0] = (gate*)malloc(sizeof(gate));
  n2[1] = (gate*)malloc(sizeof(gate));
  make_gate(n2[0], and_g, 2);
  make_gate(n2[1], or_g, 2);
  gate* n_i2[2];
  n_i2[0] = (gate*)malloc(sizeof(gate));
  n_i2[1] = (gate*)malloc(sizeof(gate));
  make_gate(n_i2[0], input_g, 1);
  make_gate(n_i2[1], input_g, 1);
  connect(n_i2[0], n2[0]);
  connect(n_i2[1], n2[1]);
  connect(n2[0], n2[1]);
  connect(n2[1], n2[0]);
  gate* o2[1];
  o2[0] = n2[1];
  network n_2 = {n2, 2, n_i2, 2, o2, 1};
  int* output_2[4];
  int* t2[4];
  int t_2_1[1] = {0},
      t_2_2[1] = {1},
      t_2_3[1] = {INDETERMINATE},
      t_2_4[1] = {1};
  t2[0] = t_2_1;
  t2[1] = t_2_2;
  t2[2] = t_2_3;
  t2[3] = t_2_4;
  eval_network_all(output_2, &n_2);
  assertMatrixEq("Network 2", output_2, t2, 4, 1);
  deg = degree(&cyclic, &n_1);
  assert(cyclic);
  free(n2[0]);
  free(n2[1]);
  free(n_i2[0]);
  free(n_i2[1]);
  for (i = 0; i < 4; i++) {
    free(output_2[i]);
  }

  gate* n3[6];
  n3[0] = (gate*)malloc(sizeof(gate));
  n3[1] = (gate*)malloc(sizeof(gate));
  n3[2] = (gate*)malloc(sizeof(gate));
  n3[3] = (gate*)malloc(sizeof(gate));
  n3[4] = (gate*)malloc(sizeof(gate));
  n3[5] = (gate*)malloc(sizeof(gate));
  make_gate(n3[0], and_g, 2);
  make_gate(n3[1], or_g, 2);
  make_gate(n3[2], and_g, 2);
  make_gate(n3[3], or_g, 2);
  make_gate(n3[4], and_g, 2);
  make_gate(n3[5], or_g, 2);
  gate* n_i3[3];
  n_i3[0] = (gate*)malloc(sizeof(gate));
  n_i3[1] = (gate*)malloc(sizeof(gate));
  n_i3[2] = (gate*)malloc(sizeof(gate));
  make_gate(n_i3[0], input_g, 1);
  make_gate(n_i3[1], input_g, 1);
  make_gate(n_i3[2], input_g, 1);
  connect(n_i3[0], n3[0]);
  connect(n_i3[0], n3[3]);
  connect(n_i3[1], n3[1]);
  connect(n_i3[1], n3[4]);
  connect(n_i3[2], n3[2]);
  connect(n_i3[2], n3[5]);
  connect(n3[0], n3[1]);
  connect(n3[1], n3[2]);
  connect(n3[2], n3[3]);
  connect(n3[3], n3[4]);
  connect(n3[4], n3[5]);
  connect(n3[5], n3[0]);
  network n_3 = {n3, 6, n_i3, 3, n3, 6};
  int* output_3[8];
  int* t3[8];
  int t_3_1[6] = {0, 0, 0, 0, 0, 0},
      t_3_2[6] = {0, 0, 0, 0, 0, 1},
      t_3_3[6] = {0, 1, 0, 0, 0, 0},
      t_3_4[6] = {0, 1, 1, 1, 1, 1},
      t_3_5[6] = {0, 0, 0, 1, 0, 0},
      t_3_6[6] = {1, 1, 1, 1, 0, 1},
      t_3_7[6] = {1, 1, 0, 1, 1, 1},
      t_3_8[6] = {1, 1, 1, 1, 1, 1};
  t3[0] = t_3_1;
  t3[1] = t_3_2;
  t3[2] = t_3_3;
  t3[3] = t_3_4;
  t3[4] = t_3_5;
  t3[5] = t_3_6;
  t3[6] = t_3_7;
  t3[7] = t_3_8;
  eval_network_all(output_3, &n_3);
  assertMatrixEq("Network 3", output_3, t3, 8, 6);
  deg = degree(&cyclic, &n_1);
  assert(cyclic);
  for (i = 0; i < 6; i++) {
    free(n3[i]);
  }
  for (i = 0; i < 3; i++) {
    free(n_i3[i]);
  }
  for (i = 0; i < 8; i++) {
    free(output_3[i]);
  }

  gate* n_4[10];
  for (i = 0; i < 10; i++) {
    n_4[i] = (gate*)malloc(sizeof(gate));
    make_gate(n_4[i], nand_g, 2);
  }
  gate* n_i_4[4];
  for (i = 0; i < 4; i++) {
    n_i_4[i] = (gate*)malloc(sizeof(gate));
    make_gate(n_i_4[i], input_g, 1);
  }
  connect(n_i_4[0], n_4[0]);
  connect(n_i_4[0], n_4[3]);
  connect(n_i_4[1], n_4[0]);
  connect(n_i_4[1], n_4[4]);
  connect(n_i_4[2], n_4[2]);
  connect(n_i_4[2], n_4[1]);
  connect(n_i_4[3], n_4[1]);
  connect(n_i_4[3], n_4[5]);
  connect(n_4[0], n_4[2]);
  connect(n_4[0], n_4[5]);
  connect(n_4[1], n_4[3]);
  connect(n_4[1], n_4[4]);
  connect(n_4[2], n_4[6]);
  connect(n_4[3], n_4[7]);
  connect(n_4[4], n_4[7]);
  connect(n_4[5], n_4[6]);
  connect(n_4[6], n_4[8]);
  connect(n_4[7], n_4[8]);
  connect(n_4[8], n_4[9]);
  connect(n_4[8], n_4[9]);

  gate* output[1];
  output[0] = n_4[9];

  network net = {n_4, 10, n_i_4, 4, output, 1};

  assertTrue("Goal 1", eval_network_fitness(&net, goal1) == 1.0);
  deg = degree(&cyclic, &net);
  assert(!cyclic);
  // printf("Degree: %d\n", degree(&net));
  for (i = 0; i < 10; i++) {
    free(n_4[i]);
  }
  for (i = 0; i < 4; i++) {
    free(n_i_4[i]);
  }
}

typedef struct {
  int* DNA;
  int DNA_length;
  network* network;
  double fitness;
} circuit;

uint32_t rand_range(sfmt_t* sfmt, uint32_t min, uint32_t max)
{
    int r;
    const uint32_t range = max - min;
    const uint32_t buckets = UINT32_MAX / range;
    const uint32_t limit = buckets * range;

    do
    {
        r = sfmt_genrand_uint32(sfmt);
    } while (r >= limit);

    return min + (r / buckets);
}

void random_dna(sfmt_t* sfmt, circuit* c) {
  int i;
  for (i = 0; i < c->DNA_length; i++) {
    c->DNA[i] = rand_range(sfmt, 0, GATES + INPUTS);
  }
}

void mutate(sfmt_t* sfmt, circuit* c) {
  int i;
  if (sfmt_genrand_real1(sfmt) < MUTATION) {
    int mutation_gate = rand_range(sfmt, 0, DNA_LENGTH);
    c->DNA[mutation_gate] = rand_range(sfmt, 0, GATES + INPUTS);
    for (i = 0; i < c->DNA_length; i++) {
      if (sfmt_genrand_real1(sfmt) < 0.001) {
        c->DNA[i] = rand_range(sfmt, 0, GATES + INPUTS);
      }
    }
  }
}

// void mutate(sfmt_t* sfmt, circuit* c) {
//   int i;
  // for (i = 0; i < c->DNA_length; i++) {
  //   if (sfmt_genrand_real1(sfmt) < 0.2) {
  //     c->DNA[i] = sfmt_genrand_uint32(sfmt) & 0xFF;
  //   }
  // }
// }

void create_circuit_network(circuit* c) {
  c->network->num_gates = GATES;
  c->network->num_inputs = INPUTS;
  c->network->num_outputs = OUTPUTS;

  gate** gates = (gate**)malloc(sizeof(gate*) * c->network->num_gates);
  gate** inputs = (gate**)malloc(sizeof(gate*) * c->network->num_inputs);
  gate** output = (gate**)malloc(sizeof(gate*) * c->network->num_outputs);

  int i;
  for (i = 0; i < c->network->num_gates; i++) {
    gates[i] = (gate*)malloc(sizeof(gate));
    make_gate(gates[i], nand_g, 2);
  }
  for (i = 0; i < c->network->num_inputs; i++) {
    inputs[i] = (gate*)malloc(sizeof(gate));
    make_gate(inputs[i], input_g, 1);
  }

  c->network->gates = gates;
  c->network->inputs = inputs;
  c->network->output = output;
}

void circuitize(circuit* c) {
  int i, j;
  for (i = 0; i < c->network->num_gates; i++) {
    reset_gate(c->network->gates[i]);
  }
  for (i = 0; i < c->network->num_inputs; i++) {
    reset_gate(c->network->inputs[i]);
  }

  int dna_pos = 0;
  for (i = 0; i < GATES; i++) {
    for (j = 0; j < INPUTS_PER_GATE; j++) {
      uint16_t address = c->DNA[dna_pos++];
      if (address >= INPUTS) {
        connect(c->network->gates[address - INPUTS], c->network->gates[i]);
      } else {
        connect(c->network->inputs[address], c->network->gates[i]);
      }
    }
  }

  // test assert code: TO DELETE
  for (i = 0; i < c->network->num_gates; i++) {
    assert(c->network->gates[i]->num_inputs == 2);
  }

  for (i = 0; i < OUTPUTS; i++) {
    uint16_t address = c->DNA[dna_pos++];
    if (address >= INPUTS) {
      c->network->output[i] = c->network->gates[address - INPUTS];
    } else {
      c->network->output[i] = c->network->gates[address];
    }
  }
}

void make_circuit(circuit* c) {
  c->network = (network*)malloc(sizeof(network));
  c->DNA_length = DNA_LENGTH;
  c->DNA = (int*)malloc(sizeof(int) * c->DNA_length);
}

void free_gate(gate* g) {
  free(g->inputs);
  free(g->outputs);
  free(g->cache);
}

void free_network(network* n) {
  int i;
  for (i = 0; i < n->num_gates; i++) {
    free_gate(n->gates[i]);
    free(n->gates[i]);
  }
  for (i = 0; i < n->num_inputs; i++) {
    free_gate(n->inputs[i]);
    free(n->inputs[i]);
  }
  free(n->gates);
  free(n->inputs);
  free(n->output);
}

int circuit_compare(const void* c1, const void* c2) {
  if (((circuit*)c1)->fitness < ((circuit*)c2)->fitness) {
    return -1;
  } else if (((circuit*)c1)->fitness == ((circuit*)c2)->fitness) {
    return 0;
  } else {
    return 1;
  }
}

int time_to_perfect(sfmt_t* sfmt) {
  circuit circuits[CIRCUITS];

  int i, j;
  for (i = 0; i < CIRCUITS; i++) {
    make_circuit(&circuits[i]);
    random_dna(sfmt, &circuits[i]);
    create_circuit_network(&circuits[i]);
  }

  double max_fitness = 0.0;

  int max_degree = -1;
  int max_cyclic = -1;

  int reached = -1;
  for (j = 0; ; j++) {
    for (i = 0; i < CIRCUITS; i++) {
      circuitize(&circuits[i]);
      circuits[i].fitness = eval_network_fitness_vector(circuits[i].network, rivest);

      if (circuits[i].fitness > max_fitness) {
        int cyclic;
        int deg = degree(&cyclic, circuits[i].network);
        max_fitness = circuits[i].fitness;
        max_degree = deg;
        max_cyclic = cyclic;
      } else if (circuits[i].fitness == max_fitness) {
        int cyclic;
        int deg = degree(&cyclic, circuits[i].network);
        if (max_cyclic && !cyclic) {
          max_cyclic = cyclic;
        }
        if (deg < max_degree) {
          max_degree = deg;
          max_cyclic = cyclic;
        }
      }
    }
    qsort(circuits, CIRCUITS, sizeof(circuit), circuit_compare);
    for (i = 0; i < CIRCUITS; i++) {
      if (i < ELITE) {
        memcpy(circuits[i].DNA, circuits[CIRCUITS - i - 1].DNA, circuits[i].DNA_length);
        mutate(sfmt, &circuits[i]);
      } else if (i < CIRCUITS - ELITE - 1) {
        mutate(sfmt, &circuits[i]);
      }
    }

    if (reached == -1 && max_fitness == 1.0) {
      reached = j;
      break;
    }
    if ((j + 1) % 100 == 0) {
      printf("%d: %f %f (deg: %d; cyclic: %d)\n", j + 1, circuits[CIRCUITS-1].fitness, max_fitness, max_degree, max_cyclic);
    }
  }

  for (i = 0; i < CIRCUITS; i++) {
    free(circuits[i].DNA);
    free_network(circuits[i].network);
    free(circuits[i].network);
  }

  return reached;
}

int main(int argc, char** argv) {
  RunTests();

  sfmt_t sfmt;
  sfmt_init_gen_rand(&sfmt, SEED);

  int i;
  int total = 0;
  for (i = 0; i < 50; i++) {
    int reached = time_to_perfect(&sfmt);
    total += reached;
    printf("%d: %d (%0.2f)\n", i+1, reached, (double)total / (i+1));
  }

  return 0;
}