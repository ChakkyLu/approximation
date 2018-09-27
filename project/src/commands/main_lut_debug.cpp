#include "global.h"
#include <bitset>
#include <ctime>

using namespace std;

// ABC headers
#include "base/abc/abc.h"
#include "base/main/main.h"
#include "base/main/mainInt.h"

// Glucose headers
#include "simp/SimpSolver.h"
#include "core/SolverTypes.h"
#include "core/Dimacs.h"
#include "mtl/Vec.h"
#include "utils/System.h"

#include "circuit.h"

using namespace nodecircuit;

#include "cnf.h"
#include "glucose_interface.h"

using namespace satcnf;

class LutDebugStep {
public:
  LutDebugStep(long insize, long outsize, long lutsize) {
    inputs.resize(insize);
    outputs.resize(outsize);
    lut_vars.resize(lutsize);
    reduced_in_impl_cir = NULL;
    removed_lut_impl_cir = NULL;
    reduced_in_impl_cnf = NULL;
    removed_lut_impl_cnf = NULL;
    offset = 0;
  };

  virtual ~LutDebugStep() {
    // TODO: delete everything!
  };

public:
  ValVector inputs;
  ValVector outputs;
  ValVector lut_vars;

  Circuit *reduced_in_impl_cir;
  Circuit *removed_lut_impl_cir;

  Cnf *reduced_in_impl_cnf;
  Cnf *removed_lut_impl_cnf;

  long offset;
};


int LutDebug(Circuit *spec_cir, Circuit *impl_cir, ValVector& lut_result, vector<ValVector> &init_inputs) {
  long i, j, k;
  vector<LutDebugStep *> debug_steps;
  LutDebugStep *cur_debug_step = NULL;
  LutDebugStep *prev_debug_step = NULL;

  lut_result.clear();

  int cur_pattern_num = 0;

  Cnf org_impl_cnf;
  org_impl_cnf.Convert2Cnf(impl_cir);

  long input_size = spec_cir->inputs.size();
  long output_size = spec_cir->outputs.size();
  long lut_size = impl_cir->inputs.size() - input_size;


  cur_debug_step = new LutDebugStep(input_size, output_size, lut_size);

  cout << "INPUTS:" << endl;
  cur_debug_step->inputs = init_inputs[cur_pattern_num];
  for (i = 0; i < input_size; i++) {
    cout << cur_debug_step->inputs[i];
  }
  cout << endl;
  spec_cir->Simulate(cur_debug_step->inputs, cur_debug_step->outputs);
  cout << "OUTPUTS:" << endl;
  for (i = 0; i < output_size; i++) {
    cout << cur_debug_step->outputs[i];
  }
  cout << endl;
  cur_pattern_num++;

  GInterface sat_candid_gen;

  while (cur_pattern_num < init_inputs.size()) {
    Circuit *dup_impl = impl_cir->GetDuplicate("", "", "");

    dup_impl->ApplyInOutSimplify(cur_debug_step->inputs, cur_debug_step->outputs);
    dup_impl->RemoveBufNot();
    //dup_impl->WriteBlif("dup-impl-reduced.blif");

    cur_debug_step->reduced_in_impl_cir = dup_impl;

    cur_debug_step->reduced_in_impl_cnf = new Cnf;
    cur_debug_step->reduced_in_impl_cnf->Convert2Cnf(dup_impl);

    // add the output constraints.
    // input constraints are already applied to simplify the circuit!
    for (j = 0; j < output_size; j++) {
      LiteralType lit = cur_debug_step->reduced_in_impl_cnf->literal_node_map[dup_impl->outputs[j]];
      LiteralVector clause;
      clause.push_back(cur_debug_step->outputs[j] ? lit : -lit);
      cur_debug_step->reduced_in_impl_cnf->clauses.push_back(clause);
    }

    sat_candid_gen.AddCnf(cur_debug_step->reduced_in_impl_cnf->clauses, cur_debug_step->offset);

    // TODO: add in/out equalities
    if (prev_debug_step) {
      for (k = 0; k < lut_size; k++)
        sat_candid_gen.SetEqual(
            cur_debug_step->reduced_in_impl_cnf->literal_node_map[dup_impl->inputs[k + input_size]] +
            cur_debug_step->offset,
            prev_debug_step->reduced_in_impl_cnf->literal_node_map[prev_debug_step->reduced_in_impl_cir->inputs[k + input_size]] +
            prev_debug_step->offset);
    }

    long new_offset = cur_debug_step->offset + cur_debug_step->reduced_in_impl_cnf->literal_str_map.size();
    long next_offset = new_offset;

    if (prev_debug_step) {
      delete prev_debug_step->reduced_in_impl_cir;
      delete prev_debug_step->reduced_in_impl_cnf;
      delete prev_debug_step->removed_lut_impl_cir;
      delete prev_debug_step->removed_lut_impl_cnf;
      prev_debug_step->reduced_in_impl_cir = NULL;
      prev_debug_step->reduced_in_impl_cnf = NULL;
      prev_debug_step->removed_lut_impl_cir = NULL;
      prev_debug_step->removed_lut_impl_cnf = NULL;
    }

    //debug_steps.push_back(cur_debug_step);
    prev_debug_step = cur_debug_step;
    cur_debug_step = new LutDebugStep(input_size, output_size, lut_size);
    cur_debug_step->offset = new_offset;

    prev_debug_step->removed_lut_impl_cir = new Circuit;
    prev_debug_step->removed_lut_impl_cnf = new Cnf;

    cout << "INPUTS:" << endl;
    cur_debug_step->inputs = init_inputs[cur_pattern_num];
    for (i = 0; i < input_size; i++) {
      cout << cur_debug_step->inputs[i];
    }
    cout << endl;
    spec_cir->Simulate(cur_debug_step->inputs, cur_debug_step->outputs);
    cout << "OUTPUTS:" << endl;
    for (i = 0; i < output_size; i++) {
      cout << cur_debug_step->outputs[i];
    }
    cout << endl;

    cur_pattern_num++;

  };

  do {
    cout << endl << "-----> " << debug_steps.size() << " <-----" << endl;
    Circuit *dup_impl = impl_cir->GetDuplicate("", "", "");

    dup_impl->ApplyInOutSimplify(cur_debug_step->inputs, cur_debug_step->outputs);
    dup_impl->RemoveBufNot();
    dup_impl->WriteBlif("dup-impl-reduced.blif");

    cur_debug_step->reduced_in_impl_cir = dup_impl;

    cur_debug_step->reduced_in_impl_cnf = new Cnf;
    cur_debug_step->reduced_in_impl_cnf->Convert2Cnf(dup_impl);

    // add the output constraints.
    // input constraints are already applied to simplify the circuit!
    for (j = 0; j < output_size; j++) {
      LiteralType lit = cur_debug_step->reduced_in_impl_cnf->literal_node_map[dup_impl->outputs[j]];
      LiteralVector clause;
      clause.push_back(cur_debug_step->outputs[j] ? lit : -lit);
      cur_debug_step->reduced_in_impl_cnf->clauses.push_back(clause);
    }

    sat_candid_gen.AddCnf(cur_debug_step->reduced_in_impl_cnf->clauses, cur_debug_step->offset);

    // TODO: add in/out equalities
    if (prev_debug_step) {
      for (k = 0; k < lut_size; k++)
        sat_candid_gen.SetEqual(
            cur_debug_step->reduced_in_impl_cnf->literal_node_map[dup_impl->inputs[k + input_size]] +
            cur_debug_step->offset,
            prev_debug_step->reduced_in_impl_cnf->literal_node_map[prev_debug_step->reduced_in_impl_cir->inputs[k + input_size]] +
            prev_debug_step->offset);
    }

    bool sat_res1 = sat_candid_gen.GetAnswer();
    if (!sat_res1) {
      MSG("No Answer!");
      return -1;
    }
    cout << "LUT Vars: " << endl;
    for (k = 0; k < lut_size; k++) {
      cur_debug_step->lut_vars[k] = sat_candid_gen.ModelNum(cur_debug_step->reduced_in_impl_cnf->literal_node_map[dup_impl->inputs[k + input_size]]) == l__True;
      cout << cur_debug_step->lut_vars[k];
    }
    cout << endl;


    dup_impl = impl_cir->GetDuplicate("", "", "");
    NodeVector shuffled_in;
    for (j = 0; j < lut_size; j++)
      shuffled_in.push_back(dup_impl->inputs[input_size + j]);
    for (k = 0; k < input_size; k++)
      shuffled_in.push_back(dup_impl->inputs[k]);
    dup_impl->inputs.clear();
    dup_impl->inputs = shuffled_in;

    dup_impl->ApplyInOutSimplify(cur_debug_step->lut_vars, cur_debug_step->outputs);
    dup_impl->RemoveBufNot();
    dup_impl->WriteBlif("dup-impl-removed.blif");

    cur_debug_step->removed_lut_impl_cir = dup_impl;

    cur_debug_step->removed_lut_impl_cnf = new Cnf;
    cur_debug_step->removed_lut_impl_cnf->Convert2Cnf(dup_impl);

    long new_offset = cur_debug_step->offset + cur_debug_step->reduced_in_impl_cnf->literal_str_map.size();
    long next_offset = new_offset + cur_debug_step->removed_lut_impl_cnf->literal_str_map.size();

    GInterface sat_next_candid_input_gen(sat_candid_gen);
    sat_next_candid_input_gen.AddCnf(cur_debug_step->removed_lut_impl_cnf->clauses, new_offset);
    sat_next_candid_input_gen.AddCnf(org_impl_cnf.clauses, next_offset);
    // TODO: add in/out equalities
    for (k = 0; k < input_size; k++)
      sat_next_candid_input_gen.SetEqual(
          cur_debug_step->removed_lut_impl_cnf->literal_node_map[cur_debug_step->removed_lut_impl_cir->inputs[k]] +
          new_offset,
          org_impl_cnf.literal_node_map[impl_cir->inputs[k]] + next_offset);
    for (k = 0; k < lut_size; k++)
      sat_next_candid_input_gen.SetEqual(
          cur_debug_step->reduced_in_impl_cnf->literal_node_map[cur_debug_step->reduced_in_impl_cir->inputs[k]] +
          cur_debug_step->offset,
          org_impl_cnf.literal_node_map[impl_cir->inputs[k + input_size]] + next_offset);
    LiteralVector lit_v1, lit_v2;
    for (k = 0; k < output_size; k++) {
      lit_v1.push_back(
          cur_debug_step->removed_lut_impl_cnf->literal_node_map[cur_debug_step->removed_lut_impl_cir->outputs[k]] +
          new_offset);
      lit_v2.push_back(org_impl_cnf.literal_node_map[impl_cir->outputs[k]] + next_offset);
    }
    sat_next_candid_input_gen.SetNEqual(lit_v1, lit_v2);
    // exclude previous inputs
    LiteralVector all_exclude_literals;
    for (k = 0; k < input_size; k++) {
      LiteralType lit =
          cur_debug_step->removed_lut_impl_cnf->literal_node_map[cur_debug_step->removed_lut_impl_cir->inputs[k]] +
          new_offset;
      all_exclude_literals.push_back(lit);
    }
    ClauseVector exclude_clauses;
    LiteralVector exclude_clause;
    exclude_clause.resize(input_size);
    for (j = 0; j < debug_steps.size(); j++) {
      for (i = 0; i < input_size; i++) {
        exclude_clause[i] = debug_steps[j]->inputs[i] ? -all_exclude_literals[i] : all_exclude_literals[i];
      }
      exclude_clauses.push_back(exclude_clause);
    }
    for (i = 0; i < input_size; i++) {
      exclude_clause[i] = cur_debug_step->inputs[i] ? -all_exclude_literals[i] : all_exclude_literals[i];
    }
    exclude_clauses.push_back(exclude_clause);
    sat_next_candid_input_gen.AddCnf(exclude_clauses, 0);

    bool sat_res2 = sat_next_candid_input_gen.GetAnswer();
    if (!sat_res2) {
      MSG("Solution found!");
      lut_result.resize(lut_size);
      lut_result = cur_debug_step->lut_vars;
      return -2;
    }
    cout << "ANOTHER LUT Vars: " << endl;
    ValVector new_lut_vars;
    new_lut_vars.resize(lut_size);
    for (k = 0; k < lut_size; k++) {
      new_lut_vars[k] = sat_next_candid_input_gen.ModelNum(
          cur_debug_step->reduced_in_impl_cnf->literal_node_map[cur_debug_step->reduced_in_impl_cir->inputs[k + input_size]] +
          cur_debug_step->offset) == l__True;
      cout << new_lut_vars[k];
    }
    cout << endl;

    if (prev_debug_step) {
      delete prev_debug_step->reduced_in_impl_cir;
      delete prev_debug_step->reduced_in_impl_cnf;
      delete prev_debug_step->removed_lut_impl_cir;
      delete prev_debug_step->removed_lut_impl_cnf;
      prev_debug_step->reduced_in_impl_cir = NULL;
      prev_debug_step->reduced_in_impl_cnf = NULL;
      prev_debug_step->removed_lut_impl_cir = NULL;
      prev_debug_step->removed_lut_impl_cnf = NULL;
    }

    debug_steps.push_back(cur_debug_step);
    prev_debug_step = cur_debug_step;
    cur_debug_step = new LutDebugStep(input_size, output_size, lut_size);
    cur_debug_step->offset = new_offset;

    cout << "NEW INPUT: " << endl;
    for (k = 0; k < input_size; k++) {
      cur_debug_step->inputs[k] = sat_next_candid_input_gen.ModelNum(prev_debug_step->removed_lut_impl_cnf->literal_node_map[prev_debug_step->removed_lut_impl_cir->inputs[k]] + new_offset) == l__True;
      cout << cur_debug_step->inputs[k];
    }
    cout << endl;

    spec_cir->Simulate(cur_debug_step->inputs, cur_debug_step->outputs);
    cout << "OUTPUTS:" << endl;
    for (i = 0; i < output_size; i++) {
      cout << cur_debug_step->outputs[i];
    }
    cout << endl;

  } while (1);

  return 0;
}