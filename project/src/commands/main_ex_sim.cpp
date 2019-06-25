#include "global.h"
#include <cmath>

#include <bitset>

#include <iomanip>

#include <thread>

#include <fstream>

#include <sstream>

#include <algorithm>

#include <string>

using namespace std;

#include "circuit.h"

using namespace nodecircuit;

#include "neural_network.h"

using namespace neuralnetwork;

const uint64_t m1 = 0x5555555555555555;
const uint64_t m2 = 0x3333333333333333;
const uint64_t m4 = 0x0f0f0f0f0f0f0f0f;
const uint64_t m8 = 0x00ff00ff00ff00ff;
const uint64_t m16 = 0x0000ffff0000ffff;
const uint64_t m32 = 0x00000000ffffffff;
const uint64_t hff = 0xffffffffffffffff;
const uint64_t h01 = 0x0101010101010101;

inline uint64_t Count64_opt(uint64_t x) {
  x -= (x >> 1) & m1;
  x = (x & m2) + ((x >> 2) & m2);
  x = (x + (x >> 4)) & m4;
  return (x * h01) >> 56;
}

inline uint64_t Count64_sh(uint64_t data) {
  uint64_t cnt = 0;
  for (int i = 0; data != 0 && i < 64; i++) {
    if (data%2)
      cnt++;
    data >>= 1;
  }
  return cnt;
}

inline bit64 Count64(bit64 data) {
  if (data == ONE64)
    return 64;
  if (data == ZERO64)
    return 0;
  return Count64_sh(data);
}


int DoExSim(string spec_filename, string output_filename) {
  Circuit new_spec_cir;

  new_spec_cir.ReadBlif(spec_filename);

  if (new_spec_cir.inputs.size() == 0) {
    MSG("Cannot simulate circuit with no primary input!");
    return 0;
  }

  if (new_spec_cir.inputs.size() > 64) {
    MSG("Currently only simulating up to 64 inputs!");
    return 0;
  }

  ofstream resfile(output_filename);
  if (!resfile.good() || !resfile.is_open()) {
    ERR("cannot open simulation result file: " + output_filename);
    return 0;
  }

  auto start_time = chrono::system_clock::now();

  Val64Vector pi_vect;
  Val64Vector po_vect;
  pi_vect.resize(new_spec_cir.inputs.size());
  po_vect.resize(new_spec_cir.outputs.size());
  int pi_index = 0;
  while (pi_index < new_spec_cir.inputs.size() && pi_index < NUM_SIG_PAR_CHANGE) {
    pi_vect[pi_index] = PAT64[pi_index];
    pi_index++;
  }
  while (pi_index < new_spec_cir.inputs.size()) {
    pi_vect[pi_index] = ZERO64;
    pi_index++;
  }

  if (pi_vect.size() <= NUM_SIG_PAR_CHANGE) {
    new_spec_cir.Simulate(pi_vect, po_vect);
    int num_patterns = pow(2, pi_vect.size());
    vector<string> res_strs;
    res_strs.resize(num_patterns);
    for (int i = 0; i < num_patterns; i++)
      res_strs[i].resize(po_vect.size(), '0');
    for (int i = 0; i < po_vect.size(); i++) {
      bit64 cur_res = po_vect[i];
      for (int j = 0; j < num_patterns; j++) {
        if (cur_res%2)
          res_strs[j][i] = '1';
        cur_res >>= 1;
      }
    }
    for (int k = 0; k < res_strs.size(); k++)
      resfile << res_strs[k] << endl;
  }
  else {
    bit64 num_sim_run = pow(2, new_spec_cir.inputs.size()-NUM_SIG_PAR_CHANGE);
    bit64 prev_sim_run = 0;
    bit64 num_sim_run_scale = num_sim_run >> 3; // 8 scales
    char progress_str[9] = {'.','.','.','.','.','.','.','.',0};
    for (bit64 cur_sim_run = 0; cur_sim_run < num_sim_run; ) {
      if (num_sim_run_scale > 1 && cur_sim_run%num_sim_run_scale == 0) {
        MSG(progress_str);
        progress_str[cur_sim_run/num_sim_run_scale] = '*';
      }
      new_spec_cir.Simulate(pi_vect, po_vect);
      vector<string> res_strs;
      res_strs.resize(64);
      for (int i = 0; i < 64; i++)
        res_strs[i].resize(po_vect.size(), '0');
      for (int i = 0; i < po_vect.size(); i++) {
        bit64 cur_res = po_vect[i];
        for (int j = 0; j < 64; j++) {
          if (cur_res%2)
            res_strs[j][i] = '1';
          cur_res >>= 1;
        }
      }
      for (int k = 0; k < res_strs.size(); k++)
        resfile << res_strs[k] << endl;

      prev_sim_run = cur_sim_run;
      cur_sim_run++;
      if (cur_sim_run < num_sim_run) {
        bit64 diff_sim_run = prev_sim_run ^cur_sim_run;
        bit64 cur_sim_run_dup = cur_sim_run;
        bit64 pi_vect_index = NUM_SIG_PAR_CHANGE;
        while (diff_sim_run != 0) {
          if (diff_sim_run%2) {
            if (cur_sim_run_dup%2) {
              pi_vect[pi_vect_index] = ONE64;
            }
            else {
              pi_vect[pi_vect_index] = ZERO64;
            }
          }
          diff_sim_run >>= 1;
          cur_sim_run_dup >>= 1;
          pi_vect_index++;
        }
      }
    }

  }

  auto end_time = chrono::system_clock::now();
chrono::duration<double> diff_time = end_time - start_time;
  cout << endl << "Simulation time is: " << setprecision(3) << diff_time.count() << " s" << endl;

  resfile.close();

  return 0;
}
int DoSimEq(Circuit &new_spec_cir, Circuit &new_impl_cir, vector<bit64> &po_diff_count, vector<bit64> &total_diff_count) {

  auto start_time = chrono::system_clock::now();

  Val64Vector pi_vect;
  Val64Vector po_vect_spec;
  Val64Vector po_vect_impl;
  pi_vect.resize(new_spec_cir.inputs.size());
  po_vect_spec.resize(new_spec_cir.outputs.size());
  po_vect_impl.resize(new_impl_cir.outputs.size());
  int pi_index = 0;
  while (pi_index < new_spec_cir.inputs.size() && pi_index < NUM_SIG_PAR_CHANGE) {
    pi_vect[pi_index] = PAT64[pi_index];
    pi_index++;
  }
  while (pi_index < new_spec_cir.inputs.size()) {
    pi_vect[pi_index] = ZERO64;
    pi_index++;
  }

  total_diff_count.clear();
  total_diff_count.resize(8, 0);
  po_diff_count.clear();
  po_diff_count.resize(po_vect_impl.size() ,0);

  if (pi_vect.size() <= NUM_SIG_PAR_CHANGE) {
    new_spec_cir.Simulate(pi_vect, po_vect_spec);
    new_impl_cir.Simulate(pi_vect, po_vect_impl);
    bit64 diff_data = 0;
    for (int i = 0; i < po_vect_spec.size(); i++) {
      diff_data |= (po_vect_spec[i] ^ po_vect_impl[i]);
    }
    int num_patterns = pow(2, pi_vect.size());
    int num_diff = 0;
    for (int i = 0; diff_data != 0 && i < num_patterns; i++) {
      if (diff_data%2)
        num_diff++;
      diff_data >>= 1;
    }
    // TODO: collect the difference count for each po
    for (int i = 0; i < 8; i++)
      total_diff_count[i] = (i+1)*(num_diff/8);
    total_diff_count[7] = num_diff;
    cout << "Number of different minterms: " << num_diff << endl;
  }
  else {
    bit64 num_sim_run = pow(2, new_spec_cir.inputs.size()-NUM_SIG_PAR_CHANGE);
    bit64 prev_sim_run = 0;
    bit64 num_diff = 0;
    bit64 num_sim_run_scale = num_sim_run >> 3; // 8 scales
    char progress_str[9] = {'.','.','.','.','.','.','.','.',0};
    for (bit64 cur_sim_run = 0; cur_sim_run < num_sim_run; ) {
      if (num_sim_run_scale > 1 && cur_sim_run%num_sim_run_scale == 0) {
        // cout << progress_str << " -> number of minterm differences so far: " << num_diff << endl;
        progress_str[cur_sim_run/num_sim_run_scale] = '*';
        total_diff_count[cur_sim_run/num_sim_run_scale] = num_diff;
      }
      new_spec_cir.Simulate(pi_vect, po_vect_spec);
      new_impl_cir.Simulate(pi_vect, po_vect_impl);
      bit64 diff_data = 0;
      for (int i = 0; i < po_vect_spec.size(); i++) {
        bit64 cur_diff_data = po_vect_spec[i] ^ po_vect_impl[i];
        diff_data |= cur_diff_data;
        po_diff_count[i] += Count64(cur_diff_data);
      }
      num_diff += Count64(diff_data);

      prev_sim_run = cur_sim_run;
      cur_sim_run++;
      if (cur_sim_run < num_sim_run) {
        bit64 diff_sim_run = prev_sim_run ^cur_sim_run;
        bit64 cur_sim_run_dup = cur_sim_run;
        bit64 pi_vect_index = NUM_SIG_PAR_CHANGE;
        while (diff_sim_run != 0) {
          if (diff_sim_run%2) {
            if (cur_sim_run_dup%2) {
              pi_vect[pi_vect_index] = ONE64;
            }
            else {
              pi_vect[pi_vect_index] = ZERO64;
            }
          }
          diff_sim_run >>= 1;
          cur_sim_run_dup >>= 1;
          pi_vect_index++;
        }
      }
    }
    total_diff_count[7] = num_diff;
    cout << "Number of different minterms: " << num_diff << endl;
    cout << "Differences of minterm for each PO:" << endl;
    for (int i = 0; i < po_vect_spec.size(); i++) {
      cout << "Number of different minterms for " << new_spec_cir.outputs[i]->name << " : " << po_diff_count[i] << endl;
    }
  }

  auto end_time = chrono::system_clock::now();
  chrono::duration<double> diff_time = end_time - start_time;
  cout << endl << "Simulation time is: " << setprecision(3) << diff_time.count() << " s" << endl;

  return 0;

}

int DoSimEq(Circuit &new_spec_cir, Circuit &new_impl_cir) {
  vector<bit64> po_diff_count;
  vector<bit64> total_diff_count;
  return DoSimEq(new_spec_cir, new_impl_cir, po_diff_count, total_diff_count);
}

// function to be used for multi-threading
void DoSimEqTh4(Circuit *new_spec_cir, Circuit *new_impl_cir, vector<bit64> *po_diff_count, vector<bit64> *total_diff_count) {
  DoSimEq(*new_spec_cir, *new_impl_cir, *po_diff_count, *total_diff_count);
}

int DoSimEq(string spec_filename, string impl_filename) {
  Circuit new_spec_cir, new_impl_cir;

  new_spec_cir.ReadBlif(spec_filename);
  new_impl_cir.ReadBlif(impl_filename);


  if (new_spec_cir.inputs.size() == 0) {
    MSG("Cannot simulate circuit with no primary input!");
    return 0;
  }

  if (new_spec_cir.inputs.size() > 64) {
    MSG("Currently only simulating up to 64 inputs!");
    return 0;
  }

  if (new_spec_cir.inputs.size() != new_impl_cir.inputs.size() || new_spec_cir.outputs.size() != new_impl_cir.outputs.size()) {
    MSG("Number of Pi/Po of two circuits are different!");
    return 0;
  }

  return DoSimEq(new_spec_cir, new_impl_cir);
}

bit64 MyRand64() {
  bit64 rnd_out = std::rand();
  rnd_out <<= 16;
  rnd_out |= (std::rand() & 0x0FFFF);
  rnd_out <<= 16;
  rnd_out |= (std::rand() & 0x0FFFF);
  rnd_out <<= 16;
  rnd_out |= (std::rand() & 0x0FFFF);
  return rnd_out;
}

int RndSimEq(string spec_filename, string impl_filename, int num_patterns) {
  Circuit new_spec_cir, new_impl_cir;

  new_spec_cir.ReadBlif(spec_filename);
  new_impl_cir.ReadBlif(impl_filename);

  if (new_spec_cir.inputs.size() == 0) {
    MSG("Cannot simulate circuit with no primary input!");
    return 0;
  }

  if (new_spec_cir.inputs.size() != new_impl_cir.inputs.size() || new_spec_cir.outputs.size() != new_impl_cir.outputs.size()) {
    MSG("Number of Pi/Po of two circuits are different!");
    return 0;
  }

  auto start_time = chrono::system_clock::now();

  std::srand(std::time(0));

  Val64Vector pi_vect;
  Val64Vector po_vect_spec;
  Val64Vector po_vect_impl;
  pi_vect.resize(new_spec_cir.inputs.size());
  po_vect_spec.resize(new_spec_cir.outputs.size());
  po_vect_impl.resize(new_impl_cir.outputs.size());

  bit64 num_sim_run = num_patterns/64;
  if (num_patterns%64)
    num_sim_run++;
  bit64 num_diff = 0;
  for (bit64 cur_sim_run = 0; cur_sim_run < num_sim_run; cur_sim_run++) {
    int pi_index = 0;
    while (pi_index < new_spec_cir.inputs.size()) {
      bit64 new_rnd = MyRand64();
      //cout << std::bitset<64>(new_rnd) << endl;
      // TODO: check it is different from previous simulation pattern(s)
      pi_vect[pi_index] = new_rnd;
      pi_index++;
    }
    new_spec_cir.Simulate(pi_vect, po_vect_spec);
    new_impl_cir.Simulate(pi_vect, po_vect_impl);
    bit64 diff_data = 0;
    for (int i = 0; i < po_vect_spec.size(); i++) {
      diff_data |= (po_vect_spec[i] ^ po_vect_impl[i]);
    }
    if (diff_data == ONE64)
      num_diff += 64;
    else
      for (int i = 0; diff_data != 0 && i < 64; i++) {
        if (diff_data % 2)
          num_diff++;
        diff_data >>= 1;
      }

  }
  cout << "Number of different minterms: " << num_diff << " (out of " << num_sim_run*64 << " runs)" << endl;

  auto end_time = chrono::system_clock::now();
  chrono::duration<double> diff_time = end_time - start_time;
  // cout << endl << "Simulation time is: " << setprecision(3) << diff_time.count() << " s" << endl;

  return 0;
}


int DoSimEqMod(string spec_filename, string target_node_name) {
  Circuit new_spec_cir;

  new_spec_cir.ReadBlif(spec_filename);

  if (new_spec_cir.all_nodes_map.find(target_node_name) == new_spec_cir.all_nodes_map.end()) {
    ERR("Node not found: " + target_node_name);
    return -1;
  }
  Node* target_node = new_spec_cir.all_nodes_map[target_node_name];
    NodeVector target_nodes;
    target_nodes.push_back(target_node);


  if (new_spec_cir.inputs.size() == 0) {
    MSG("Cannot simulate circuit with no primary input!");
    return 0;
  }

  if (new_spec_cir.inputs.size() > 64) { // TODO: do we need to limit the inputs (because of high runtime?)
    MSG("Currently only simulating up to 64 inputs!");
    return 0;
  }


  if (target_node->inputs.size() == 2) {
    const int num_changes = 2+(2*2)+5+4;
    Circuit* spec_cirs[num_changes];
    Circuit* impl_cirs[num_changes];
    std::thread sim_threads[num_changes];
    vector<bit64> po_diff_count[num_changes];
    vector<bit64> total_diff_count[num_changes];

    int cur_change_num = 0;
    Node* cur_node = NULL;

    spec_cirs[cur_change_num] = new_spec_cir.GetDuplicate("", "", "");
    impl_cirs[cur_change_num] = new_spec_cir.GetDuplicate("", "", "");
    cur_node = impl_cirs[cur_change_num]->all_nodes_map[target_node->name];
    cur_node->type = NODE_ZERO;
    if (cur_node->inputs[0]->outputs.size() == 1)
      cur_node->inputs[0]->outputs.clear();
    else {
      NodeVector::iterator it = cur_node->inputs[0]->outputs.begin();
      while (it != cur_node->inputs[0]->outputs.end() && *it != cur_node) {
        ++it;
      }
      if (it != cur_node->inputs[0]->outputs.end())
        cur_node->inputs[0]->outputs.erase(it);
    }
    if (cur_node->inputs[1]->outputs.size() == 1)
      cur_node->inputs[1]->outputs.clear();
    else {
      NodeVector::iterator it = cur_node->inputs[1]->outputs.begin();
      while (it != cur_node->inputs[1]->outputs.end() && *it != cur_node) {
        ++it;
      }
      if (it != cur_node->inputs[1]->outputs.end())
        cur_node->inputs[1]->outputs.erase(it);
    }
    cur_node->inputs.clear();
    impl_cirs[cur_change_num]->Simplify();
    impl_cirs[cur_change_num]->SetIndexes();
    cout << "Creating thread " << cur_change_num << " for node type: " << cur_node->type << endl;
    sim_threads[cur_change_num] = std::thread(DoSimEqTh4, spec_cirs[cur_change_num], impl_cirs[cur_change_num], &po_diff_count[cur_change_num], &total_diff_count[cur_change_num]);

    cur_change_num++;

    spec_cirs[cur_change_num] = new_spec_cir.GetDuplicate("", "", "");
    impl_cirs[cur_change_num] = new_spec_cir.GetDuplicate("", "", "");
    cur_node = impl_cirs[cur_change_num]->all_nodes_map[target_node->name];
    cur_node->type = NODE_ONE;
    if (cur_node->inputs[0]->outputs.size() == 1)
      cur_node->inputs[0]->outputs.clear();
    else {
      NodeVector::iterator it = cur_node->inputs[0]->outputs.begin();
      while (it != cur_node->inputs[0]->outputs.end() && *it != cur_node) {
        ++it;
      }
      if (it != cur_node->inputs[0]->outputs.end())
        cur_node->inputs[0]->outputs.erase(it);
    }
    if (cur_node->inputs[1]->outputs.size() == 1)
      cur_node->inputs[1]->outputs.clear();
    else {
      NodeVector::iterator it = cur_node->inputs[1]->outputs.begin();
      while (it != cur_node->inputs[1]->outputs.end() && *it != cur_node) {
        ++it;
      }
      if (it != cur_node->inputs[1]->outputs.end())
        cur_node->inputs[1]->outputs.erase(it);
    }
    cur_node->inputs.clear();
    impl_cirs[cur_change_num]->Simplify();
    impl_cirs[cur_change_num]->SetIndexes();
    cout << "Creating thread " << cur_change_num << " for node type: " << cur_node->type << endl;
    sim_threads[cur_change_num] = std::thread(DoSimEqTh4, spec_cirs[cur_change_num], impl_cirs[cur_change_num], &po_diff_count[cur_change_num], &total_diff_count[cur_change_num]);

    cur_change_num++;

    spec_cirs[cur_change_num] = new_spec_cir.GetDuplicate("", "", "");
    impl_cirs[cur_change_num] = new_spec_cir.GetDuplicate("", "", "");
    cur_node = impl_cirs[cur_change_num]->all_nodes_map[target_node->name];
    cur_node->type = NODE_BUF;
    if (cur_node->inputs[1]->outputs.size() == 1)
      cur_node->inputs[1]->outputs.clear();
    else {
      NodeVector::iterator it = cur_node->inputs[1]->outputs.begin();
      while (it != cur_node->inputs[1]->outputs.end() && *it != cur_node) {
        ++it;
      }
      if (it != cur_node->inputs[1]->outputs.end())
        cur_node->inputs[1]->outputs.erase(it);
    }
    cur_node->inputs.resize(1);
    impl_cirs[cur_change_num]->Simplify();
    impl_cirs[cur_change_num]->SetIndexes();
    cout << "Creating thread " << cur_change_num << " for node type: " << cur_node->type << endl;
    sim_threads[cur_change_num] = std::thread(DoSimEqTh4, spec_cirs[cur_change_num], impl_cirs[cur_change_num], &po_diff_count[cur_change_num], &total_diff_count[cur_change_num]);

    cur_change_num++;

    spec_cirs[cur_change_num] = new_spec_cir.GetDuplicate("", "", "");
    impl_cirs[cur_change_num] = new_spec_cir.GetDuplicate("", "", "");
    cur_node = impl_cirs[cur_change_num]->all_nodes_map[target_node->name];
    cur_node->type = NODE_BUF;
    if (cur_node->inputs[0]->outputs.size() == 1)
      cur_node->inputs[0]->outputs.clear();
    else {
      NodeVector::iterator it = cur_node->inputs[0]->outputs.begin();
      while (it != cur_node->inputs[0]->outputs.end() && *it != cur_node) {
        ++it;
      }
      if (it != cur_node->inputs[0]->outputs.end())
        cur_node->inputs[0]->outputs.erase(it);
    }
    cur_node->inputs[0] = cur_node->inputs[1];
    cur_node->inputs.resize(1);
    impl_cirs[cur_change_num]->Simplify();
    impl_cirs[cur_change_num]->SetIndexes();
    cout << "Creating thread " << cur_change_num << " for node type: " << cur_node->type << endl;
    sim_threads[cur_change_num] = std::thread(DoSimEqTh4, spec_cirs[cur_change_num], impl_cirs[cur_change_num], &po_diff_count[cur_change_num], &total_diff_count[cur_change_num]);

    cur_change_num++;

    spec_cirs[cur_change_num] = new_spec_cir.GetDuplicate("", "", "");
    impl_cirs[cur_change_num] = new_spec_cir.GetDuplicate("", "", "");
    cur_node = impl_cirs[cur_change_num]->all_nodes_map[target_node->name];
    cur_node->type = NODE_NOT;
    if (cur_node->inputs[1]->outputs.size() == 1)
      cur_node->inputs[1]->outputs.clear();
    else {
      NodeVector::iterator it = cur_node->inputs[1]->outputs.begin();
      while (it != cur_node->inputs[1]->outputs.end() && *it != cur_node) {
        ++it;
      }
      if (it != cur_node->inputs[1]->outputs.end())
        cur_node->inputs[1]->outputs.erase(it);
    }
    cur_node->inputs.resize(1);
    impl_cirs[cur_change_num]->Simplify();
    impl_cirs[cur_change_num]->SetIndexes();
    cout << "Creating thread " << cur_change_num << " for node type: " << cur_node->type << endl;
    sim_threads[cur_change_num] = std::thread(DoSimEqTh4, spec_cirs[cur_change_num], impl_cirs[cur_change_num], &po_diff_count[cur_change_num], &total_diff_count[cur_change_num]);

    cur_change_num++;

    spec_cirs[cur_change_num] = new_spec_cir.GetDuplicate("", "", "");
    impl_cirs[cur_change_num] = new_spec_cir.GetDuplicate("", "", "");
    cur_node = impl_cirs[cur_change_num]->all_nodes_map[target_node->name];
    cur_node->type = NODE_NOT;
    if (cur_node->inputs[0]->outputs.size() == 1)
      cur_node->inputs[0]->outputs.clear();
    else {
      NodeVector::iterator it = cur_node->inputs[0]->outputs.begin();
      while (it != cur_node->inputs[0]->outputs.end() && *it != cur_node) {
        ++it;
      }
      if (it != cur_node->inputs[0]->outputs.end())
        cur_node->inputs[0]->outputs.erase(it);
    }
    cur_node->inputs[0] = cur_node->inputs[1];
    cur_node->inputs.resize(1);
    impl_cirs[cur_change_num]->Simplify();
    impl_cirs[cur_change_num]->SetIndexes();
    cout << "Creating thread " << cur_change_num << " for node type: " << cur_node->type << endl;
    sim_threads[cur_change_num] = std::thread(DoSimEqTh4, spec_cirs[cur_change_num], impl_cirs[cur_change_num], &po_diff_count[cur_change_num], &total_diff_count[cur_change_num]);

    cur_change_num++;

    for (int nn = NODE_AND; nn < NODE_BLIF; nn++) {
      if (nn == target_node->type)
        continue;
      spec_cirs[cur_change_num] = new_spec_cir.GetDuplicate("", "", "");
      impl_cirs[cur_change_num] = new_spec_cir.GetDuplicate("", "", "");
      cur_node = impl_cirs[cur_change_num]->all_nodes_map[target_node->name];
      cur_node->type = (NodeType)nn;
      cout << "Creating thread " << cur_change_num << " for node type: " << cur_node->type << endl;
      sim_threads[cur_change_num] = std::thread(DoSimEqTh4, spec_cirs[cur_change_num], impl_cirs[cur_change_num], &po_diff_count[cur_change_num], &total_diff_count[cur_change_num]);

      cur_change_num++;
    }

    for (int t = 0; t < num_changes; t++)
      sim_threads[t].join();

    for (int t = 0; t < num_changes; t++) {
      cout << endl << "Summary of final result:" << endl;
      cout << "---------- change " << t << " ----------" << endl;
      for (int i = 0; i < 8; i++)
        cout << total_diff_count[t][i] << "\t";
      cout << endl;
      for (int j = 0; j < new_spec_cir.outputs.size(); j++)
        cout << new_spec_cir.outputs[j]->name << " -> " << po_diff_count[t][j] << endl;
    }
  }
  else
    MSG("currently only working on 2 input gates!");

  return 0;

  //Circuit* new_impl_cir = new_spec_cir.GetDuplicate("", "", "");
  //target_node->type = NODE_XOR;
  // just for debug!
  //new_spec_cir.WriteBlif("new_spec.blif");
  //new_impl_cir->WriteBlif("new_impl.blif");

  //return DoSimEq(new_spec_cir, *new_impl_cir);
}


#include "cnf.h"
int MyTest(string filename) {
  Circuit test_cir1;
  Circuit test_cir2;
  Circuit test_cir3;

  //test_cir1.ReadBlif(filename, false, false, false);
  //test_cir1.Clear();

  auto start_time = chrono::system_clock::now();

  test_cir1.ReadBlif(filename, false, false, false);

  auto end_time = chrono::system_clock::now();
  chrono::duration<double> diff_time = end_time - start_time;
  cout << "read time1 is: " << setprecision(3) << diff_time.count() << " s" << endl;

  start_time = chrono::system_clock::now();

  test_cir2.ReadBlif(filename, false, true, false);

  end_time = chrono::system_clock::now();
  diff_time = end_time - start_time;
  cout << "read time2 is: " << setprecision(3) << diff_time.count() << " s" << endl;

  test_cir2.WriteBlif("sort-new.blif");

  start_time = chrono::system_clock::now();

  test_cir3.ReadBlif(filename, true, true, false);

  end_time = chrono::system_clock::now();
  diff_time = end_time - start_time;
  cout << "read time3 is: " << setprecision(3) << diff_time.count() << " s" << endl;

  test_cir3.WriteBlif("sort-old.blif");

  return 0;

  Circuit test_cir;

  test_cir.ReadVerilogGL(filename, false, false, false);

  test_cir.WriteVerilog("testout.v");
  test_cir.WriteBlif("testout.blif");

  test_cir.LevelizeSortTopological(false);
  test_cir.RemoveBufNot();
  test_cir.Simplify();
  test_cir.LevelizeSortTopological(false);

  test_cir.SetIndexes();

  test_cir.WriteVerilog("testout.simple.v");
  test_cir.WriteBlif("testout.simple.blif");

  //satcnf::Cnf test_cnf;
  //test_cnf.Convert2Cnf(&test_cir, NULL, NULL, false); // do not optimize as it is for test!

  //test_cnf.WriteDimacs("test.cnf");

  return 0;
}



// For big data analysis
/*
  notfouts(NodeVector v1, Node* node): whether node is fanout
  debug_msg; cout message

*/

/*
  *c1 = c2.GetDuplicate, when change c1, c2 not change
  c1 = c2, when change c1, c2 will change, same to c2
  *c1 = c2->GetDuplicate, when change c1, c2 not change

*/


//****************************//
int DoSimEq(Circuit &new_spec_cir, Circuit &new_impl_cir, double &WCAE,  double &MRE, bit64 &total_diff_count) {

  WCAE = 0;
  MRE = 0;
  double WRE = 0;
  vector<double> relative_err;
  Val64Vector pi_vect;
  Val64Vector po_vect_spec;
  Val64Vector po_vect_impl;
  pi_vect.resize(new_spec_cir.inputs.size());
  po_vect_spec.resize(new_spec_cir.outputs.size());
  po_vect_impl.resize(new_impl_cir.outputs.size());
  total_diff_count = 0;
  int pi_index = 0;
  while (pi_index < new_spec_cir.inputs.size() && pi_index < NUM_SIG_PAR_CHANGE) {
    pi_vect[pi_index] = PAT64[pi_index];
    pi_index++;
  }
  while (pi_index < new_spec_cir.inputs.size()) {
    pi_vect[pi_index] = ZERO64;
    pi_index++;
  }


  if (pi_vect.size() <= NUM_SIG_PAR_CHANGE) {
    new_spec_cir.Simulate(pi_vect, po_vect_spec);
    new_impl_cir.Simulate(pi_vect, po_vect_impl);
    int num_patterns = pow(2, pi_vect.size());
    double worst_case_ae = 0;
    bit64 diff_data = 0;

    vector<double> val_spec;
    vector<double> val_impl;

    val_spec.resize(num_patterns);
    val_impl.resize(num_patterns);


    for (int i = 0; i < po_vect_spec.size(); i++) {
      diff_data |= (po_vect_spec[i] ^ po_vect_impl[i]);
      bit64 cur_res_spec = po_vect_spec[i];
      bit64 cur_res_impl = po_vect_impl[i];
      for (int j = 0; j < num_patterns; j++) {
        if (cur_res_spec%2)
          val_spec[j] += pow(2, i);
        if (cur_res_impl%2)
          val_impl[j] += pow(2, i);
        cur_res_spec >>= 1;
        cur_res_impl >>= 1;
      }
    }

    int num_diff = 0;
    for (int i = 0; diff_data != 0 && i < num_patterns; i++) {
      if (diff_data%2)
        num_diff++;
      diff_data >>= 1;
    }

    total_diff_count = num_diff;


    for (int i = 0; i < val_spec.size(); i++) {
      double cur_res_spec = val_spec[i];
      double cur_res_impl = val_impl[i];
      worst_case_ae = abs(cur_res_spec - cur_res_impl) / pow(2, new_spec_cir.outputs.size());
      WCAE = max(worst_case_ae, WCAE);

      if (cur_res_spec == 0) {
        // relative_err.push_back(cur_res_impl);
        MRE += cur_res_impl;
        // WRE = max(WRE,cur_res_impl);
      } else {
        double this_err = abs(cur_res_spec - cur_res_impl) / cur_res_spec;
        MRE += this_err;
        // relative_err.push_back(this_err);
      }

    }
    // MAE = MAE / num_patterns;
    // std::sort(relative_err.begin(), relative_err.end());
    MRE = MRE / num_patterns;

  }
  else {
    bit64 num_sim_run = pow(2, new_spec_cir.inputs.size()-NUM_SIG_PAR_CHANGE);
    bit64 prev_sim_run = 0;
    bit64 num_sim_run_scale = num_sim_run >> 3; // 8 scales
    char progress_str[9] = {'.','.','.','.','.','.','.','.',0};
    double worst_case_ae = 0;
    bit64 num_diff = 0;

    for (bit64 cur_sim_run = 0; cur_sim_run < num_sim_run; ) {

      // if (num_sim_run_scale > 1 && cur_sim_run%num_sim_run_scale == 0) {
      //   cout << "up to now, the worst case absolute error rate is " << WCAE << endl;
      // }


      new_spec_cir.Simulate(pi_vect, po_vect_spec);
      new_impl_cir.Simulate(pi_vect, po_vect_impl);
      bit64 diff_data = 0;

      vector<double> val_spec;
      vector<double> val_impl;

      val_spec.resize(64);
      val_impl.resize(64);



      for (int i=0; i<po_vect_impl.size(); i++) {
        bit64 cur_diff_data = po_vect_spec[i] ^ po_vect_impl[i];
        diff_data |= cur_diff_data;
        bit64 cur_res_spec = po_vect_spec[i];
        bit64 cur_res_impl = po_vect_impl[i];
        for (int j = 0; j < 64; j++) {
          if (cur_res_spec%2)
            val_spec[j] += pow(2, i);
          if (cur_res_impl%2)
            val_impl[j] += pow(2, i);
          cur_res_spec >>= 1;
          cur_res_impl >>= 1;
        }
      }

      num_diff += Count64(diff_data);

      for (int i = 0; i < val_spec.size(); i++) {
        double cur_res_spec = val_spec[i];
        double cur_res_impl = val_impl[i];
        worst_case_ae = abs(cur_res_spec - cur_res_impl) / pow(2, new_spec_cir.outputs.size());
        WCAE = max(worst_case_ae, WCAE);
        if (cur_res_spec == 0) {
          // relative_err.push_back(cur_res_impl);
          MRE += cur_res_impl;
          // WRE = max(WRE,cur_res_impl);
        } else {
          double this_err = abs(cur_res_spec - cur_res_impl) / cur_res_spec;
          MRE += this_err;
          // relative_err.push_back(this_err);
        }
      }


      prev_sim_run = cur_sim_run;
      cur_sim_run++;
      if (cur_sim_run < num_sim_run) {
        bit64 diff_sim_run = prev_sim_run ^cur_sim_run;
        bit64 cur_sim_run_dup = cur_sim_run;
        bit64 pi_vect_index = NUM_SIG_PAR_CHANGE;
        while (diff_sim_run != 0) {
          if (diff_sim_run%2) {
            if (cur_sim_run_dup%2) {
              pi_vect[pi_vect_index] = ONE64;
            }
            else {
              pi_vect[pi_vect_index] = ZERO64;
            }
          }
          diff_sim_run >>= 1;
          cur_sim_run_dup >>= 1;
          pi_vect_index++;
        }
      }

    }

    double num_patterns = pow(2, new_spec_cir.inputs.size());
    std::sort(relative_err.begin(), relative_err.end());
    // MRE = (relative_err[relative_err.size()/2 - 1] + relative_err[relative_err.size()/2])/2;
    MRE /= num_patterns;
    total_diff_count = num_diff;
  }


  // cout << "Worst Relative Distance is: " << WRE << endl;




  return 0;

}

int WCAEcheck(Circuit &new_spec_cir, Circuit &new_impl_cir) {
  double WCAE,MRE;
  bit64 total_diff_count;
  DoSimEq(new_spec_cir, new_impl_cir, WCAE, MRE, total_diff_count);
  cout << "Worst Case Absolute Error Rate: " << WCAE << endl;
  cout << "Mean Relative Relative Error Rate: " << MRE << endl;

  long int size1 = 0;
  long int size2 = 0;
  for (int i=0; i<new_spec_cir.all_nodes.size(); i++) {
    if (new_spec_cir.all_nodes[i]->inputs.size()==2) size1++;
  }
  new_impl_cir.Simplify();
  // new_impl_cir.Simplify2();
  new_impl_cir.SetIndexes();
  for (int i=0; i<new_impl_cir.all_nodes.size(); i++) {
    if (new_impl_cir.all_nodes[i]->inputs.size()==2) size2++;
  }
  double red = (double)(1.0*size2)/ (double) size1;
  cout << "original size: " << size1 << " current size: " << size2 <<" reduced to be: " << red << endl;

  return 0;
}

void DoSimEqTh2(Circuit *new_spec_cir, Circuit *new_impl_cir, double *WCAE, double *MAE, bit64 *total_diff) {
  DoSimEq(*new_spec_cir, *new_impl_cir, *WCAE, *MAE, *total_diff);
}

int RndSimEq(Circuit &new_spec_cir, Circuit &new_impl_cir, int num_patterns) {
  if (new_spec_cir.inputs.size() == 0) {
    MSG("Cannot simulate circuit with no primary input!");
    return 0;
  }

  if (new_spec_cir.inputs.size() != new_impl_cir.inputs.size() || new_spec_cir.outputs.size() != new_impl_cir.outputs.size()) {
    MSG("Number of Pi/Po of two circuits are different!");
    return 0;
  }

  auto start_time = chrono::system_clock::now();

  std::srand(std::time(0));

  Val64Vector pi_vect;
  Val64Vector po_vect_spec;
  Val64Vector po_vect_impl;
  pi_vect.resize(new_spec_cir.inputs.size());
  po_vect_spec.resize(new_spec_cir.outputs.size());
  po_vect_impl.resize(new_impl_cir.outputs.size());
  double WCAE = 0;
  double MRE = 0;
  double WRE = 0;
  double ED = 0;
  vector<double> relative_err;

  bit64 num_sim_run = num_patterns/64;
  if (num_patterns%64)
    num_sim_run++;
  bit64 num_diff = 0;
  for (bit64 cur_sim_run = 0; cur_sim_run < num_sim_run; cur_sim_run++) {
    int pi_index = 0;
    while (pi_index < new_spec_cir.inputs.size()) {
      bit64 new_rnd = MyRand64();
      //cout << std::bitset<64>(new_rnd) << endl;
      // TODO: check it is different from previous simulation pattern(s)
      pi_vect[pi_index] = new_rnd;
      pi_index++;
    }
    new_spec_cir.Simulate(pi_vect, po_vect_spec);
    new_impl_cir.Simulate(pi_vect, po_vect_impl);
    // bit64 diff_data = 0;

    double worst_case_ae = 0;

    vector<double> val_spec;
    vector<double> val_impl;

    val_spec.resize(num_patterns);
    val_impl.resize(num_patterns);


    for (int i = 0; i < po_vect_spec.size(); i++) {
      bit64 cur_res_spec = po_vect_spec[i];
      bit64 cur_res_impl = po_vect_impl[i];
      for (int j = 0; j < num_patterns; j++) {
        if (cur_res_spec%2)
          val_spec[j] += pow(2, i);
        if (cur_res_impl%2)
          val_impl[j] += pow(2, i);
        cur_res_spec >>= 1;
        cur_res_impl >>= 1;
      }
    }

    for (int i = 0; i < val_spec.size(); i++) {
      double cur_res_spec = val_spec[i];
      double cur_res_impl = val_impl[i];
      worst_case_ae = abs(cur_res_spec - cur_res_impl) / pow(2, new_spec_cir.outputs.size());
      WCAE = max(worst_case_ae, WCAE);
      ED += abs(cur_res_spec - cur_res_impl);
      if (cur_res_spec == 0) {
        MRE += cur_res_impl;
        // relative_err.push_back(cur_res_impl);
        // WRE = max(WRE,cur_res_impl);
      } else {
        double this_err = abs(cur_res_spec - cur_res_impl) / cur_res_spec;
        MRE += this_err;
        // WRE = max(WRE,this_err);
        // relative_err.push_back(this_err);
      }
    }

  }
  MRE = MRE / num_patterns;
  ED = ED / num_patterns;
  cout << "The worst case absolute error rate is: " << WCAE << " (out of " << num_sim_run*64 << " runs)" << endl;
  cout << "Mean Relative error rate is: " << MRE << " (out of " << num_sim_run*64 << " runs)" << endl;
  cout << "Average Error Distance is: " << ED << " (out of " << num_sim_run*64 << " runs)" << endl;
  // cout << "Worst Relative distance is: " << WRE << " (out of " << num_sim_run*64 << " runs)" << endl;
  auto end_time = chrono::system_clock::now();
  chrono::duration<double> diff_time = end_time - start_time;
  //

  long int size1 = 0;
  long int size2 = 0;
  for (int i=0; i<new_spec_cir.all_nodes.size(); i++) {
    if (new_spec_cir.all_nodes[i]->inputs.size()==2) size1++;
  }
  new_impl_cir.Simplify();
  // new_impl_cir.Simplify2();
  new_impl_cir.SetIndexes();
  for (int i=0; i<new_impl_cir.all_nodes.size(); i++) {
    if (new_impl_cir.all_nodes[i]->inputs.size()==2) size2++;
  }
  double red = (double)(1.0*size2)/ (double) size1;
  cout << "original size: " << size1 << " current size: " << size2 <<" reduced to be: " << red << endl;
  cout << endl << "Simulation time is: " << setprecision(3) << diff_time.count() << " s" << endl;
  return 0;
}

bool isInVector(vector<string> v1, string s1) {
  vector<string>::iterator it;
  for(it=v1.begin(); it!=v1.end(); it++) {
    if (*it == s1) return true;
  }
  return false;
}

void fanouts(Node* n1, set<string> &effect_nodes){
  long int n1_size = n1->outputs.size();
  for (int i = 0; i < n1_size; i++) {
    effect_nodes.insert(n1->outputs[i]->name);
    fanouts(n1->outputs[i], effect_nodes);
  }
}

bool notfouts(set<string> effect_nodes, string node_name) {
  std::set<string>::iterator it;
  it = effect_nodes.find(node_name);
  if (it == effect_nodes.end()) return true;
  return false;

}

void debug_msg(string kk){
  cout << kk << endl;
}

int whichOutput(Node* n1, vector<string> &outputName, vector<string> outputs){
  for (int i=0; i<n1->outputs.size(); i++) {
    if (isInVector(outputs, n1->outputs[i]->name)) {
      outputName.push_back(n1->outputs[i]->name);
    } else {
      whichOutput(n1->outputs[i], outputName, outputs);
    }
  }
  return 0;
}

void findEffectNodes(Circuit *spec_cir, string outName, set<string> &effect_nodes) {
  Node *outNode = spec_cir->all_nodes_map[outName];
  Node *cur_node = new Node;
  if (outNode->is_output) {
    effect_nodes.insert(cur_node->name);
  }
  else {
    for(int i=0; i<outNode->outputs.size(); i++) {
      cur_node = outNode->outputs[i];
      if (cur_node->is_output) effect_nodes.insert(cur_node->name);
      else findEffectNodes(spec_cir, cur_node->name, effect_nodes);

    }
  }
}

void findEffectNodes(Circuit *spec_cir, string outName, vector<string> &effects) {
  effects.clear();
  set<string> tmp;
  findEffectNodes(spec_cir, outName, tmp);
  set<string>::iterator it;

  for (it=tmp.begin(); it!=tmp.end(); it++) {
    effects.push_back(*it);
  }
}

void GetPossibleNodes(Circuit* new_spec_cir, string first_node_name, NodeVector &second_nodes, int mode){
  Node* cur_node = NULL;
  Node* c1 = NULL;
  Circuit* spec_cir = new_spec_cir->GetDuplicate("", "", "");
  second_nodes.clear();
  if (first_node_name !=""){
    NodeVector tar_node;
    cur_node = spec_cir->all_nodes_map[first_node_name];
    tar_node.push_back(cur_node);
    spec_cir->Simplify(tar_node);
    spec_cir->Simplify2();
    spec_cir->ResetLevels();
    spec_cir->LevelizeSortTopological(false);
    spec_cir->SetIndexes();
  } else {
    spec_cir->Simplify();
    spec_cir->Simplify2();
    spec_cir->SetIndexes();
  }
  switch (mode) {
    case 0:
      for (long i=0; i<spec_cir->all_nodes.size(); i++){
        c1 = spec_cir->all_nodes[i];
        if(c1->inputs.size()==2 && !c1->is_input && c1->name!=first_node_name){
          second_nodes.push_back(c1);
        }
      }
      break;
    case 1:
      for (long i=0; i<spec_cir->all_nodes.size(); i++){
        c1 = spec_cir->all_nodes[i];
        if(first_node_name==""){
          if(c1->inputs.size()==2 && !c1->is_input){
            second_nodes.push_back(c1);
          }
        }
        if(first_node_name!="" && c1->name!=first_node_name) continue;
        if(c1->name==first_node_name) first_node_name="";
      }
      break;
  }
  // delete spec_cir;
}

void RanSelectNodes(NodeVector source_nodes, int amount, NodeVector &select_nodes){
  int total_nodes_amount = source_nodes.size();
  select_nodes.clear();
  srand((unsigned)time( NULL));
  int seed = 0;
  int seedup = 0;
  int seeddown = 0;
  Node* cur_node;
  int level_size = 3;
  if (total_nodes_amount<amount) amount = total_nodes_amount;
  int level_amount = ceil(total_nodes_amount/3);
  for (int i=0; i<level_size; i++){
    level_amount = ceil(total_nodes_amount/3);
    seedup = level_amount * (i+1) > total_nodes_amount ? total_nodes_amount-1 : level_amount * (i+1);
    seeddown = seedup - level_amount;
    while (amount>0 and level_amount>0){
      seed = rand()%(seedup - seeddown) + seeddown;
      cur_node = source_nodes[seed];
      select_nodes.push_back(cur_node);
      amount --;
      level_amount--;
    }
  }
}

void GetPossibleConnection(Circuit *spec_cir, string current_node_name, NodeVector &connections){
  NodeVector possible_connections;
  connections.clear();
  Node* current_node = spec_cir->all_nodes_map[current_node_name];
  Node* n1 = NULL;
  NodeVector tar_node;
  tar_node.push_back(current_node);
  spec_cir->Simplify(tar_node);
  spec_cir->ResetLevels();
  spec_cir->LevelizeSortTopological(false);
  spec_cir->SetIndexes();

  std::set<string> effect_nodes;
  fanouts(current_node, effect_nodes);

  for (long i=0; i<spec_cir->all_nodes.size(); i++){
    n1 = spec_cir->all_nodes[i];
    if(notfouts(effect_nodes, n1->name) and current_node->name != n1->name and current_node->inputs[0]->name != n1->name and n1->name != current_node->inputs[1]->name and !n1->is_input and !n1->is_output){
      possible_connections.push_back(n1);
    }
  }
  connections = possible_connections;
}

void GetWrongCircuit(Circuit &new_spec_cir, string target_node_name, int type){
  Node* target_node = new_spec_cir.all_nodes_map[target_node_name];
  if (target_node->inputs.size() == 2) {
    Node* cur_node = new Node;
    switch(type){
      case 0:
        cur_node = new_spec_cir.all_nodes_map[target_node->name];
        cur_node->type = NODE_ZERO;
        if (cur_node->inputs[0]->outputs.size() == 1)
          cur_node->inputs[0]->outputs.clear();
        else {
          NodeVector::iterator it = cur_node->inputs[0]->outputs.begin();
          while (it != cur_node->inputs[0]->outputs.end() && *it != cur_node) {
            ++it;
          }
          if (it != cur_node->inputs[0]->outputs.end())
            cur_node->inputs[0]->outputs.erase(it);
        }
        if (cur_node->inputs[1]->outputs.size() == 1)
          cur_node->inputs[1]->outputs.clear();
        else {
          NodeVector::iterator it = cur_node->inputs[1]->outputs.begin();
          while (it != cur_node->inputs[1]->outputs.end() && *it != cur_node) {
            ++it;
          }
          if (it != cur_node->inputs[1]->outputs.end())
            cur_node->inputs[1]->outputs.erase(it);
        }
        cur_node->inputs.clear();
        break;
      case 1:
        cur_node = new_spec_cir.all_nodes_map[target_node->name];
        cur_node->type = NODE_ONE;
        if (cur_node->inputs[0]->outputs.size() == 1)
          cur_node->inputs[0]->outputs.clear();
        else {
          NodeVector::iterator it = cur_node->inputs[0]->outputs.begin();
          while (it != cur_node->inputs[0]->outputs.end() && *it != cur_node) {
            ++it;
          }
          if (it != cur_node->inputs[0]->outputs.end())
            cur_node->inputs[0]->outputs.erase(it);
        }
        if (cur_node->inputs[1]->outputs.size() == 1)
          cur_node->inputs[1]->outputs.clear();
        else {
          NodeVector::iterator it = cur_node->inputs[1]->outputs.begin();
          while (it != cur_node->inputs[1]->outputs.end() && *it != cur_node) {
            ++it;
          }
          if (it != cur_node->inputs[1]->outputs.end())
            cur_node->inputs[1]->outputs.erase(it);
        }
        cur_node->inputs.clear();
        break;
      case 2:
        cur_node = new_spec_cir.all_nodes_map[target_node->name];
        cur_node->type = NODE_BUF;
        if (cur_node->inputs[1]->outputs.size() == 1)
          cur_node->inputs[1]->outputs.clear();
        else {
          NodeVector::iterator it = cur_node->inputs[1]->outputs.begin();
          while (it != cur_node->inputs[1]->outputs.end() && *it != cur_node) {
            ++it;
          }
          if (it != cur_node->inputs[1]->outputs.end())
            cur_node->inputs[1]->outputs.erase(it);
        }
        cur_node->inputs.resize(1);
        break;
      case 3:
        cur_node = new_spec_cir.all_nodes_map[target_node->name];
        cur_node->type = NODE_BUF;
        if (cur_node->inputs[0]->outputs.size() == 1)
          cur_node->inputs[0]->outputs.clear();
        else {
          NodeVector::iterator it = cur_node->inputs[0]->outputs.begin();
          while (it != cur_node->inputs[0]->outputs.end() && *it != cur_node) {
            ++it;
          }
          if (it != cur_node->inputs[0]->outputs.end())
            cur_node->inputs[0]->outputs.erase(it);
        }
        cur_node->inputs[0] = cur_node->inputs[1];
        cur_node->inputs.resize(1);
        break;
      case 4:
        cur_node = new_spec_cir.all_nodes_map[target_node->name];
        cur_node->type = NODE_NOT;
        if (cur_node->inputs[1]->outputs.size() == 1)
          cur_node->inputs[1]->outputs.clear();
        else {
          NodeVector::iterator it = cur_node->inputs[1]->outputs.begin();
          while (it != cur_node->inputs[1]->outputs.end() && *it != cur_node) {
            ++it;
          }
          if (it != cur_node->inputs[1]->outputs.end())
            cur_node->inputs[1]->outputs.erase(it);
        }
        cur_node->inputs.resize(1);
        break;
      case 5:
        cur_node = new_spec_cir.all_nodes_map[target_node->name];
        cur_node->type = NODE_NOT;
        if (cur_node->inputs[0]->outputs.size() == 1)
          cur_node->inputs[0]->outputs.clear();
        else {
          NodeVector::iterator it = cur_node->inputs[0]->outputs.begin();
          while (it != cur_node->inputs[0]->outputs.end() && *it != cur_node) {
            ++it;
          }
          if (it != cur_node->inputs[0]->outputs.end())
            cur_node->inputs[0]->outputs.erase(it);
        }
        cur_node->inputs[0] = cur_node->inputs[1];
        cur_node->inputs.resize(1);
        break;
      default:
        type--;
        for (int nn = NODE_AND; nn < NODE_BLIF; nn++) {
          if(nn == type){
            cur_node = new_spec_cir.all_nodes_map[target_node->name];
            cur_node->type = (NodeType)nn;
            break;
          }
        }
        break;
    }
  }else MSG("currently only working on 2 input gates!");
}

void GetWrongCircuit(Circuit *spec_cir, vector<Circuit*> &impl_cirs, string target_node_name){
  Node *target_node = spec_cir->all_nodes_map[target_node_name];
  int cur_num = 0;
  impl_cirs.clear();
  impl_cirs.resize(15);

  for (int jj = 0; jj < 16; jj++){
    if (jj-1 == target_node->type) continue;
    impl_cirs[cur_num] = spec_cir->GetDuplicate("","","");
    GetWrongCircuit(*impl_cirs[cur_num], target_node_name, jj);
    cur_num++;
  }
}

void GetHeuristic(Circuit *spec_cir, vector<Circuit*> &impl_cirs, string target_node_name){
  Node *target_node = spec_cir->all_nodes_map[target_node_name];
  int cur_num = 0;
  impl_cirs.clear();
  impl_cirs.resize(6);

  for (int jj = 0; jj < 6; jj++){
    if (jj-1 == target_node->type) continue;
    impl_cirs[cur_num] = spec_cir->GetDuplicate("","","");
    GetWrongCircuit(*impl_cirs[cur_num], target_node_name, jj);
    cur_num++;
  }
}

int RndSimEq(Circuit *new_spec_cir, Circuit *new_impl_cir, int num_patterns, double &total_diff) {

  total_diff = 0;

  if (new_spec_cir->inputs.size() == 0) {
    MSG("Cannot simulate circuit with no primary input!");
    return 0;
  }

  if (new_spec_cir->inputs.size() != new_impl_cir->inputs.size() || new_spec_cir->outputs.size() != new_impl_cir->outputs.size()) {
    MSG("Number of Pi/Po of two circuits are different!");
    return 0;
  }

  auto start_time = chrono::system_clock::now();

  std::srand(std::time(0));

  Val64Vector pi_vect;
  Val64Vector po_vect_spec;
  Val64Vector po_vect_impl;
  pi_vect.resize(new_spec_cir->inputs.size());
  po_vect_spec.resize(new_spec_cir->outputs.size());
  po_vect_impl.resize(new_impl_cir->outputs.size());

  bit64 num_sim_run = num_patterns/64;
  if (num_patterns%64)
    num_sim_run++;
  bit64 num_diff = 0;
  for (bit64 cur_sim_run = 0; cur_sim_run < num_sim_run; cur_sim_run++) {
    int pi_index = 0;
    while (pi_index < new_spec_cir->inputs.size()) {
      bit64 new_rnd = MyRand64();
      pi_vect[pi_index] = new_rnd;
      pi_index++;
    }
    new_spec_cir->Simulate(pi_vect, po_vect_spec);
    new_impl_cir->Simulate(pi_vect, po_vect_impl);
    bit64 diff_data = 0;
    for (int i = 0; i < po_vect_spec.size(); i++) {
      diff_data |= (po_vect_spec[i] ^ po_vect_impl[i]);
    }
    if (diff_data == ONE64)
      num_diff += 64;
    else
      for (int i = 0; diff_data != 0 && i < 64; i++) {
        if (diff_data % 2)
          num_diff++;
        diff_data >>= 1;
      }

  }
  cout << "Number of different minterms: " << num_diff << " (out of " << num_sim_run*64 << " runs)" << endl;
  auto end_time = chrono::system_clock::now();
  chrono::duration<double> diff_time = end_time - start_time;
  cout << endl << "Random Simulation time is: " << setprecision(3) << diff_time.count() << " s" << endl;
  total_diff = num_diff;

  return 0;
}

void doRndTh4(Circuit *spec_cir_1, Circuit *spec_cir_2, int num_patterns, double *num_diff) {
  RndSimEq(spec_cir_1, spec_cir_2, num_patterns, *num_diff);
}

int DoSimEqModRes(Circuit &new_spec_cir, NodeVector target_nodes, string spec_filename, int mode) {
  ofstream outFile;
  ofstream outText;

  string filename = spec_filename;
  size_t dotp = filename.find('.');
  filename = filename.substr(0,dotp);

  switch(mode) {
    case 0: {
      cout << "this mode is for testeq case 5" << endl;

    }break;
    case 1: {
      cout << "Now this mode is in maintence" << endl;
    }break;
    case 2: {
      outFile.open(filename+".csv",  ios::out);

      const int thread_num = 15;
      int numpatterns = 50000;
      double total_diff[thread_num];

      Circuit* spec_cir = new_spec_cir.GetDuplicate("","","");
      std::vector<Circuit*> impl_cirs;
      std::thread sim_threads[thread_num];

      string target_node_name = "";

      long int node_size = target_nodes.size();

      for (int i = 0; i < node_size; i++) {
        target_node_name = target_nodes[i]->name;
        spec_cir = new_spec_cir.GetDuplicate("","","");

        GetWrongCircuit(spec_cir, impl_cirs, target_node_name);

        for (int j = 0; j < thread_num; j++){
          sim_threads[j] = std::thread(doRndTh4, spec_cir, impl_cirs[j], numpatterns, &total_diff[j]);
        }

        for (int j = 0; j < thread_num; j++){
          sim_threads[j].join();
        }

        for (int j = 0; j < thread_num; j++){
          // cout << target_node_name << "," << impl_cirs[j]->all_nodes_map[target_node_name]->level << "," << spec_cir->all_nodes_map[target_node_name]->type << "," << impl_cirs[j]->all_nodes_map[target_node_name]->type << "," << total_diff[i] << "," << numpatterns << endl;
          outFile << target_node_name << "," << impl_cirs[j]->all_nodes_map[target_node_name]->level << "," << spec_cir->all_nodes_map[target_node_name]->type << "," << impl_cirs[j]->all_nodes_map[target_node_name]->type << "," << total_diff[j] << "," << numpatterns << endl;
        }
        delete spec_cir;
        for (int i = 0; i < thread_num; i++) {
          delete impl_cirs[i];
        }
      }

      outFile.close();
    }break;
    default: {
      cout << "this mode is in maintence" << endl;
    }break;

    return 0;
  }





  // if(mode==1 || mode==4){
  //   outFile.open(filename+"withRes.csv",  ios::out | ios::app);
  //   outText.open(filename+"withRes.txt", ios::out | ios::app);
  // }
  // if(mode==2){
  //   outFile.open(filename+".csv",  ios::out | ios::app);
  // }
  // if(mode==3){
  //   outFile.open(filename+"random_analysis.csv",  ios::out | ios::app);
  //   outText.open(filename+"random_analysis.csv",  ios::out | ios::app);
  // }
  //
  //
  // const int num_changes = 15;
  // Circuit* spec_cirs[num_changes];
  // Circuit* impl_cirs[num_changes];
  // std::thread sim_threads[num_changes];
  // vector<bit64> po_diff_count[num_changes];
  // vector<bit64> total_diff_count[num_changes];
  //
  // int numpatterns;
  // int pattern_nums[5] = {1000,5000,10000,20000,50000};
  // int numberofwrong[5] = {0, 0, 0, 0, 0};
  // int threadnum = 0;
  // std::thread rand2ex_thread[15];
  // string nodename[15]={"","","","","","","","","","","","","","","",};
  // long nodelevel[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  // int nodechange[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  // double inputsizechange[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  // vector<bit64> total_node_count[15];
  // Node* cur_node = NULL;
  // int cur_change_num;
  // Circuit* spec_cir;
  // vector<int> real_sim_num;
  // string target_node_name;
  // Node* target_node = NULL;
  // long int size1;



  // for(int m = 0; m < targets_nodes.size(); m++){
  //
  //   target_node_name = targets_nodes[m]->name;
  //   target_node = new_spec_cir.all_nodes_map[target_node_name];
  //   NodeVector target_nodes;
  //   target_nodes.push_back(target_node);
  //
  //   size1 = new_spec_cir.inputs.size();
  //
  //   if(mode==1 || mode==4) {
  //     new_spec_cir.Simplify(target_nodes);
  //     new_spec_cir.ResetLevels();
  //     new_spec_cir.LevelizeSortTopological(false);
  //     new_spec_cir.SetIndexes();
  //   }
  //
  //
  //   double needtoadd = 1;
  //
  //   if(mode==1){
  //     outText<<"change node is "<<target_node->name<<" node level is "<<target_node->level<<" PIs "<< new_spec_cir.inputs.size()<<endl;
  //     needtoadd = pow(2, size1 - new_spec_cir.inputs.size());
  //   }
  //
  //
  //   for (int i=0; i<num_changes; i++){
  //     spec_cirs[i] = new_spec_cir.GetDuplicate("","","");
  //     impl_cirs[i] = new_spec_cir.GetDuplicate("","","");
  //     GetWrongCircuit(*impl_cirs[i], target_node_name, i);
  //     // sim_threads[i] = std::thread(DoSimEqTh4, spec_cirs[i], impl_cirs[i], &po_diff_count[i], &total_diff_count[i]);
  //   }


  //
  //   if (mode==1 || mode==4){
  //     for(int t=0; t<num_changes; t++){
  //       numpatterns = 10000;
  //       int h = RndSimEq(spec_cirs[t], impl_cirs[t], numpatterns);
  //       outText<<" change "<< t << " wrong cases "<< h << "in "<<numpatterns<<" 's input pattern cases"<<endl;
  //       if(h<=100 && threadnum<15) {
  //         total_node_count[threadnum] = total_diff_count[t];
  //         rand2ex_thread[threadnum] = std::thread(DoSimEqTh4, spec_cirs[t], impl_cirs[t], &po_diff_count[t], &total_node_count[threadnum]);
  //         nodename[threadnum] = target_node->name;
  //         nodelevel[threadnum] = target_node->level;
  //         nodechange[threadnum] = t;
  //         inputsizechange[threadnum] = needtoadd;
  //         threadnum ++;
  //       }
  //       if(threadnum==15 || m==targets_nodes.size()-1){
  //         cout<<"start thread"<<endl;
  //         for (int t=0; t<threadnum; t++){
  //           cout<<"node "<<nodename[t]<<" level "<<nodelevel[t]<<" change "<<nodechange[t]<<" has be sent to accurate simulation"<<endl;
  //           outText<<"node "<<nodename[t]<<" level "<<nodelevel[t]<<" change "<<nodechange[t]<<" has be sent to accurate simulation"<<endl;
  //         }
  //         for (int t=0; t<threadnum; t++){
  //           rand2ex_thread[t].join();
  //         }
  //         if(mode==1){
  //           for (int t = 0; t < threadnum; t++) {
  //             cout << endl << "Summary of final result:" << endl;
  //             cout << "-------"<<nodename[t]<<"-----"<<nodelevel[t] << "--------"<<nodechange[t]<<"----------" << endl;
  //             for (int i = 0; i < 8; i++){
  //               cout << total_node_count[t][i] << "\t";
  //             }
  //             cout<<"total differnces is "<< total_node_count[t][7]<<endl;
  //             outFile<< nodename[t] <<","<<nodelevel[t]<<","<<nodechange[t]<<","<<total_node_count[t][7]*inputsizechange[t]<<endl;
  //           }
  //         }
  //         threadnum = 0;
  //       }
  //     }
  //   }
  //
  //   if(mode==2 ){
  //     for(int t=0; t<num_changes; t++){ //ALL RANDOM
  //       numpatterns = 50000;
  //       int h = RndSimEq(spec_cirs[t], impl_cirs[t], numpatterns);
  //       cout<<"node is: "<<target_node->name <<" node level is "<<target_node->level<<" change "<< t<< " wrong cases "<< h << "in "<<numpatterns<<" 's input pattern cases"<<endl;
  //       outText<<"node is: "<<target_node->name <<" node level is "<<target_node->level<<" change "<< t<< " wrong cases "<< h << "in "<<numpatterns<<" 's input pattern cases"<<endl;
  //       outFile<<target_node->name <<","<<target_node->level<<","<<t<<","<<h<<","<<numpatterns<<endl;
  //     }
  //   }
  //
  //   if(mode==3){ //ALL RANDOM UNDER DIFFERENT PATTERNS
  //     for(int t=0; t<num_changes; t++){
  //       int h;
  //       cout<<"node is: "<<target_node->name <<" node level is "<<target_node->level<<" change "<< t;
  //       outText<<"node is: "<<target_node->name <<" node level is "<<target_node->level<<" change "<< t;
  //       outFile<<target_node->name <<","<<target_node->level<<","<<t<<",";
  //       for(int j=0; j<5; j++){
  //         h = RndSimEq(spec_cirs[t], impl_cirs[t], pattern_nums[j]);
  //         if(h==0) numberofwrong[j] += 1;
  //         cout<< " wrong cases "<< h << "in "<<pattern_nums[j]<<" 's input pattern cases"<<"\t";
  //         outText<< " wrong cases "<< h << "in "<<pattern_nums[j]<<" 's input pattern cases"<<"\t";
  //         outFile<<h<<",";
  //       }
  //       cout<<endl;
  //       outText<<endl;
  //       outFile<<endl;
  //     }
  //   }
  // }
  // if(mode==3){
  //   outFile<<""<<","<<""<<","<<""<<",";
  //   for(int j=0; j<5; j++){
  //     outFile<<numberofwrong[j]<<",";
  //   }
  //   outFile<<endl;
  // }

  outText.close();
  outFile.close();
  return 0;

}

int checkPattern(vector<bit64> total_diff_count){
  long standard = total_diff_count[0];
  for(int i=2; i<total_diff_count.size(); i++){
    if (Count64(total_diff_count[i] - total_diff_count[i-1])!= Count64(standard) || Count64(total_diff_count[i] - total_diff_count[i-1])!= Count64(standard*2)){
      return 0;
    }
  }
  return 1;
}

void DoTestEq(Circuit &new_spec_cir, double &min_err) {

  min_err = 256;

  const int thread_num = 15;

  Circuit* spec_cir = new_spec_cir.GetDuplicate("","","");
  std::vector<Circuit*> impl_cirs;
  std::thread sim_threads[thread_num];
  vector<bit64> po_diff_count[thread_num];
  vector<bit64> total_diff_count[thread_num];

  NodeVector target_nodes;
  GetPossibleNodes(spec_cir, "", target_nodes, 0);

  long int original_size = new_spec_cir.inputs.size();

  cout << pow(2, original_size) << endl;



  cout << min_err << endl;

  double needtoadd = 1;
  long int node_length = target_nodes.size();


  string target_node_name;

  NodeVector change_nodes;

  Node* target_node = new Node;


  for(long int i = 0; i < node_length; i++){

    target_node_name = target_nodes[i]->name;
    spec_cir = new_spec_cir.GetDuplicate("","","");

    target_node = spec_cir->all_nodes_map[target_node_name];

    change_nodes.clear();
    change_nodes.push_back(target_node);

    spec_cir->Simplify(change_nodes);
    spec_cir->ResetLevels();
    spec_cir->LevelizeSortTopological(false);
    spec_cir->SetIndexes();

    needtoadd = pow(2, (original_size - spec_cir->inputs.size()));

    GetWrongCircuit(spec_cir, impl_cirs, target_node_name);

    for (int ii = 0; ii < thread_num; ii++){
      sim_threads[ii] = std::thread(DoSimEqTh4, spec_cir, impl_cirs[ii], &po_diff_count[ii], &total_diff_count[ii]);
    }


    for (int ii = 0; ii < thread_num; ii++) {
      sim_threads[ii].join();
    }

    for (int t = 0; t < thread_num; t++) {
      cout << target_node->name << "," << target_node->level << "," << target_node->type << "," << impl_cirs[t]->all_nodes_map[target_node_name]->type << "," << total_diff_count[t][7]*needtoadd << endl;
      if (total_diff_count[t][7]>0) min_err = min(min_err, total_diff_count[t][7]*needtoadd);
    }
  }

  delete spec_cir;
  for (int ii = 0; ii < thread_num; ii++) {
    delete impl_cirs[ii];
  }
}

void DoSimEqMod(Circuit &new_spec_cir, NodeVector target_nodes, string spec_filename) {

  ofstream outFile;
  string filename = spec_filename;
  size_t dotp = filename.find('.');
  filename = filename.substr(0,dotp);
  outFile.open(filename+".csv", ios::out);


  const int thread_num = 15;

  Circuit* spec_cir;
  std::vector<Circuit*> impl_cirs;
  std::thread sim_threads[thread_num];
  vector<bit64> po_diff_count[thread_num];
  vector<bit64> total_diff_count[thread_num];

  long int original_size = new_spec_cir.inputs.size();
  double needtoadd = 1;
  long int node_length = target_nodes.size();


  string target_node_name;

  NodeVector change_nodes;

  Node* target_node = new Node;


  for(long int i = 0; i < node_length; i++){

    target_node_name = target_nodes[i]->name;
    spec_cir = new_spec_cir.GetDuplicate("","","");

    target_node = spec_cir->all_nodes_map[target_node_name];

    change_nodes.clear();
    change_nodes.push_back(target_node);


    spec_cir->Simplify(change_nodes);
    spec_cir->ResetLevels();
    spec_cir->LevelizeSortTopological(false);
    spec_cir->SetIndexes();

    needtoadd = pow(2, (original_size - spec_cir->inputs.size()));

    // spec_cir->WriteBlif("spec.blif");

    GetWrongCircuit(spec_cir, impl_cirs, target_node_name);


    // for (int ii = 0; ii < thread_num; ii++) {
    //   cout << target_node->name << "," << target_node->level << "," << target_node->type << "," << impl_cirs[ii]->all_nodes_map[target_node_name]->type << endl;
    //   cout << impl_cirs[ii]->all_nodes.size() << endl;
    // }


    for (int ii = 0; ii < thread_num; ii++){
      // DoSimEq(*spec_cir, *impl_cirs[ii], po_diff_count[ii], total_diff_count[ii]);
      // cout << total_diff_count[ii][7] << endl;
      sim_threads[ii] = std::thread(DoSimEqTh4, spec_cir, impl_cirs[ii], &po_diff_count[ii], &total_diff_count[ii]);
    }




    for (int ii = 0; ii < thread_num; ii++) {
      sim_threads[ii].join();
    }

    for (int t = 0; t < thread_num; t++) {
      cout << target_node->name << "," << target_node->level << "," << target_node->type << "," << impl_cirs[t]->all_nodes_map[target_node_name]->type << "," << total_diff_count[t][7]*needtoadd << endl;
      outFile << target_node->name << "," << target_node->level << "," << target_node->type << "," << impl_cirs[t]->all_nodes_map[target_node_name]->type << "," << total_diff_count[t][7]*needtoadd << endl;
    }

    delete spec_cir;

    for (int ii = 0; ii < thread_num; ii++) {
      // impl_cirs[ii]->WriteBlif("spec"+to_string(ii)+".blif");
      delete impl_cirs[ii];
    }
  }
}

int DoSimEqMulGate(Circuit &new_spec_cir, string spec_filename, NodeVector first_nodes, int mode){
  ofstream outFile;
  ofstream outMean;
  string filename = spec_filename;
  size_t dotp = filename.find('.');
  filename = filename.substr(0,dotp);
  if(first_nodes.size()==1) {
    outFile.open("test2.csv",  ios::out | ios::app);
    outMean.open("testMulG.csv", ios::out | ios::app);
  }

  else {
    outFile.open(filename+"_2G.csv",  ios::out | ios::app);
    outMean.open(filename+"wished.csv",  ios::out | ios::app);
  }

  const int num_changes = 15;

  std::vector<Circuit*> impl_cirs;
  Circuit* origin_cir = new_spec_cir.GetDuplicate("", "", "");
  Circuit* impl_cir;
  Circuit* spec_cir;
  std::vector<Circuit*> spec_cirs;
  std::thread sim_threads[num_changes];
  vector<bit64> po_diff_count[num_changes];
  vector<bit64> total_diff_count[num_changes];

  Node* first_node = NULL;
  Node* second_node = NULL;
  string first_node_name;
  string second_node_name;

  double needtoadd = 1;
  double max_error = 0;
  NodeVector second_nodes;
  NodeVector target_nodes;


  //mode 0
  // simplifying previous code, use GetWrongCircuit function to get 15 circuits at the same time

  switch(mode) {
    case 0:{

      for (int i=0; i<first_nodes.size()-2; i++){

        first_node_name = first_nodes[i]->name;
        // GetPossibleNodes(origin_cir, first_node_name, second_nodes, 1);

        for(int j=i+1; j<first_nodes.size()-1; j++){

          second_node_name = first_nodes[j]->name;

          spec_cir = new_spec_cir.GetDuplicate("","","");

          first_node = spec_cir->all_nodes_map[first_node_name];
          second_node = spec_cir->all_nodes_map[second_node_name];

          target_nodes.clear();
          target_nodes.push_back(first_node);
          target_nodes.push_back(second_node);

          spec_cir->Simplify(target_nodes);
          spec_cir->ResetLevels();
          spec_cir->LevelizeSortTopological(false);
          spec_cir->SetIndexes();

          if (spec_cir->all_nodes_map.find(first_node_name) != spec_cir->all_nodes_map.end() && spec_cir->all_nodes_map.find(second_node_name) != spec_cir->all_nodes_map.end()){

            needtoadd = pow(2, new_spec_cir.inputs.size() - spec_cir->inputs.size());
            max_error = ceil(pow(2, spec_cir->inputs.size()) * 0.01);

            GetWrongCircuit(spec_cir, spec_cirs, first_node_name);

            for (int ii=0; ii < num_changes; ii++) {
              impl_cir = spec_cirs[ii]->GetDuplicate("","","");
              if (spec_cir->all_nodes_map[first_node_name]->type==5 && impl_cir->all_nodes_map[first_node_name]->type>5) {
                continue;
              }
              GetWrongCircuit(impl_cir, impl_cirs, second_node_name);


              for (int jj=0; jj < num_changes; jj++) {
                sim_threads[jj] = std::thread(DoSimEqTh4, spec_cir, impl_cirs[jj], &po_diff_count[jj], &total_diff_count[jj]);
              }
              for (int jj=0; jj < num_changes; jj++) {
                sim_threads[jj].join();
              }

              for (int t=0; t < num_changes; t++) {
                if (spec_cir->all_nodes_map[second_node_name]->type ==5 && impl_cirs[t]->all_nodes_map[second_node_name]->type>5) continue;
                cout << first_node_name << "," << first_node->type << "," << impl_cirs[t]->all_nodes_map[first_node_name]->type << "," << second_node_name << "," << second_node->type << "," << impl_cirs[t]->all_nodes_map[second_node_name]->type << "," << total_diff_count[t][7]*needtoadd << endl;
                outFile << first_node_name << "," << first_node->type << "," << impl_cirs[t]->all_nodes_map[first_node_name]->type << "," << second_node_name << "," << second_node->type << "," << impl_cirs[t]->all_nodes_map[second_node_name]->type << "," << total_diff_count[t][7]*needtoadd << endl;
                if(total_diff_count[t][7] < max_error && total_diff_count[t][7] > 0){
                  outMean << first_node_name << "," << first_node->type << "," << impl_cirs[t]->all_nodes_map[first_node_name]->type << "," << second_node_name << "," << second_node->type << "," << impl_cirs[t]->all_nodes_map[second_node_name]->type << "," << total_diff_count[t][7]*needtoadd << endl;
                }
                delete impl_cirs[t];
              }
              delete impl_cir;
            }
            for (int ii=0; ii<num_changes; ii++) {
              delete spec_cirs[ii];
            }
          }
          delete spec_cir;
        }
      }
      // for (int i = 0; i < impl_cirs.size(); i++) {
      //   delete impl_cirs[i];
      // }
    }break;
    case 1: {
      for (int i=0; i<first_nodes.size(); i++){
        first_node_name = first_nodes[i]->name;

        GetPossibleNodes(origin_cir, first_node_name, second_nodes, 1);

        for(int j=0; j<second_nodes.size(); j++){
          second_node_name = second_nodes[j]->name;

          spec_cir = new_spec_cir.GetDuplicate("","","");
          first_node = spec_cir->all_nodes_map[first_node_name];
          second_node = spec_cir->all_nodes_map[second_node_name];

          target_nodes.clear();
          target_nodes.push_back(first_node);
          target_nodes.push_back(second_node);
          spec_cir->Simplify(target_nodes);
          spec_cir->ResetLevels();
          spec_cir->LevelizeSortTopological(false);
          spec_cir->SetIndexes();

          if (spec_cir->all_nodes_map.find(first_node_name) != spec_cir->all_nodes_map.end() && spec_cir->all_nodes_map.find(second_node_name) != spec_cir->all_nodes_map.end()){

            needtoadd = pow(2, new_spec_cir.inputs.size() - spec_cir->inputs.size());
            max_error = ceil(pow(2, spec_cir->inputs.size()) * 0.01);

            GetWrongCircuit(spec_cir, spec_cirs, first_node_name);

            for (int ii=0; ii < num_changes; ii++) {
              impl_cir = spec_cirs[ii]->GetDuplicate("","","");
              GetWrongCircuit(impl_cir, impl_cirs, second_node_name);
              for (int jj=0; jj < num_changes; jj++) {
                sim_threads[jj] = std::thread(DoSimEqTh4, spec_cir, impl_cirs[jj], &po_diff_count[jj], &total_diff_count[jj]);
              }
              for (int jj=0; jj < num_changes; jj++) {
                sim_threads[jj].join();
              }
              for (int t=0; t < num_changes; t++) {
                cout << first_node_name << "," << first_node->type << "," << impl_cirs[t]->all_nodes_map[first_node_name]->type << "," << second_node_name << "," << second_node->type << "," << impl_cirs[t]->all_nodes_map[second_node_name]->type << "," << total_diff_count[t][7]*needtoadd << endl;
                outFile << first_node_name << "," << first_node->type << "," << impl_cirs[t]->all_nodes_map[first_node_name]->type << "," << second_node_name << "," << second_node->type << "," << impl_cirs[t]->all_nodes_map[second_node_name]->type << "," << total_diff_count[t][7]*needtoadd << endl;
                if(0 < total_diff_count[t][7] < max_error){
                  outMean << first_node_name << "," << first_node->type << "," << impl_cirs[t]->all_nodes_map[first_node_name]->type << "," << second_node_name << "," << second_node->type << "," << impl_cirs[t]->all_nodes_map[second_node_name]->type << "," << total_diff_count[t][7]*needtoadd << endl;
                }
              }
            }
          }
        }
      }
    }break;
  }
  return 0;
}

void GetWrongConnect(Circuit &new_spec_cir, string gate_name, string connection_name, int change2which){
  Node* nodechange = new_spec_cir.all_nodes_map[gate_name];
  Node* connectchange = new_spec_cir.all_nodes_map[connection_name];
  Node* previouscon = nodechange->inputs[change2which];
  connectchange->outputs.push_back(nodechange);
  nodechange->inputs[change2which] = connectchange;
  NodeVector::iterator it = previouscon->outputs.begin();
  while (it != previouscon->outputs.end()) {
    if(*it==nodechange){
      previouscon->outputs.erase(it);
      break;
    }
    ++it;
  }
  if(nodechange->level <= connectchange->level) {
    new_spec_cir.ResetLevels();
    new_spec_cir.LevelizeSortTopological(true);
    new_spec_cir.SetIndexes();
  }
}

void GetWrongConnect(Circuit *spec_cir, vector<Circuit*> &impl_cirs, vector<string>target_nodes, string con_name, string first_con_name) {
  impl_cirs.clear();
  impl_cirs.resize(target_nodes.size());

  string target_node_name;
  Circuit *impl_cir = spec_cir->GetDuplicate("","","");
  GetWrongConnect(*impl_cir, con_name, first_con_name, 0);

  for (int i = 0; i < target_nodes.size(); i++) {
    target_node_name = target_nodes[i];
    impl_cirs[i] = impl_cir->GetDuplicate("","","");
    GetWrongConnect(*impl_cirs[i], con_name, target_node_name, 1);
  }

  delete impl_cir;
}

int DoSimEqConnect(string spec_filename, int mode, int amount){
  // Circuit new_spec_cir;
  // string filename = spec_filename;
  // size_t dotp = filename.find('.');
  // ofstream outFile;
  // ofstream outText;
  // filename = filename.substr(0,dotp);
  // // if(mode==1){
  // //   outFile.open(filename+"conRes.csv", ios::out | ios::app);
  // //   outText.open(filename+"conRes.txt", ios::out | ios::app);
  // // }
  // if(mode==0){
  //   // outFile.open(filename+"test.csv", ios::out | ios::app);
  //   // outText.open(filename+"test.txt", ios::out | ios::app);
  //   outFile.open(filename+"con.csv", ios::out | ios::app);
  //   outText.open(filename+"con.txt", ios::out | ios::app);
  // }
  // if(mode==2){
  //   outFile.open(filename+"ConRandom.csv", ios::out | ios::app);
  //   outText.open(filename+"ConRandom.txt", ios::out | ios::app);
  // }
  // // if(mode==3){
  // //   outFile.open(filename+"conRadomAnaly.csv", ios::out | ios::app);
  // //   outText.open(filename+"conRandomAnaly.txt", ios::out | ios::app);
  // // }
  //
  // new_spec_cir.ReadBlif(spec_filename);
  // if (new_spec_cir.inputs.size() == 0) {
  //   MSG("Cannot simulate circuit with no primary input!");
  //   return 0;
  // }
  //
  // if (new_spec_cir.inputs.size() > 64) { // TODO: do we need to limit the inputs (because of high runtime?)
  //   MSG("Currently only simulating up to 64 inputs!");
  //   return 0;
  // }
  //
  //
  // // srand((unsigned)time( NULL));
  // // int pattern_nums[5]={1000,5000,10000,20000,50000};
  //
  // if(mode==4){
  //   NodeVector select_nodes;
  //   Circuit* spec_cir = new_spec_cir.GetDuplicate("","","");
  //   NodeVector source_nodes;
  //   GetPossibleNodes(spec_cir, "", source_nodes, 0);
  //   RanSelectNodes(source_nodes, 100, select_nodes);
  //   for (int i=0; i<select_nodes.size(); i++){
  //     cout<<select_nodes[i]->name<<endl;
  //   }
  //   return 0;
  // }
  //
  // if(amount==-1){
  //     int pattern_num = 10000;
  //     Circuit *this_spec_cir = new_spec_cir.GetDuplicate("","","");
  //     Circuit *pre_spec_cir = new_spec_cir.GetDuplicate("","","");
  //     NodeVector nodes;
  //     Node* this_node;
  //     Node* that_node;
  //     string this_node_name;
  //     string that_node_name;
  //     const int thread_num = 2;
  //     Circuit* spec_cirs[thread_num];
  //     Circuit* impl_cirs[thread_num];
  //     vector<bit64> po_diff_count[thread_num];
  //     vector<bit64> total_diff_count[thread_num];
  //     std::thread sim_threads[thread_num];
  //     GetPossibleNodes(pre_spec_cir, "", nodes, 0);
  //     NodeVector connections;
  //     NodeVector this_nodes;
  //     for (long i=0; i<nodes.size(); i++){
  //       this_node_name = nodes[i]->name;
  //       pre_spec_cir = new_spec_cir.GetDuplicate("","","");
  //       connections.clear();
  //       GetPossibleConnection(pre_spec_cir, this_node_name, connections);
  //       for (long j=0; j<connections.size(); j++){
  //         this_spec_cir = new_spec_cir.GetDuplicate("","","");
  //         that_node_name = connections[j]->name;
  //         this_node = this_spec_cir->all_nodes_map[this_node_name];
  //         that_node = this_spec_cir->all_nodes_map[that_node_name];
  //         this_nodes.clear();
  //         this_nodes.push_back(this_node);
  //         this_nodes.push_back(that_node);
  //         this_spec_cir->Simplify(this_nodes);
  //         this_spec_cir->SetIndexes();
  //         if (this_spec_cir->all_nodes_map.find(that_node_name) != this_spec_cir->all_nodes_map.end() and this_spec_cir->all_nodes_map.find(this_node_name) != this_spec_cir->all_nodes_map.end())  {
  //           spec_cirs[0] = this_spec_cir->GetDuplicate("","","");
  //           impl_cirs[0] = this_spec_cir->GetDuplicate("","","");
  //           GetWrongConnect(*impl_cirs[0], this_node_name, that_node_name, 0);
  //           if (mode==0) sim_threads[0] = std::thread(DoSimEqTh4, spec_cirs[0], impl_cirs[0], &po_diff_count[0], &total_diff_count[0]);
  //           spec_cirs[1] = this_spec_cir->GetDuplicate("","","");
  //           impl_cirs[1] = this_spec_cir->GetDuplicate("","","");
  //           GetWrongConnect(*impl_cirs[1], this_node_name, that_node_name, 1);
  //           if (mode==0) sim_threads[1] = std::thread(DoSimEqTh4, spec_cirs[1], impl_cirs[1], &po_diff_count[1], &total_diff_count[1]);
  //           if (mode==0){
  //             for(int k=0; k<2; k++){
  //               sim_threads[k].join();
  //             }
  //             for(int t=0; t<2; t++){
  //               double max_error = ceil(pow(2, impl_cirs[t]->inputs.size()) * 0.01);
  //               double needtoadd = pow(2, new_spec_cir.inputs.size() - this_spec_cir->inputs.size());
  //               if (total_diff_count[t][7] <= max_error){
  //                 outText << "----- changed node " << this_node->name << " ----"<< "connection " << that_node->name<<"---"<<"which input-----"<<t<<"----";
  //                 outText <<"wrong cases "<< total_diff_count[t][7]*needtoadd<<endl;
  //                 outFile << this_node->name<<","<<this_node->level<<","<< that_node->name<<","<<that_node->level<<","<<notfouts(that_node->outputs, this_node)<<","<<t<<","<<total_diff_count[t][7]*needtoadd <<endl;
  //                 for (int ii=0; ii<impl_cirs[t]->outputs.size(); ii++){
  //                   outText<<impl_cirs[t]->outputs[ii]->name<<"--";
  //                   cout<<endl;
  //                   outText<<endl;
  //                 }
  //               }
  //             }
  //           }
  //           if (mode==2){
  //             int h;
  //             for (int t=0; t<2; t++){
  //               cout << "------change node----- " << this_node_name << " ----connection name-----" << that_node_name << "----- "<<endl;
  //
  //               // for(int t=0; t<2; t++){
  //               //   h = RndSimEq(spec_cirs[t], impl_cirs[t], pattern_num, );
  //               //   cout << "------change node----- " << this_node_name << " ----connection name-----" << that_node_name << "-----which input---- " << t << "---wrong case----" << h << endl;
  //               //   outText << "------change node----- " << this_node_name << " ----connection name-----" << that_node_name << "-----which input---- " << t << "---wrong case----" << h << endl;
  //               //   outFile << this_node_name << "," << that_node_name << "," << t << "," << notfouts(that_node->outputs, this_node) <<","<< h <<endl;
  //               // }
  //             }
  //           }
  //         }
  //       }
  //     }
  //     return 0;
  //   }
  //
  //
  // else if(amount>0){
  //   const int threads_num = 16;
  //   Circuit* spec_cirs[threads_num];
  //   Circuit* impl_cirs[threads_num];
  //   vector<bit64> po_diff_count[threads_num];
  //   vector<bit64> total_diff_count[threads_num];
  //   NodeVector outer_nodes;
  //   NodeVector inner_nodes;
  //   std::thread sim_threads[threads_num];
  //   int inner_amount = 20;
  //   Node* outer_node = NULL;
  //   Node* inner_node = NULL;
  //   Circuit* this_spec_cir;
  //   string outer_node_name;
  //   string inner_node_name;
  //   int cur_thread_num = 0;
  //   int pattern_num = 10000;
  //   NodeVector select_inner_nodes;
  //   NodeVector source_nodes;
  //
  //
  //   Circuit* spec_cir = new_spec_cir.GetDuplicate("","","");
  //   NodeVector total_nodes;
  //   GetPossibleNodes(spec_cir, "", total_nodes, 0);
  //   NodeVector select_outer_nodes;
  //   RanSelectNodes(total_nodes, amount, select_outer_nodes);
  //   for (int i=0; i<select_outer_nodes.size(); i++){
  //     outer_node_name = select_outer_nodes[i]->name;
  //     source_nodes.clear();
  //     select_inner_nodes.clear();
  //     inner_amount = 20;
  //     GetPossibleConnection(spec_cir, outer_node_name, source_nodes);
  //     if(source_nodes.size()<inner_amount) inner_amount = source_nodes.size();
  //     RanSelectNodes(source_nodes, inner_amount, select_inner_nodes);
  //     for (int j=0; j<select_inner_nodes.size(); j++){
  //       inner_node_name = select_inner_nodes[j]->name;
  //       this_spec_cir = new_spec_cir.GetDuplicate("","","");
  //       outer_node = this_spec_cir->all_nodes_map[outer_node_name];
  //       inner_node = this_spec_cir->all_nodes_map[inner_node_name];
  //       NodeVector this_nodes;
  //       this_nodes.push_back(outer_node);
  //       this_nodes.push_back(inner_node);
  //       this_spec_cir->Simplify(this_nodes);
  //       this_spec_cir->SetIndexes();
  //       if (this_spec_cir->all_nodes_map.find(outer_node_name) != this_spec_cir->all_nodes_map.end() and this_spec_cir->all_nodes_map.find(inner_node_name) != this_spec_cir->all_nodes_map.end())  {
  //         spec_cirs[cur_thread_num] = this_spec_cir->GetDuplicate("","","");
  //         impl_cirs[cur_thread_num] = this_spec_cir->GetDuplicate("","","");
  //         GetWrongConnect(*impl_cirs[cur_thread_num], outer_node_name, inner_node_name, 0);
  //         if(mode==0) sim_threads[cur_thread_num] = std::thread(DoSimEqTh4, spec_cirs[cur_thread_num], impl_cirs[cur_thread_num], &po_diff_count[cur_thread_num], &total_diff_count[cur_thread_num]);
  //         outer_nodes.push_back(outer_node);
  //         inner_nodes.push_back(inner_node);
  //         cur_thread_num++;
  //         spec_cirs[cur_thread_num] = this_spec_cir->GetDuplicate("","","");
  //         impl_cirs[cur_thread_num] = this_spec_cir->GetDuplicate("","","");
  //         GetWrongConnect(*impl_cirs[cur_thread_num], outer_node_name, inner_node_name, 1);
  //         if(mode==0) sim_threads[cur_thread_num] = std::thread(DoSimEqTh4, spec_cirs[cur_thread_num], impl_cirs[cur_thread_num], &po_diff_count[cur_thread_num], &total_diff_count[cur_thread_num]);
  //         cur_thread_num++;
  //       }
  //       if (mode==2){
  //         int h;
  //         cout << "------change node----- "<< outer_node_name << " ----connection name-----"<<inner_node_name<<"----- "<<endl;
  //         for(int t=0; t<2; t++){
  //           // h=RndSimEq(spec_cirs[t], impl_cirs[t], pattern_num);
  //           cout << "------change node----- " << outer_node_name << " ----connection name-----" << inner_node_name << "-----which input---- " << t << "---wrong case----" << h << endl;
  //           outText << "------change node----- " << outer_node_name << " ----connection name-----" << inner_node_name << "-----which input---- " << t << "---wrong case----" << h << endl;
  //           outFile << outer_node_name << "," << inner_node_name << "," << t << "," << notfouts(inner_node->outputs, outer_node) <<","<< h <<endl;
  //         }
  //         cur_thread_num = 0;
  //       }
  //       if (mode == 0 and cur_thread_num == threads_num){
  //         for (int t=0; t<cur_thread_num; t++){
  //           sim_threads[t].join();
  //         }
  //         for (int t=0; t<cur_thread_num; t++){
  //           double needtoadd = pow(2, new_spec_cir.inputs.size() - spec_cirs[t]->inputs.size());
  //           double max_error = ceil(pow(2, spec_cirs[t]->inputs.size()) * 0.01);
  //           // if (total_diff_count[t][7] <= max_error){
  //           //   outFile << outer_nodes[floor(t/2)] << "," << inner_nodes[floor(t/2)] << "," << t%2 << "," << notfouts(inner_nodes[floor(t/2)]->outputs, outer_nodes[floor(t/2)]) <<","<< total_diff_count[t][7] <<endl;
  //           // }
  //           cout << " change node " << outer_nodes[floor(t/2)]->name << " connection name " << inner_nodes[floor(t/2)]->name << " which input " << t%2 << " wrong case " << total_diff_count[t][7]*needtoadd << endl;
  //           outFile << outer_nodes[floor(t/2)]->name << "," << inner_nodes[floor(t/2)]->name << "," << t%2 << "," << notfouts(inner_nodes[floor(t/2)]->outputs, outer_nodes[floor(t/2)]) <<","<< total_diff_count[t][7]*needtoadd <<endl;
  //         }
  //         cur_thread_num = 0;
  //         inner_nodes.clear();
  //         outer_nodes.clear();
  //       }
  //     }
  //   }
  //   return 0;
  // }
  //
  // outFile.close();
  // outText.close();
  return 0;
}

int DosimEqMulCon(Circuit &new_spec_cir, string spec_filename, NodeVector gate_nodes, int mode){
  ofstream outFile;
  ofstream outMean;
  string filename = spec_filename;
  size_t dotp = filename.find('.');
  filename = filename.substr(0,dotp);

  outFile.open(filename+"mulCon.csv",  ios::out);
  outMean.open(filename+"wished.csv",  ios::out | ios::app);

  const int thread_num = 15;


  std::vector<Circuit*> impl_cirs;
  Circuit* spec_cir = new_spec_cir.GetDuplicate("", "", "");
  Circuit* impl_cir;
  std::thread sim_threads[thread_num];
  vector<bit64> po_diff_count[thread_num];
  vector<bit64> total_diff_count[thread_num];

  vector<string> connections;
  vector<string> nv2;
  vector<string> nv3;

  Node* first_node = new Node;
  Node* second_node = new Node;
  Node* gate_node = new Node;
  string first_node_name;
  string second_node_name;
  string gate_node_name;

  NodeVector connection_nodes;
  NodeVector simplify_nodes;

  double needtoadd = 1;
  double max_error = 0;
  long int cn_size;
  long int thread_size;

  switch(mode) {
    case 0:{
      for (int i = 0; i < gate_nodes.size(); i++){

        spec_cir = new_spec_cir.GetDuplicate("", "", "");
        gate_node_name = gate_nodes[i]->name;

        GetPossibleConnection(spec_cir, gate_node_name, connection_nodes);
        cn_size = connection_nodes.size();

        for (int j = 0; j < cn_size; j++){
          first_node_name = connection_nodes[j]->name;

          for (int jj=0; jj<cn_size; jj++){
            if (j == jj) continue;

            second_node_name = connection_nodes[jj]->name;
            impl_cir = new_spec_cir.GetDuplicate("","","");
            gate_node = impl_cir->all_nodes_map[gate_node_name];
            first_node = impl_cir->all_nodes_map[first_node_name];
            second_node = impl_cir->all_nodes_map[second_node_name];
            simplify_nodes.clear();
            simplify_nodes.push_back(gate_node);
            simplify_nodes.push_back(first_node);
            simplify_nodes.push_back(second_node);
            impl_cir->Simplify(simplify_nodes);
            impl_cir->SetIndexes();
            needtoadd = pow(2, new_spec_cir.inputs.size() - impl_cir->inputs.size());
            max_error = ceil(pow(2, impl_cir->inputs.size())*0.02);

            if (impl_cir->all_nodes_map.find(gate_node_name) != impl_cir->all_nodes_map.end() && impl_cir->all_nodes_map.find(first_node_name) != impl_cir->all_nodes_map.end() && impl_cir->all_nodes_map.find(second_node_name) != impl_cir->all_nodes_map.end()){
              connections.push_back(second_node_name);

              if (jj == cn_size - 1 || connections.size() >= thread_num) {

                GetWrongConnect(impl_cir, impl_cirs, connections, gate_node_name, first_node_name);
                thread_size = connections.size();

                cout << thread_size << endl;

                for (int k = 0; k < thread_size; k++) {
                  // cout << impl_cirs[k]->all_nodes_map[gate_node_name]->inputs[0]->name << "=====" << impl_cirs[k]->all_nodes_map[gate_node_name]->inputs[1]->name  << "====" << k << "====" << thread_size << endl;
                  sim_threads[k] = std::thread(DoSimEqTh4, impl_cir, impl_cirs[k], &po_diff_count[k], &total_diff_count[k]);
                }
                for (int k = 0; k < thread_size; k++) {
                  sim_threads[k].join();
                }
                for (int k = 0; k < thread_size; k++) {
                  cout << gate_node_name << "," << first_node_name << "," << connections[k] << "," << total_diff_count[k][7]*needtoadd << endl;
                  outFile << gate_node_name << "," << first_node_name << "," << connections[k] << "," << total_diff_count[k][7]*needtoadd << endl;
                  if (total_diff_count[k][7] < max_error) {
                    outMean << gate_node_name << "," << first_node_name << "," << connections[k] << "," << total_diff_count[k][7]*needtoadd << endl;
                  }
                }
                connections.clear();
              }
            }
          }
        }
      }
      delete impl_cir;
      delete spec_cir;
      for (int i = 0; i < impl_cirs.size(); i++) {
        delete impl_cirs[i];
      }
    }break;
    // case 1: {
    //   Circuit *temp_cir;
    //   NodeVector type_nodes;
    //   string type_node_name;
    //
    //   for (int i = 0; i < gate_nodes.size(); i++){
    //
    //     spec_cir = new_spec_cir.GetDuplicate("", "", "");
    //     gate_node_name = gate_nodes[i]->name;
    //
    //     GetPossibleConnection(spec_cir, gate_node_name, connection_nodes);
    //     cn_size = connection_nodes.size();
    //
    //     for (int j = 0; j < cn_size; j++){
    //       first_node_name = connection_nodes[j]->name;
    //
    //       for (int jj=0; jj<cn_size; jj++){
    //         if (j == jj) continue;
    //         second_node_name = connection_nodes[jj]->name;
    //         impl_cir = new_spec_cir.GetDuplicate("","","");
    //         gate_node = impl_cir->all_nodes_map[gate_node_name];
    //         first_node = impl_cir->all_nodes_map[first_node_name];
    //         second_node = impl_cir->all_nodes_map[second_node_name];
    //
    //         simplify_nodes.clear();
    //         simplify_nodes.push_back(gate_node);
    //         simplify_nodes.push_back(first_node);
    //         simplify_nodes.push_back(second_node);
    //         impl_cir->Simplify(simplify_nodes);
    //         impl_cir->SetIndexes();
    //         needtoadd = pow(2, new_spec_cir.inputs.size() - impl_cir->inputs.size());
    //         max_error = ceil(pow(2, impl_cir->inputs.size())*0.02);
    //
    //         if (impl_cir->all_nodes_map.find(gate_node_name) != impl_cir->all_nodes_map.end() && impl_cir->all_nodes_map.find(first_node_name) != impl_cir->all_nodes_map.end() && impl_cir->all_nodes_map.find(second_node_name) != impl_cir->all_nodes_map.end()){
    //           temp_cir = impl_cir->GetDuplicate("","","");
    //           GetWrongConnect(*temp_cir, gate_node_name, first_node_name, 0);
    //           GetWrongConnect(*temp_cir, gate_node_name, second_node_name, 1);
    //           GetPossibleNodes(temp_cir, "", type_nodes, 0);
    //
    //           for (int kk = 0; kk < type_nodes.size(); kk++) {
    //             type_node_name = type_nodes[kk]->name;
    //             GetWrongCircuit(temp_cir, impl_cirs, type_node_name);
    //
    //
    //             for (int k = 0; k < thread_num; k++) {
    //               sim_threads[k] = std::thread(DoSimEqTh4, impl_cir, impl_cirs[k], &po_diff_count[k], &total_diff_count[k]);
    //             }
    //
    //             for (int k = 0; k < thread_num; k++) {
    //               sim_threads[k].join();
    //             }
    //             for (int k = 0; k < thread_num; k++) {
    //               cout << gate_node_name << "," << first_node_name << "," << second_node_name << "," << type_node_name << "," << impl_cirs[k]->all_nodes_map[type_node_name]->type<< "," << total_diff_count[k][7]*needtoadd << endl;
    //               outFile << gate_node_name << "," << first_node_name << "," << second_node_name << "," << type_node_name << "," << impl_cirs[k]->all_nodes_map[type_node_name]->type<< "," << total_diff_count[k][7]*needtoadd << endl;
    //               if (total_diff_count[k][7] < max_error) {
    //                 outMean << gate_node_name << "," << first_node_name << "," << second_node_name << "," << type_node_name << "," << impl_cirs[k]->all_nodes_map[type_node_name]->type<< "," << total_diff_count[k][7]*needtoadd << endl;
    //               }
    //             }
    //           }
    //         }
    //       }
    //     }
    //   }
    //   delete impl_cir;
    //   delete spec_cir;
    //   delete temp_cir;
    //   for (int i = 0; i < impl_cirs.size(); i++) {
    //     delete impl_cirs[i];
    //   }
    //
    // }break;
  }

  return 0;
}

int DoMixedChange(Circuit &new_spec_cir, string spec_filename, int amount, int mode){
  ofstream outFile;
  ofstream outMean;
  string filename = spec_filename;
  size_t dotp = filename.find('.');
  filename = filename.substr(0,dotp);
  outMean.open(filename+"wished.csv", ios::out | ios::app);

  srand((unsigned)time( NULL));

  const int thread_num = 15;

  // amount = 10;

  // mode: 0, change type and connection on one node
  // mode: 1, change type and connection on two nodes


  if (mode==0){

    /*
      gate_nodes: all nodes of the circuits (2 input);
      connection_nodes: for one node, those connections can be re-connected to
      thread_num: threads used in parallel
    */
    outFile.open(filename+"mixed.csv",  ios::out);

    double max_error = 0;
    double needtoadd = 1;

    Circuit* impl_cir = new_spec_cir.GetDuplicate("","","");
    Circuit* spec_cir = new_spec_cir.GetDuplicate("","","");
    std::vector<Circuit*> impl_cirs;
    vector<bit64> po_diff_count[thread_num];
    vector<bit64> total_diff_count[thread_num];
    std::thread sim_threads[thread_num];


    NodeVector simplify_nodes;
    NodeVector gate_nodes;
    NodeVector connection_nodes;

    string gate_node_name;
    string connection_name;
    string inner_node_name;

    Node* connection_node = new Node;
    Node* gate_node = new Node;

    long int node_size = gate_nodes.size();
    long int connection_node_size;

    GetPossibleNodes(spec_cir, "", gate_nodes, 0);

    for (int i = 0; i < node_size; i++){
      gate_node_name = gate_nodes[i]->name;
      spec_cir = new_spec_cir.GetDuplicate("","","");
      GetPossibleConnection(spec_cir, gate_node_name, connection_nodes);
      connection_node_size = connection_nodes.size();

      for (int j=0; j<connection_node_size; j++){
        impl_cir = new_spec_cir.GetDuplicate("","","");
        connection_name = connection_nodes[j]->name;
        connection_node = impl_cir->all_nodes_map[connection_name];
        gate_node = impl_cir->all_nodes_map[gate_node_name];

        simplify_nodes.clear();
        simplify_nodes.push_back(gate_node);
        simplify_nodes.push_back(connection_node);
        impl_cir->Simplify(simplify_nodes);
        impl_cir->ResetLevels();
        impl_cir->LevelizeSortTopological(false);
        impl_cir->SetIndexes();


        if (impl_cir->all_nodes_map.find(gate_node_name) != impl_cir->all_nodes_map.end() && impl_cir->all_nodes_map.find(connection_name) != impl_cir->all_nodes_map.end()) {

          GetWrongCircuit(impl_cir, impl_cirs, gate_node_name);
          needtoadd = pow(2, new_spec_cir.inputs.size() - impl_cir->inputs.size());
          max_error = ceil(pow(2, impl_cir->inputs.size()) * 0.02);

          for(int ii=0; ii<impl_cirs.size(); ii++) {
            GetWrongConnect(*impl_cirs[ii], gate_node_name, connection_name, rand()%1);
            sim_threads[ii] = std::thread(DoSimEqTh4, impl_cir, impl_cirs[ii], &po_diff_count[ii], &total_diff_count[ii]);
          }

          for(int jj=0; jj<thread_num; jj++) {
            sim_threads[jj].join();
          }

          for (int k=0; k<thread_num; k++){
            cout << "connection from " << connection_name <<" to gate " << gate_node_name << " gate type change to " << impl_cirs[k]->all_nodes_map[gate_node_name]->type << " wrong is " << total_diff_count[k][7] << endl;
            outFile << connection_name << "," << gate_node_name << "," << impl_cirs[k]->all_nodes_map[gate_node_name]->type  << "," <<total_diff_count[k][7] * needtoadd << endl;
            if (total_diff_count[k][7] <= max_error) {
              outMean << connection_name << "," << gate_node_name << "," << impl_cirs[k]->all_nodes_map[gate_node_name]->type  << "," <<total_diff_count[k][7] * needtoadd << endl;
            }
          }
        }
      }
    }
    delete impl_cir;
    delete spec_cir;
    for(int jj=0; jj<thread_num; jj++) delete impl_cirs[jj];
  }

  if (mode==1) {

    /*
      gate_nodes: change its type
      con_nodes: change its connection
      connection_nodes: reconnect it to con_nodes
      thread_num: threads
    */
    outFile.open(filename+"mixed.csv",  ios::out);

    Circuit* spec_cir = new_spec_cir.GetDuplicate("","","");
    Circuit* impl_cir;
    std::vector<Circuit*> impl_cirs;
    std::vector<bit64> po_diff_count[thread_num];
    std::vector<bit64> total_diff_count[thread_num];
    std::thread sim_threads[thread_num];

    NodeVector gate_nodes;
    NodeVector con_nodes;
    NodeVector connection_nodes;
    NodeVector simplify_nodes;

    Node* gate_node = new Node;
    Node* connection_node = new Node;
    Node* con_node = new Node;

    string gate_node_name;
    string connection_node_name;
    string con_node_name;

    GetPossibleNodes(spec_cir, "", gate_nodes, 0);

    double needtoadd =  1;
    double max_error;
    long int node_size = gate_nodes.size();
    long int connection_node_size;

    for (int i = 0; i < node_size; i++){
      gate_node_name = gate_nodes[i]->name;
      cout << endl << gate_node_name << endl;

      for (int j = 0; j < node_size; j++){

        con_node_name = gate_nodes[j]->name;

        cout << endl << gate_node_name << "," << con_node_name  << endl;
        spec_cir = new_spec_cir.GetDuplicate("", "", "");

        GetPossibleConnection(spec_cir, con_node_name, connection_nodes);

        connection_node_size = connection_nodes.size();
        for (int p = 0; p < connection_node_size; p++) {
          connection_node_name = connection_nodes[p]->name;


          impl_cir = new_spec_cir.GetDuplicate("","","");
          gate_node = impl_cir->all_nodes_map[gate_node_name];
          con_node = impl_cir->all_nodes_map[con_node_name];
          connection_node = impl_cir->all_nodes_map[connection_node_name];

          simplify_nodes.clear();
          simplify_nodes.push_back(gate_node);
          simplify_nodes.push_back(con_node);
          simplify_nodes.push_back(connection_node);
          impl_cir->Simplify();
          impl_cir->SetIndexes();

          needtoadd = pow(2, new_spec_cir.inputs.size() - impl_cir->inputs.size());
          max_error = ceil(pow(2, impl_cir->inputs.size()) * 0.02);

          GetWrongCircuit(impl_cir, impl_cirs, gate_node_name);

          for (int ii = 0; ii < thread_num; ii++) {
            GetWrongConnect(*impl_cirs[ii], con_node_name, connection_node_name, rand()%1);
            sim_threads[ii] = std::thread(DoSimEqTh4, impl_cir, impl_cirs[ii], &po_diff_count[ii], &total_diff_count[ii]);
          }

          for (int k = 0; k < thread_num; k++){
            sim_threads[k].join();
          }

          for (int k = 0; k < thread_num; k++){
            cout << gate_node_name << "," << gate_node->type << "," << impl_cirs[k]->all_nodes_map[gate_node_name]->type << "," << connection_node_name <<"," << con_node_name << "," << total_diff_count[k][7] * needtoadd << endl;
            outFile << gate_node_name << "," << gate_node->type << "," << impl_cirs[k]->all_nodes_map[gate_node_name]->type << "," << connection_node_name << "," << con_node_name <<"," << total_diff_count[k][7] * needtoadd << endl;
            if ( total_diff_count[k][7] < max_error ) {
              outMean << gate_node_name << "," << impl_cirs[k]->all_nodes_map[gate_node_name]->type << "," << connection_node_name << "," << con_node_name <<"," << total_diff_count[k][7] * needtoadd << endl;
            }
          }
          delete impl_cir;
          for (int k = 0; k < thread_num; k++){
            delete impl_cirs[k];
          }
        }
        delete spec_cir;
      }
      // delete spec_cir;
      // delete impl_cir;
      // for (int k = 0; k < impl_cirs.size(); k++){
      //   delete impl_cirs[k];
      // }
    }
    cout << endl;


  }

  outFile.close();
  outMean.close();
  return 0;
}

int DoHueristic(Circuit &new_spec_cir, NodeVector target_nodes, int mode=1) {

  ofstream outFile;
  outFile.open("static.csv", ios::out | ios::app);


  long int node_size = target_nodes.size();
  const int nums_start = 2;
  const int nums_end = 6;
  const int nums = 8;
  std::thread sim_threads[nums];
  Circuit *spec_cir;
  vector<Circuit*> impl_cirs;
  string target_node_name;
  vector<bit64> po_diff_count[nums];
  vector<bit64> total_diff_count[nums];
  NodeVector connections;
  Node* cur_node;
  long int con_size;
  string con_name;
  int real_thread;
  int j = 0;
  Circuit *c1 = new_spec_cir.GetDuplicate("","","");
// GetPossibleConnection(Circuit *spec_cir, string current_node_name, NodeVector &connections
  switch (mode) {
    case 1: {
      for (long int i=0; i < node_size; i++) {
        target_node_name = target_nodes[i]->name;
        c1 = new_spec_cir.GetDuplicate("","","");
        connections.clear();
        GetPossibleConnection(c1, target_node_name, connections);
        con_size = connections.size();
        for (long k=0; k<con_size-1; k+=2) {
          if (k+1>=con_size) continue;
          spec_cir = new_spec_cir.GetDuplicate("","","");
          con_name = connections[k]->name;
          impl_cirs.clear();
          impl_cirs.resize(nums);
          j = 0;
          impl_cirs[j] = spec_cir->GetDuplicate("","","");
          GetWrongCircuit(*impl_cirs[j], target_node_name, 2);
          GetWrongConnect(*impl_cirs[j],target_node_name, con_name, 0);
          sim_threads[j] = std::thread(DoSimEqTh4, spec_cir, impl_cirs[j], &po_diff_count[j], &total_diff_count[j]);
          j = 1;
          impl_cirs[j] = spec_cir->GetDuplicate("","","");
          GetWrongCircuit(*impl_cirs[j], target_node_name, 4);
          GetWrongConnect(*impl_cirs[j],target_node_name, con_name, 0);
          sim_threads[j] = std::thread(DoSimEqTh4, spec_cir, impl_cirs[j], &po_diff_count[j], &total_diff_count[j]);
          con_name = connections[k+1]->name;
          j = 2;
          impl_cirs[j] = spec_cir->GetDuplicate("","","");
          GetWrongCircuit(*impl_cirs[j], target_node_name, 2);
          GetWrongConnect(*impl_cirs[j],target_node_name, con_name, 0);
          sim_threads[j] = std::thread(DoSimEqTh4, spec_cir, impl_cirs[j], &po_diff_count[j], &total_diff_count[j]);
          j = 3;
          impl_cirs[j] = spec_cir->GetDuplicate("","","");
          GetWrongCircuit(*impl_cirs[j], target_node_name, 4);
          GetWrongConnect(*impl_cirs[j],target_node_name, con_name, 0);
          sim_threads[j] = std::thread(DoSimEqTh4, spec_cir, impl_cirs[j], &po_diff_count[j], &total_diff_count[j]);
          real_thread = 4;
          if (spec_cir->all_nodes_map[target_node_name]->type == NODE_NOR) {
            real_thread = 8;
            j = 4;
            con_name = connections[k]->name;
            impl_cirs[j] = spec_cir->GetDuplicate("","","");
            GetWrongConnect(*impl_cirs[j],target_node_name, con_name, 0);
            GetWrongCircuit(*impl_cirs[j], target_node_name, 6);
            sim_threads[j] = std::thread(DoSimEqTh4, spec_cir, impl_cirs[j], &po_diff_count[j], &total_diff_count[j]);
            j = 5;
            impl_cirs[j] = spec_cir->GetDuplicate("","","");
            GetWrongConnect(*impl_cirs[j],target_node_name, con_name, 1);
            GetWrongCircuit(*impl_cirs[j], target_node_name, 6);
            sim_threads[j] = std::thread(DoSimEqTh4, spec_cir, impl_cirs[j], &po_diff_count[j], &total_diff_count[j]);
            j = 6;
            con_name = connections[k+1]->name;
            impl_cirs[j] = spec_cir->GetDuplicate("","","");
            GetWrongConnect(*impl_cirs[j],target_node_name, con_name, 0);
            GetWrongCircuit(*impl_cirs[j], target_node_name, 6);
            sim_threads[j] = std::thread(DoSimEqTh4, spec_cir, impl_cirs[j], &po_diff_count[j], &total_diff_count[j]);
            j = 7;
            impl_cirs[j] = spec_cir->GetDuplicate("","","");
            GetWrongConnect(*impl_cirs[j],target_node_name, con_name, 1);
            GetWrongCircuit(*impl_cirs[j], target_node_name, 6);
            sim_threads[j] = std::thread(DoSimEqTh4, spec_cir, impl_cirs[j], &po_diff_count[j], &total_diff_count[j]);

          }

          for (j=0; j<real_thread; j++) {
            sim_threads[j].join();
          }

          for (j=0; j<real_thread; j++) {
           cout  << target_node_name << "," << impl_cirs[j]->all_nodes_map[target_node_name]->inputs[0]->name << "," << impl_cirs[j]->all_nodes_map[target_node_name]->type << ","<< total_diff_count[j][7] << endl;
           outFile  << target_node_name << "," << impl_cirs[j]->all_nodes_map[target_node_name]->inputs[0]->name  << "," << impl_cirs[j]->all_nodes_map[target_node_name]->type << "," << total_diff_count[j][7] << endl;
          }

         cout << endl;

         delete spec_cir;
         for (j=0; j<real_thread; j++) {
           delete impl_cirs[j];
         }
       }
       delete c1;
      }
    }break;
    case 2: {

      for (long int i = 0; i < node_size; i++) {
        real_thread = 6;
        target_node_name = target_nodes[i]->name;
        spec_cir = new_spec_cir.GetDuplicate("","","");
        impl_cirs.clear();
        impl_cirs.resize(real_thread);
        for (int j=0; j<real_thread; j++) {
          impl_cirs[j] = spec_cir->GetDuplicate("","","");
          GetWrongCircuit(*impl_cirs[j], target_node_name, j);
          sim_threads[j] = std::thread(DoSimEqTh4, spec_cir, impl_cirs[j], &po_diff_count[j], &total_diff_count[j]);
        }
        for (int j=0; j<real_thread; j++) {
          sim_threads[j].join();
        }
        // outFile << target_node_name;
        for (int j=0; j<real_thread; j++) {
          cout  << target_node_name << "," << impl_cirs[j]->all_nodes_map[target_node_name]->type << "," << total_diff_count[j][7] << endl;
          outFile  << target_node_name << "," << impl_cirs[j]->all_nodes_map[target_node_name]->type << "," << total_diff_count[j][7] << endl;
        }
        // outFile << endl;
        // for (int j=0; j<nums; j++) {
        //   for (int jj=1; jj<8; jj++) {
        //     outFile << total_diff_count[j][jj] - total_diff_count[j][jj-1] << ",";
        //   }
        //   outFile << endl;
        // }
        cout << endl;
        // outFile << endl;
        delete spec_cir;
        for (int j=0; j<real_thread; j++) {
          delete impl_cirs[j];
        }
      }

    }break;
    default: break;
  }


  outFile.close();
  return 0;
}

int DoSimEqLim(string spec_filename, int mode ,int amount, int submode=0){
  Circuit new_spec_cir;
  new_spec_cir.ReadBlif(spec_filename);

  if (new_spec_cir.inputs.size() == 0) {
    MSG("Cannot simulate circuit with no primary input!");
    return 0;
  }

  if (new_spec_cir.inputs.size() > 64) { // TODO: do we need to limit the inputs (because of high runtime?)
    MSG("Currently only simulating up to 64 inputs!");
    return 0;
  }


  NodeVector first_nodes;
  Circuit* spec_cir = new_spec_cir.GetDuplicate("", "", "");

  if(amount!=-1){
    NodeVector all_nodes;
    GetPossibleNodes(spec_cir, "", all_nodes, 0);
    RanSelectNodes(all_nodes, amount, first_nodes);
  }
  else{
    GetPossibleNodes(spec_cir, "", first_nodes, 0);
  }

  switch (mode) {
    case 0:
      DoSimEqMod(new_spec_cir, first_nodes, spec_filename);
      break;
    case 1:
      DoSimEqModRes(new_spec_cir, first_nodes, spec_filename, 1);
      break;
    case 2:
      DoSimEqModRes(new_spec_cir, first_nodes, spec_filename, 2);
      break;
    case 3:
      DoSimEqModRes(new_spec_cir, first_nodes, spec_filename, 2);
      break;
    case 4:
      DoSimEqMulGate(new_spec_cir, spec_filename, first_nodes, 0);
      break;
    case 5:
      DoMixedChange(new_spec_cir, spec_filename, amount, 1);
      break;
    case 6:
      DosimEqMulCon(new_spec_cir, spec_filename, first_nodes, submode);
      break;
    case 7:
      DoHueristic(new_spec_cir, first_nodes, 2);
      break;
    case 8:
      DoHueristic(new_spec_cir, first_nodes, 1);
      break;
    default:
      break;
  }
  delete spec_cir;
  return 0;
}

int DoConnectionTwo(Circuit &spec_cir_1, Circuit &spec_cir_2, vector<bit64> &total_diff_count_1, vector<bit64> &total_diff_count_2, string gate_name, string connection_name) {
  Circuit *impl_cir_1 = spec_cir_1.GetDuplicate("", "", "");
  Circuit *impl_cir_2 = spec_cir_2.GetDuplicate("", "", "");
  GetWrongConnect(*impl_cir_1, gate_name, connection_name, 1);
  GetWrongConnect(*impl_cir_2, gate_name, connection_name, 1);
  vector<bit64> unused;
  DoSimEq(spec_cir_1, *impl_cir_1, total_diff_count_1, unused);
  DoSimEq(spec_cir_2, *impl_cir_2, total_diff_count_2, unused);
  return 0;
}

void doConTwoTh4(Circuit *spec_cir_1, Circuit *spec_cir_2, vector<bit64> *toatal_diff_count_1, vector<bit64> *toatal_diff_count_2, string gate_name, string connection_name) {
  DoConnectionTwo(*spec_cir_1, *spec_cir_2, *toatal_diff_count_1, *toatal_diff_count_2, gate_name, connection_name);
}

int issamenode(Node* n1, Node* n2){
  if ((n1->inputs[0]->name == n2->inputs[0]->name && n1->inputs[1]->name == n2->inputs[1]->name) || (n1->inputs[1]->name == n2->inputs[0]->name && n1->inputs[0]->name == n2->inputs[1]->name)){
    return 1;
  }
  return 0;
}

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, std::back_inserter(elems));
    return elems;
}

/*
dosimeqmodtwo: do simeq on gates to two circuits
*/
int unit4mod(Circuit &new_spec_cir, string target_node_name, vector<bit64> (&po_diff_count)[15], vector<bit64> (&total_diff_count)[15]){

  const int num_changes = 15;
  Circuit* impl_cirs[num_changes];
  Circuit* spec_cir;
  std::thread sim_threads[num_changes];



  for (int i = 0; i < 15; i++) {
    po_diff_count[i].clear();
    po_diff_count[i].resize(new_spec_cir.outputs.size(), 0);
    total_diff_count[i].clear();
    total_diff_count[i].resize(8,0);
  }

  Node* target_node = new Node;

  int cur_num = 0;
  double needtoadd;

  spec_cir = new_spec_cir.GetDuplicate("","","");
  target_node = spec_cir->all_nodes_map[target_node_name];


  if (spec_cir->inputs.size() < new_spec_cir.inputs.size())
    needtoadd = pow(2, (new_spec_cir.inputs.size() - spec_cir->inputs.size()));
  else needtoadd = 1;

  cur_num=0;


  for (int i=0; i<num_changes+1; i++){
    if (i-1 == target_node->type)  continue;
    impl_cirs[cur_num] = spec_cir->GetDuplicate("","","");
    GetWrongCircuit(*impl_cirs[cur_num], target_node_name, i);
    sim_threads[cur_num] = std::thread(DoSimEqTh4, spec_cir, impl_cirs[cur_num], &po_diff_count[cur_num], &total_diff_count[cur_num]);
    cur_num++;
  }


  for (int t = 0; t < cur_num; t++){
    sim_threads[t].join();
  }


  for (int i = 0; i < cur_num; i++) {
    for (int j=0; j < 8; j++) {
      po_diff_count[i][j] *= needtoadd;
      total_diff_count[i][j] *= needtoadd;
    }
    delete impl_cirs[i];
  }
  delete spec_cir;
  return 0;
}

void DoSimEqModTwo(string spec_filename_1, string spec_filename_2, Circuit &new_spec_cir_1, Circuit &new_spec_cir_2, int mode, int amount) {

  NodeVector target_nodes;
  Circuit* spec_cir;

  if (amount == -1){
    if(new_spec_cir_1.all_nodes.size() < new_spec_cir_2.all_nodes.size()){
      spec_cir = new_spec_cir_1.GetDuplicate("","","");
      GetPossibleNodes(spec_cir, "", target_nodes, 0);
    }else{
      spec_cir = new_spec_cir_2.GetDuplicate("","","");
      GetPossibleNodes(spec_cir, "", target_nodes, 0);
    }
  }
  else{
    if(new_spec_cir_1.all_nodes.size()<new_spec_cir_2.all_nodes.size()){
      spec_cir = new_spec_cir_1.GetDuplicate("","","");
      NodeVector all_nodes;
      GetPossibleNodes(spec_cir, "", all_nodes, 0);
      RanSelectNodes(all_nodes, amount, target_nodes);
    }else{
      spec_cir = new_spec_cir_2.GetDuplicate("","","");
      NodeVector all_nodes;
      GetPossibleNodes(spec_cir, "", all_nodes, 0);
      RanSelectNodes(all_nodes, amount, target_nodes);
    }
  }

  if (mode==0) {

    ofstream outText;
    ofstream outDetail;
    string filename = spec_filename_1;
    size_t dotp = filename.find('.');
    filename = filename.substr(0,dotp);

    outText.open(filename+"diff.csv", ios::out | ios::app);
    outDetail.open(filename+"detail.csv", ios::out | ios::app);

    vector<bit64> po_diff_count;
    vector<bit64> total_diff_count;
    DoSimEq(new_spec_cir_1, new_spec_cir_2, po_diff_count, total_diff_count);

    outText << total_diff_count[7] << ",";
    outDetail << spec_filename_2 << "\t" << total_diff_count[7] << endl;

    for (int i=0; i< po_diff_count.size(); i++) {
      if (po_diff_count[i]!=0)
        outDetail << "port " << i << "\t" << po_diff_count[i] << endl;
    }
    outDetail << "detail" << endl;

    const int num_changes = 15;

    vector<bit64> po_diff_count_1[num_changes];
    vector<bit64> total_diff_count_1[num_changes];
    vector<bit64> po_diff_count_2[num_changes];
    vector<bit64> total_diff_count_2[num_changes];

    string target_node_name;

    Circuit *spec_cir_1;
    Circuit *spec_cir_2;

    NodeVector nv1;
    Node* n1 = NULL;

    int b = 0;
    int count = 0;
    long int node_size = target_nodes.size();
    double validNode = 0;

    for (int m = 0; m < node_size; m++){

      target_node_name = target_nodes[m]->name;

      spec_cir_1 = new_spec_cir_1.GetDuplicate("","","");
      spec_cir_2 = new_spec_cir_2.GetDuplicate("","","");

      nv1.clear();
      n1 = spec_cir_1->all_nodes_map[target_node_name];
      nv1.push_back(n1);
      spec_cir_1->Simplify(nv1);
      spec_cir_1->SetIndexes();

      nv1.clear();
      n1 = spec_cir_2->all_nodes_map[target_node_name];
      nv1.push_back(n1);
      spec_cir_2->Simplify(nv1);
      spec_cir_2->SetIndexes();

      if (spec_cir_1->all_nodes_map.find(target_node_name) != spec_cir_1->all_nodes_map.end() && spec_cir_2->all_nodes_map.find(target_node_name) != spec_cir_2->all_nodes_map.end())
      {
        unit4mod(*spec_cir_1, target_node_name, po_diff_count_1, total_diff_count_1);
        unit4mod(*spec_cir_2, target_node_name, po_diff_count_2, total_diff_count_2);

        cout << endl << "Summary of final result:" << endl;

        for (int t = 0; t < num_changes; t++) {
          cout <<  target_node_name << "," << total_diff_count_1[t][7] << "," << total_diff_count_2[t][7] << endl;
          if (total_diff_count_1[t][7]!=total_diff_count_2[t][7]) {
            b = 0;
          }
          else b = 1;
          if (b==0) {
            count += 1;
            outDetail << target_node_name << "\t" << total_diff_count_1[t][7] << "\t" << total_diff_count_2[t][7] << endl;
          }
        }
        validNode ++;
      }
      delete spec_cir_1;
      delete spec_cir_2;
    }

    cout << "all have " << count <<endl;
    outText << count << "," << validNode * num_changes << endl;

    outText.close();
    outDetail.close();

  }
}

int DoSimEqModTwo(string spec_filename_1, string spec_filename_2, int mode, int amount){

  Circuit new_spec_cir_1, new_spec_cir_2;
  new_spec_cir_1.ReadBlif(spec_filename_1);
  new_spec_cir_2.ReadBlif(spec_filename_2);

  mode = 0;


  if (new_spec_cir_1.inputs.size() == 0 || new_spec_cir_2.inputs.size() == 0 ) {
    MSG("Cannot simulate circuit with no primary input!");
    return 0;
  }

  if (new_spec_cir_1.inputs.size() > 64 || new_spec_cir_2.inputs.size() > 64) { // TODO: do we need to limit the inputs (because of high runtime?)
    MSG("Currently only simulating up to 64 inputs!");
    return 0;
  }


  DoSimEqModTwo(spec_filename_1, spec_filename_2, new_spec_cir_1, new_spec_cir_2, 0, amount);
  return 0;
}

void changeTypeToLoop(int &real, Node *tmp, Node *tmp2) {
  if (real == 1 || real == 2) { real --;}
  else {
    if (real == 3) {
      if (tmp->inputs[0]->name == tmp2->inputs[0]->name) {real = 2;}
      else {real=3;}
    }
    else {
      if (real == 4) {
        if (tmp->inputs[0]->name == tmp2->inputs[0]->name) {real = 4;}
        else { real=5; }
      } else { real ++; }
    }
  }
}

void findPrimaryOutput(Circuit *spec_cir, map<string,int>&fin_len, map<string, set<string>> &effect_map) {
  NodeVector outputs = spec_cir->outputs;
  effect_map.clear();
  fin_len.clear();
  Node* node;

  for (long i=0; i<outputs.size(); i++) {
    node = outputs[i];
    effect_map[node->name].insert(node->name);
    for (int j=0; j<node->inputs.size(); j++) {
      effect_map[node->inputs[j]->name].insert(node->name);
      fin_len[node->inputs[j]->name] = effect_map[node->inputs[j]->name].size();
    }
  }

  for (long i=spec_cir->all_nodes.size()-1; i>=0; i--) {
    node = spec_cir->all_nodes[i];
    for (int j=0; j<node->inputs.size(); j++) {
      if (!node->inputs[j]->is_input) {
        // fin_len[node->inputs[j]->name] = min(stoi(),fin_len[node->inputs[j]->name]);
        effect_map[node->inputs[j]->name].insert(effect_map[node->name].begin(), effect_map[node->name].end());
        fin_len[node->inputs[j]->name] = effect_map[node->inputs[j]->name].size();
      }
    }
  }
}

void findPrimaryInput(Circuit *spec_cir, map<string,int>&fin_len) {
  NodeVector inputs = spec_cir->inputs;
  map<string, set<string>> effect_map;
  fin_len.clear();
  Node* node;

  for (long i=0; i<inputs.size(); i++) {
    node = inputs[i];
    for (int j=0; j<node->outputs.size(); j++) {
      effect_map[node->outputs[j]->name].insert(node->name);
      fin_len[node->outputs[j]->name] = effect_map[node->outputs[j]->name].size();
    }
  }

  for (long i=0; i<spec_cir->all_nodes.size(); i++) {
    node = spec_cir->all_nodes[i];
    for (int j=0; j<node->outputs.size(); j++) {
      effect_map[node->outputs[j]->name].insert(effect_map[node->name].begin(), effect_map[node->name].end());
      fin_len[node->outputs[j]->name] = effect_map[node->outputs[j]->name].size();
    }
  }
}

int TwoSimEq(string spec_filename, string impl_filename, int method) {
  Circuit new_spec_cir, new_impl_cir;
  new_spec_cir.ReadBlif(spec_filename);
  new_impl_cir.ReadBlif(impl_filename);

  switch (method) {
    case 0: {
      cout << ">>>>>>Random simulation on circuits!" << endl;
      RndSimEq(new_spec_cir, new_impl_cir, 50000);
    }
    break;
    case 1:
    case 2: {
      if (method==1) {
        cout << ">>>>Random simulation on circuits! " << endl;
        RndSimEq(new_spec_cir, new_impl_cir, 1000000);
      }
      if (method==2) WCAEcheck(new_spec_cir, new_impl_cir);

      long int size1 = 0;
      long int size2 = 0;
      long int max_level1 = 0;
      long int max_level2 = 0;
      for (int i=0; i<new_spec_cir.all_nodes.size(); i++) {
        if (new_spec_cir.all_nodes[i]->inputs.size()==2) size1++;
	max_level1 = max(max_level1, new_spec_cir.all_nodes[i]->level);
      }
      new_impl_cir.Simplify();
      // new_impl_cir.Simplify2();
      new_impl_cir.SetIndexes();
      for (int i=0; i<new_impl_cir.all_nodes.size(); i++) {
        if (new_impl_cir.all_nodes[i]->inputs.size()==2) size2++;
	      max_level2 = max(max_level2, new_impl_cir.all_nodes[i]->level);
      }
      cout << "Total Level: " << max_level1 << " Total level: " << max_level2 << endl;
      cout << "original size: " << size1 << " current size: " << size2 << endl;
    }
    break;
    default:
    break;
  }

  return 0;
}

int returnPattern(NodeType nt, int which) {
  int rt = -1;
  switch(nt) {
    case NODE_AND:
    rt = 1;
    break;
    case NODE_NAND:
    if (which == -1) rt = 0;
    else rt = 1;
    break;
    case NODE_OR:
    rt = 0;
    break;
    case NODE_NOR:
    if (which == -1) rt = 1;
    else rt = 0;
    break;
    case NODE_XOR:
    rt = 1;
    break;
    case NODE_XNOR:
    rt = 0;
    break;
    case NODE_AND2_PN:
    if (which != 1) rt = 1;
    else rt = 0;
    //10 1 w=1 0, w=others, 1
    break;
    case NODE_AND2_NP:
    if (which!=0) rt = 1;
    else rt = 0;
    // 01 1, whic
    break;
    case NODE_NAND2_PN:
    if (which!=0) rt = 0;
    else rt = 1;
    // 10 0, which =
    break;
    case NODE_NAND2_NP:
    if (which!=1) rt = 0;
    else rt = 1;
    break;
    default:
    break;
  }
  return rt;
}

int getSize(Circuit *spec_cir) {
  int size1 = 0;
  for (int i=0; i<spec_cir->all_nodes.size(); i++) {
    if (spec_cir->all_nodes[i]->inputs.size()==2) size1++;
  }
  return size1;
}

int getLeastSigBit(set<string> fanout) {
  int least = 1000000;
  set<string>::iterator it;
  for (it=fanout.begin(); it!=fanout.end(); it++) {
    least = min(least, stoi((*it).substr(1, (*it).length())));
  }
  return least;
}

int TestSimEq(string spec_filename, int method) {
  Circuit new_spec_cir;
  string filename = spec_filename;
  size_t dotp = filename.find('.');
  filename = filename.substr(0,dotp);

  new_spec_cir.ReadBlif(spec_filename);

  if (new_spec_cir.inputs.size() == 0) {
    MSG("Cannot simulate circuit with no primary input!");
    return 0;
  }

  // if (new_spec_cir.inputs.size() > 64) { // TODO: do we need to limit the inputs (because of high runtime?)
  //   MSG("Currently only simulating up to 64 inputs!");
  //   return 0;
  // }

  switch (method) {
    case 0: {
      Circuit* c1 = new_spec_cir.GetDuplicate("","","");
      cout << "given circuit has " << c1->all_nodes.size() << " 's nodes'" << endl;
      cout << "please select one node" << endl;
      int i;
      cin >> i;
      std::set<string> effect_nodes;
      Node *n1 = c1->all_nodes[i];
      cout << "You select node " << n1->name << endl << "Effect nodes has " << endl;
      fanouts(n1, effect_nodes);
      int everyline = 5;
      set<string>::iterator it;
      for (it=effect_nodes.begin(); it!=effect_nodes.end(); it++) {
        cout << *it << "\t";
        everyline --;
        if (everyline<0) {
          cout << endl;
          everyline = 5;
        }
      }
      cout << endl;
      cout << "please input a node" << endl;
      string test_node;
      cin >> test_node;
      cout << "the test node is on effect_nodes " << notfouts(effect_nodes, test_node) << endl;

    }break;
    case 3: {
      string nn;
      Circuit *c1 = new_spec_cir.GetDuplicate("","","");
      cout << "Please input the target node name " << endl;
      cin >> nn;
      cout << "Please input the type" << endl;
      cout << "ZERO 0, ONE: 1, BUF0: 2, BUF1: 3, NOT0: 4, NOT1: 5" << endl;
      int t;
      cin >> t;
      GetWrongCircuit(*c1, nn, t);
      cout << "Target Node " << nn << " TYPE TO " << c1->all_nodes_map[nn]->type << endl;
      DoSimEq(new_spec_cir, *c1);
    }break;
    case 1: {
      Circuit *c1 = new_spec_cir.GetDuplicate("","","");
      c1->Simplify();
      c1->Simplify2();
      c1->SetIndexes();
      c1->WriteBlif("see.blif");
      delete c1;

    }break;
    case 2: {
      Circuit* spec_cir = new_spec_cir.GetDuplicate("", "", "");
      cout << "===All outputs for the specified circuit are====" << endl;
      for (int i=0; i<spec_cir->outputs.size(); i++) {
        cout << spec_cir->outputs[i]->name << "\t";
      }
      cout << endl;
      string outputName;
      cin>>outputName;
      cout << "Target Node " << outputName << endl;
      if (spec_cir->all_nodes_map.find(outputName) == spec_cir->all_nodes_map.end()) {
        cout << "---Cannot find the target node in the circuit, Please verify your node input ------" << endl;
        return 0;
      }
      std::set<string> effect_nodes;
      findEffectNodes(spec_cir, outputName, effect_nodes);
      int everyline = 5;
      set<string>::iterator it;
      for (it=effect_nodes.begin(); it!=effect_nodes.end(); it++) {
        cout << *it << "\t";
        everyline --;
        if (everyline<0) {
          cout << endl;
          everyline = 5;
        }
      }
      cout << endl;




    }break;
    case 4: {
      ifstream file(filename+"type2.csv");
      std::vector<string> typeofchange;
      std::string line;
      Circuit *impl_cir;
      std::vector<std::string> types;
      int whichtype = 0;
      while (file.good()) {
        getline(file, line);
        if (line.length() > 4) {
          impl_cir = new_spec_cir.GetDuplicate("","","");
          types.clear();
          types = split(line, ',');
          whichtype = stoi(types[2]);
          if (whichtype > 4) whichtype++;
          else {
            switch (whichtype) {
              case 1:
                whichtype --;
                break;
              case 2:
                whichtype --;
                break;
              case 3:
                whichtype --;
                break;
              case 4:
                break;
            }
          }
          GetWrongCircuit(*impl_cir, types[0], whichtype);
          whichtype = stoi(types[5]);
          if (whichtype > 4) whichtype++;
          else {
            switch (whichtype) {
              case 1:
                whichtype --;
                break;
              case 2:
                whichtype --;
                break;
              case 3:
                whichtype --;
                break;
              case 4:
                break;
            }
          }
          GetWrongCircuit(*impl_cir, types[3], whichtype);

          DoSimEqModTwo("mmtest.blif", "mmtest.blif", new_spec_cir, *impl_cir, 0, -1);
        }
      }
    }break;
    case 5: {
      ifstream file(filename+".csv");
      ofstream outFile;
      outFile.open(filename+"attempt0.csv",  ios::out | ios::app);

      std::vector<string> typeofchange;
      std::string line;
      Circuit *impl_cir;
      std::vector<std::string> types;
      int whichtype;

      double min_err;

      while (file.good()) {
        getline(file, line);
        if (line.length() > 4) {
          impl_cir = new_spec_cir.GetDuplicate("","","");
          types.clear();
          types = split(line, ',');
          whichtype = stoi(types[2]);

          cout << types[0] << " -> " << whichtype << endl;

          if (whichtype > 4) whichtype++;
          else {
            switch (whichtype) {
              case 1:
                whichtype --;
                break;
              case 2:
                whichtype --;
                break;
              case 3:
                whichtype --;
                break;
              case 4:
                break;
            }
          }

          GetWrongCircuit(*impl_cir, types[0], whichtype);

          whichtype = stoi(types[5]);

          if (whichtype > 4) whichtype++;
          else {
            switch (whichtype) {
              case 1:
                whichtype --;
                break;
              case 2:
                whichtype --;
                break;
              case 3:
                whichtype --;
                break;
              case 4:
                break;
            }
          }
          GetWrongCircuit(*impl_cir, types[3], whichtype);

          cout << types[0] << " -> " << whichtype << endl;

          DoTestEq(*impl_cir, min_err);
          outFile << types[0] << "," << types[3] << "," << min_err << endl;
          delete impl_cir;
        }
      }
    } break;
    case 6: {
      ifstream file(filename+".csv");
      ofstream outFile;
      outFile.open(filename+"attempt0.csv",  ios::out|ios::app);

      std::vector<string> typeofchange;
      std::string line;
      Circuit *impl_cir;
      std::vector<std::string> types;
      int whichtype;

      double min_err;

      while (file.good()) {
        getline(file, line);
        if (line.length() > 4) {
          impl_cir = new_spec_cir.GetDuplicate("","","");
          types.clear();
          types = split(line, ',');
          whichtype = stoi(types[1]);

          cout << types[0] << " -> " << whichtype << endl;

          if (whichtype > 4) whichtype++;
          else {
            switch (whichtype) {
              case 1:
                whichtype --;
                break;
              case 2:
                whichtype --;
                break;
              case 3:
                whichtype --;
                break;
              case 4:
                break;
            }
          }

          GetWrongCircuit(*impl_cir, types[0], whichtype);

          GetWrongConnect(*impl_cir, types[3], types[2], 0);

          cout << types[2] << " -> " << types[3] << endl;

          DoTestEq(*impl_cir, min_err);
          outFile << types[2] << "," << types[3] << "," << min_err << endl;

          delete impl_cir;
        }
      }
    } break;
    case 7: {
      ifstream file(filename+".csv");
      ofstream outFile;
      outFile.open(filename+"attempt0.csv",  ios::out | ios::app);

      std::vector<string> typeofchange;
      std::string line;
      Circuit *impl_cir;
      std::vector<std::string> types;
      int whichtype;

      double min_err;

      while (file.good()) {
        getline(file, line);
        if (line.length() > 4) {
          impl_cir = new_spec_cir.GetDuplicate("","","");
          types.clear();
          types = split(line, ',');
          whichtype = stoi(types[3]);

          if (whichtype > 4) whichtype++;
          else {
            switch (whichtype) {
              case 1:
                whichtype --;
                break;
              case 2:
                whichtype --;
                break;
              case 3:
                whichtype --;
                break;
              case 4:
                break;
            }
          }

          GetWrongCircuit(*impl_cir, types[0], whichtype);
          cout << types[0] << " -> " << whichtype << endl;

          DoTestEq(*impl_cir, min_err);

          outFile << types[0] << "," << types[3] << "," << min_err << endl;
          delete impl_cir;
        }
      }
    } break;
    case 8: {
      /*
        try to change each ndoe to output
      */
      // Circuit *spec_cir = new_spec_cir.GetDuplicate("","","");
      // long int node_size = spec_cir->all_nodes.size();
      // NodeVector output_nodes = spec_cir->outputs;
      // string output_name = output_nodes[output_nodes.size()-1]->name;
      // Circuit *impl_cir;
      // Node *cur_node = new Node;
      // string cur_node_name = "";
      // Node *impl_node = new Node;
      // Node *tmp_node1 = new Node;
      // Node *tmp_node2 = new Node;
      //
      // for (long int i = 0; i < node_size; i++) {
      //   cur_node = spec_cir->all_nodes[i];
      //   cur_node_name = cur_node->name;
      //   if (cur_node->is_output or cur_node->is_input or cur_node_name == output_name) continue;
      //   impl_cir = spec_cir->GetDuplicate("","","");
      //   impl_node = impl_cir->all_nodes[i];
      //   tmp_node1 = impl_cir->all_nodes[i]->inputs[0];
      //   tmp_node2 = impl_cir->all_nodes[i]->inputs[1];
      //   impl_cir->outputs[output_nodes.size()-1]->inputs[0] = tmp_node1;
      //   impl_cir->outputs[output_nodes.size()-1]->inputs[1] = tmp_node2;
      //   cout << "ORIGINAL INPUTS: " << spec_cir->all_nodes_map[output_name]->inputs[0]->name << " " << spec_cir->all_nodes_map[output_name]->inputs[1]->name << endl;
      //   cout << "CURRENT INPUTS: " << impl_cir->all_nodes_map[output_name]->inputs[0]->name << " " << impl_cir->all_nodes_map[output_name]->inputs[1]->name << endl;
      //   impl_cir->Simplify();
      //   impl_cir->SetIndexes();
      //   DoSimEq(*spec_cir, *impl_cir);
      //   delete impl_cir;
      // }
    }break;
    case 9: {
      Circuit *spec_cir = new_spec_cir.GetDuplicate("","","");
      string node_name;
      cin >> node_name;
      Node *n1 = spec_cir->all_nodes_map[node_name];
      NodeVector nv1;

      nv1.clear();
      nv1.push_back(n1);

      spec_cir->Simplify(nv1);
      spec_cir->SetIndexes();

      cout << "After Simplified, circuit size: " << spec_cir->all_nodes.size() << endl;

      spec_cir->WriteBlif("spec.blif");

    }break;
    case 10: {
      cout << "This mode is to get specific nodes" << endl;
      Circuit *spec_cir = new_spec_cir.GetDuplicate("","","");
      NodeVector all_nodes;
      GetPossibleNodes(spec_cir, "", all_nodes, 0);
      std::set<string> effect_nodes;
      std::vector<string> target_nodes;
      std::vector<vector<string>> target_detail;
      std::vector<string> tmp_detail;
      int pIs = 0;
      set<string>::iterator it;
      for (int i=0; i<all_nodes.size(); i++) {
        effect_nodes.clear();
        findEffectNodes(spec_cir, all_nodes[i]->name, effect_nodes);
        // cout << all_nodes[i]->name << "," << effect_nodes.size() << endl;
        pIs = 0;
        tmp_detail.clear();
        for (it=effect_nodes.begin(); it!=effect_nodes.end(); it++) {
          if (spec_cir->all_nodes_map[*it]->is_input) {
            tmp_detail.push_back(*it);
            pIs ++;
          }
        }
        if (pIs>=9) {
          target_detail.push_back(tmp_detail);
          target_nodes.push_back(all_nodes[i]->name);
        }
      }
      for (int i=0; i<target_nodes.size(); i++) {
        cout << target_nodes[i] << "," << target_detail[i].size() << endl;
        for (int j=0; j<target_detail[i].size(); j++) {
          cout << target_detail[i][j] << "\t";
        }
        cout << endl;
      }
      cout << endl;
    }break;
    case 13: {

      cout << "Exhaustively Picking up Gates " << endl;
      ofstream outFile;

      outFile.open(filename+".csv", ios::out);

      Circuit *spec_cir = new_spec_cir.GetDuplicate("","","");

      std::vector<Circuit*> impl_cirs;
      NodeVector gates;
      const int thread_num = 6;
      std::thread sim_threads[thread_num];
      vector<bit64> total_diff_count;
      int real;
      Node *tmp;
      Node *tmp2;
      string cur_node_name;
      vector<double> WCAEs;
      vector<double> MAEs;
      std::vector<string> effect_nodes;
      map<string, set<string>> effect_map;


      map<string, int>po_counts;
      map<string, int>pi_counts;
      findPrimaryOutput(spec_cir, po_counts, effect_map);
      findPrimaryInput(spec_cir, pi_counts);

      GetPossibleNodes(spec_cir, "", gates, 0);


      for (int i=0; i<gates.size(); i++) {
        cur_node_name = gates[i]->name;
        spec_cir = new_spec_cir.GetDuplicate("","","");
        GetWrongCircuit(spec_cir, impl_cirs, cur_node_name);
        WCAEs.clear();
        WCAEs.resize(thread_num);

        MAEs.clear();
        MAEs.resize(thread_num);

        total_diff_count.clear();
        total_diff_count.resize(thread_num);

        effect_nodes.clear();
        findEffectNodes(spec_cir, cur_node_name, effect_nodes);


        for (int k=0; k<thread_num; k++) {
          sim_threads[k] = std::thread(DoSimEqTh2, spec_cir, impl_cirs[k], &WCAEs[k], &MAEs[k], &total_diff_count[k]);
        }

        for (int k=0; k<thread_num; k++) {
          sim_threads[k].join();
        }

        int flag = 0;
        for (int k=0; k<thread_num; k++) {
          real = impl_cirs[k]->all_nodes_map[cur_node_name]->type;
          changeTypeToLoop(real, spec_cir->all_nodes_map[cur_node_name], impl_cirs[k]->all_nodes_map[cur_node_name]);
          outFile <<cur_node_name << "," << real << "," <<spec_cir->all_nodes_map[cur_node_name]->type << "," << spec_cir->all_nodes_map[cur_node_name]->level<< "," << total_diff_count[k] << "," << WCAEs[k] << "," << MAEs[k] << "," << getLeastSigBit(effect_map[cur_node_name]) << "," << pi_counts[cur_node_name];
          // for (int m = 0; m < effect_nodes.size(); m++) {
          //   outFile << "," << effect_nodes[m];
          // }
          outFile << endl;
          delete impl_cirs[k];
        }

        delete spec_cir;
      }
    }
    break;
    case 14: {
      ofstream outFile;

      outFile.open(filename+"_2G.csv", ios::out);

      Circuit *spec_cir = new_spec_cir.GetDuplicate("","","");

      std::vector<Circuit*> impl_cirs;
      std::vector<Circuit*> real_cirs;
      NodeVector gates;
      const int thread_num = 6;
      std::thread sim_threads[thread_num];
      vector<bit64> po_diff_count[thread_num];
      vector<bit64> total_diff_count[thread_num];
      int real;
      int real2;
      Node *tmp;
      Node *tmp2;
      string cur_node_name;
      vector<double> WCAEs;
      vector<double> MAEs;
      Circuit *impl_cir;
      NodeVector second_nodes;
      Circuit* real_cir;
      string second_node_name;
      vector<bit64> total_diff;


      GetPossibleNodes(spec_cir, "", gates, 0);


      for (int i=0; i<gates.size(); i++) {
        cur_node_name = gates[i]->name;
        spec_cir = new_spec_cir.GetDuplicate("","","");
        GetWrongCircuit(spec_cir, impl_cirs, cur_node_name);

        for (int j=0; j<thread_num; j++) {
          impl_cir = impl_cirs[j];
          real = impl_cir->all_nodes_map[cur_node_name]->type;
          changeTypeToLoop(real, spec_cir->all_nodes_map[cur_node_name], impl_cir->all_nodes_map[cur_node_name]);

          GetPossibleNodes(impl_cir, cur_node_name, second_nodes, 1);

          for (int m=0; m < gates.size(); m++) {
            if ( gates[m]->name == cur_node_name) continue;
            real_cir = impl_cir->GetDuplicate("","","");
            second_node_name = gates[m]->name;
            GetWrongCircuit(real_cir, real_cirs, gates[m]->name);
            WCAEs.clear();
            WCAEs.resize(thread_num);

            MAEs.clear();
            MAEs.resize(thread_num);

            total_diff.clear();
            total_diff.resize(thread_num);


            for (int k=0; k<thread_num; k++) {
              sim_threads[k] = std::thread(DoSimEqTh2, spec_cir, real_cirs[k], &WCAEs[k], &MAEs[k], &total_diff[k]);
            }

            for (int k=0; k<thread_num; k++) {
              sim_threads[k].join();
            }

            for (int k=0; k<thread_num; k++) {

              real2 = real_cirs[k]->all_nodes_map[second_node_name]->type;
              changeTypeToLoop(real2, spec_cir->all_nodes_map[second_node_name], real_cirs[k]->all_nodes_map[second_node_name]);
              cout <<cur_node_name << "," << real << "," << second_node_name<< "," << real2 << "," << total_diff[k] << "," << WCAEs[k] << "," << MAEs[k] << endl;

              outFile <<cur_node_name << "," << real << "," << second_node_name<< "," << real2 << "," << total_diff[k] << "," << WCAEs[k] << "," << MAEs[k] << endl;
              delete real_cirs[k];
            }
            delete real_cir;

          }
          delete impl_cirs[j];

        }

        delete spec_cir;
      }
    }
    break;
    case 15: {
      ofstream outFile;

      outFile.open(filename+"_1G.csv", ios::out);

      Circuit *spec_cir = new_spec_cir.GetDuplicate("","","");
      Circuit *correct_cir = new_spec_cir.GetDuplicate("","","");

      std::vector<Circuit*> impl_cirs;
      NodeVector gates;
      const int thread_num = 6;
      std::thread sim_threads[thread_num];
      vector<bit64> total_diff_count;
      int real;
      Node *tmp;
      Node *tmp2;
      string cur_node_name;
      vector<double> WCAEs;
      vector<double> MAEs;
      Circuit *tmp_cir = new_spec_cir.GetDuplicate("","","");

      cout << "Please input the node and the type" << endl;
      string n1, n2;
      int t1,t2;

      cin>>n1;
      cin>>t1;
      cin>>n2;
      cin>>t2;

      GetWrongCircuit(*tmp_cir, n1, t1);
      GetWrongCircuit(*tmp_cir, n2, t2);

      tmp_cir->Simplify();
      tmp_cir->Simplify2();
      tmp_cir->SetIndexes();

      GetPossibleNodes(tmp_cir, "", gates, 0);


      for (int i=0; i< gates.size(); i++) {
        cur_node_name = gates[i]->name;
        spec_cir = tmp_cir->GetDuplicate("","","");
        if (spec_cir->all_nodes_map[cur_node_name]->inputs.size() != 2) continue;
        GetWrongCircuit(spec_cir, impl_cirs, cur_node_name);
        WCAEs.clear();
        WCAEs.resize(thread_num);

        MAEs.clear();
        MAEs.resize(thread_num);

        total_diff_count.clear();
        total_diff_count.resize(thread_num);

        for (int k=0; k<thread_num; k++) {
          sim_threads[k] = std::thread(DoSimEqTh2, correct_cir, impl_cirs[k], &WCAEs[k], &MAEs[k], &total_diff_count[k]);
        }

        for (int k=0; k<thread_num; k++) {
          sim_threads[k].join();
        }

        for (int k=0; k<thread_num; k++) {
          real = impl_cirs[k]->all_nodes_map[cur_node_name]->type;
          changeTypeToLoop(real, spec_cir->all_nodes_map[cur_node_name], impl_cirs[k]->all_nodes_map[cur_node_name]);
          outFile <<cur_node_name << "," << real << "," << total_diff_count[k] << "," << WCAEs[k] << "," << MAEs[k] << endl;
          delete impl_cirs[k];
        }

        delete spec_cir;
      }
      delete tmp_cir;
      delete correct_cir;
    }
    break;
    case 16: {
      ofstream out;
      out.open(filename+"_unit.csv", ios::out);

      Circuit *spec_cir = new_spec_cir.GetDuplicate("","","");
      NodeVector gates;
      GetPossibleNodes(spec_cir, "", gates, 0);
      map<string, int>po_counts;
      map<string, set<string>> effect_map;
      findPrimaryOutput(spec_cir, po_counts, effect_map);

      // write preliminary circuit analysis
      for (long int i=0; i<gates.size(); i++) {
        double err_em;
        if (gates[i]->is_output) err_em = 0;
        else err_em = 1 / pow(2,po_counts[gates[i]->name]);
        out << gates[i]->name << "," << gates[i]->level << ","<< po_counts[gates[i]->name]<< "," << err_em << endl;
      }

      out.close();
    }
    break;
    case 17: {
      // cout << "Approximate Circuits for simple pattern" << endl;
      double target_err = 0.01 * 0.01;
      cout << "Please input target error rate " << endl;
      cin >> target_err;

      auto start_time = chrono::system_clock::now();

      Circuit *spec_cir = new_spec_cir.GetDuplicate("","","");
      NodeVector gates;
      GetPossibleNodes(spec_cir, "", gates, 0);
      map<string, int>po_counts;
      map<string, double>pre_wcae;
      map<string, set<string>> effect_map;
      findPrimaryOutput(spec_cir, po_counts, effect_map);
      vector<vector<string>> std_wcae;

      std_wcae.resize(spec_cir->outputs.size()+1);



      for (long int i=0; i<gates.size(); i++) {
        double err_em;
        if (gates[i]->is_output) {
          err_em = 1 / pow(2, spec_cir->outputs.size() - stoi(gates[i]->name.substr(1,gates[i]->name.length())));
          std_wcae[spec_cir->outputs.size()-stoi(gates[i]->name.substr(1,gates[i]->name.length()))-1].push_back(gates[i]->name);
        }
        else {
          if (po_counts[gates[i]->name] == 1) {
            string po_name = *effect_map[gates[i]->name].begin();
            err_em = 1 /  pow(2, spec_cir->outputs.size() - stoi(po_name.substr(1,po_name.length())));
            std_wcae[spec_cir->outputs.size()-stoi(po_name.substr(1,po_name.length()))-1].push_back(gates[i]->name);
          } else {
            err_em = 1 / pow(2,po_counts[gates[i]->name]);
            std_wcae[po_counts[gates[i]->name]].push_back(gates[i]->name);
          }
        }
        pre_wcae[gates[i]->name] = err_em;
      }



      Circuit* answer2 = spec_cir->GetDuplicate("","","");
      Circuit *answer = spec_cir->GetDuplicate("","","");
      vector<double> answers_err;
      vector<string> answer_nodes;

      double cur_err = 0;




      long int origin_size = 0;
      NodeVector base;

      cout << "TAGRET ERROR RATE: " << target_err << endl;

      gates.clear();

      for (int i=0; i<spec_cir->all_nodes.size(); i++) {
        if (spec_cir->all_nodes[i]->inputs.size()==2) {
          origin_size++;
          if (pre_wcae[spec_cir->all_nodes[i]->name] < target_err) {
            gates.push_back(spec_cir->all_nodes[i]);
          }
        }
      }

      vector<string> cur_std_wcae;
      Node* cur_node;
      int flag = 0;
      double ans2_err = 0;


      for (int i = std_wcae.size() - 1; i >=0; i--) {
        if (flag == 1) break;
        cur_std_wcae.clear();
        cur_std_wcae = std_wcae[i];

        for (int j = 0; j < cur_std_wcae.size(); j++) {
          if (answer->all_nodes_map.find(cur_std_wcae[j]) == answer->all_nodes_map.end()) continue;
          cur_node = answer->all_nodes_map[cur_std_wcae[j]];
          if (cur_node->inputs.size() !=2 ) continue;
          cur_err += pre_wcae[cur_node->name];
          ans2_err += pre_wcae[cur_node->name];
          // cout << cur_err << endl;

          if ( cur_err >= target_err) {
            GetWrongCircuit(*answer2, cur_node->name, 0);
            cur_err -= pre_wcae[cur_node->name];
            flag = 1;
            break;
          }
          int s1,s2;
          Circuit *tmp_cir = answer->GetDuplicate("","","");
          GetWrongCircuit(*tmp_cir, cur_node->name, 0);
          tmp_cir->Simplify();
          tmp_cir->SetIndexes();
          s1 = getSize(tmp_cir);
          delete tmp_cir;
          tmp_cir = answer->GetDuplicate("","","");
          GetWrongCircuit(*tmp_cir, cur_node->name, 1);
          tmp_cir->Simplify();
          tmp_cir->SetIndexes();
          s2 = getSize(tmp_cir);
          delete tmp_cir;
          if (s1 <= s2) {
            if (cur_node->type == NODE_XOR) cur_err -= pre_wcae[cur_node->name];
            GetWrongCircuit(*answer, cur_node->name, 0);
            GetWrongCircuit(*answer2, cur_node->name, 0);
          } else {
            GetWrongCircuit(*answer, cur_node->name, 1);
            GetWrongCircuit(*answer2, cur_node->name, 1);
          }
          // cout << cur_node->name << endl;
          answer->Simplify();
          answer->SetIndexes();
          answer2->Simplify();
          answer2->SetIndexes();
          // double WCAE,MAE;
          // bit64 total_diff_count;
          // DoSimEq(*spec_cir, *answer, WCAE, MAE, total_diff_count);
          // cout << cur_node->name << "," << WCAE << "," << cur_err << endl;
        }
      }

      cout << "Estimation is: " << cur_err<< endl;

      auto end_time = chrono::system_clock::now();
      chrono::duration<double> diff_time = end_time - start_time;
      cout << endl << "Generation time is: " << setprecision(3) << diff_time.count() << " s" << endl;


      if (spec_cir->inputs.size() > 30) RndSimEq(*spec_cir, *answer, 50048);
      else WCAEcheck(*spec_cir, *answer);

      // WCAEcheck(*spec_cir, *answer);
      // WCAEcheck(*spec_cir, *answer2);
      answer->WriteBlif(filename+"ap.blif");

    }
    break;
    case 18: {
      auto start_time = chrono::system_clock::now();
      int flag = 0;

      Circuit *spec_cir = new_spec_cir.GetDuplicate("","","");
      NodeVector gates;
      GetPossibleNodes(spec_cir, "", gates, 0);
      map<string, int>po_counts;
      map<string, int>pi_counts;
      map<string, double>pre_wcae;
      map<string, set<string>> effect_map;

      findPrimaryOutput(spec_cir, po_counts, effect_map);
      findPrimaryInput(spec_cir, pi_counts);


      vector<Circuit*> answers;
      Circuit *answer = spec_cir->GetDuplicate("","","");


      double cur_err = 0;
      double target_err = 0.1 * 0.01;

      long int origin_size = 0;
      NodeVector start_nodes;
      vector<NodeVector> start_vector;
      NodeVector start_nodes2;

      cout << "TAGRET ERROR RATE: " << target_err/10 << endl;

      target_err *= 5;


      for (int i=0; i<spec_cir->all_nodes.size(); i++) {
        if (spec_cir->all_nodes[i]->inputs.size()==2) {
          origin_size++;
          Node* cur_node = spec_cir->all_nodes[i];
          if (cur_node->inputs[0]->is_input && cur_node->inputs[1]->is_input) {
            if ( cur_node->inputs[0]->name.substr(0,1) != cur_node->inputs[1]->name.substr(0,1)) {
              // cout<< cur_node->inputs[0]->name.substr(0)<< "," << cur_node->inputs[1]->name.substr(0) << "," << cur_node->name <<endl;
              if (stoi(cur_node->inputs[1]->name.substr(1)) == (spec_cir->inputs.size()/2-1) ) continue;
              start_nodes.push_back(cur_node);
            }
          }
          for (int j=0; j<spec_cir->inputs.size()/2; j++) {
            string h = "a" + to_string(j);
            if (cur_node->inputs[0]->is_input && cur_node->inputs[0]->name.substr(0)==h && !cur_node->inputs[1]->is_input ) {
              start_nodes2.push_back(cur_node);
            }
          }
        }
      }

      double cur_max = 0;

      int count = 0;
      int last_target = -1;
      Node* cur_node;
      NodeVector last_nodes;
      double last_max = 0;
      double this_err;

      int ct;



      // while (cur_err<target_err) {
      start_nodes2.clear();
      for (int i = 0; i < start_nodes.size(); i++) {
        cur_node = start_nodes[i];
        ct = -1;
        if (cur_node->is_output) {
          this_err = 1 / pow(2, spec_cir->outputs.size() - stoi(cur_node->name.substr(1)));
          ct = 0;
        }
        else {
          if (cur_node->type == NODE_AND) {
            if (po_counts[cur_node->name] == 1) {
              this_err = 1 / pow(2, spec_cir->outputs.size() - stoi(cur_node->outputs[0]->name.substr(1)));
            } else this_err = 1 / pow(2, po_counts[cur_node->name]);
            ct = 0;
          }
        }
        if (ct!=-1) {
          cur_err += this_err;
          cout <<  cur_node->name << "," << ct << "," << cur_err << endl;
          GetWrongCircuit(*answer, cur_node->name, ct);
          pre_wcae[cur_node->name] = this_err;
          for (int j=0; j<cur_node->outputs.size(); j++) {
            start_nodes2.push_back(cur_node->outputs[j]);
          }

        }
      }


      for (int i = 0; i < start_nodes2.size(); i++) {
        cur_node = start_nodes2[i];
        ct = -1;
        if (cur_node->is_output) {
          this_err = 1 / pow(2, spec_cir->outputs.size() - stoi(cur_node->name.substr(1)));
          ct = 0;
        }
        else {
          this_err = 1 / pow(2, po_counts[cur_node->name]);
          ct = 1;
        }
        if (ct!=-1 && cur_node->type==NODE_NOR) {
          if (cur_node->inputs[1]->outputs.size() >= 2) continue;
          cur_err += this_err;
          cur_err -= pre_wcae[cur_node->inputs[0]->name];
          cout <<  cur_node->name << "," << ct << "," << cur_err << endl;
          GetWrongCircuit(*answer, cur_node->name, ct);
        }
      }
      //
      //
      // for (int i = 0; i < start_nodes3.size(); i++) {
      //   cur_node = start_nodes3[i];
      //   ct = -1;
      //   if (cur_node->is_output) {
      //     this_err = 1 / pow(2, spec_cir->outputs.size() - stoi(cur_node->name.substr(1)));
      //     ct = 0;
      //   }
      //   else {
      //     this_err = 1 / pow(2, po_counts[cur_node->name]);
      //     ct = 2;
      //   }
      //   if (ct!=-1) {
      //     cur_err += this_err;
      //     cout <<  cur_node->name << "," << ct << "," << cur_err << endl;
      //     GetWrongCircuit(*answer, cur_node->name, ct);
      //   }
      // }
      //
      //
      //
      // for (int i = 0; i < start_nodes4.size(); i++) {
      //   cur_node = start_nodes4[i];
      //   ct = -1;
      //   if (cur_node->is_output) {
      //     this_err = 1 / pow(2, spec_cir->outputs.size() - stoi(cur_node->name.substr(1)));
      //     ct = 0;
      //   }
      //   else {
      //     this_err = 1 / pow(2, po_counts[cur_node->name]);
      //     ct = 2;
      //   }
      //   if (ct!=-1) {
      //     cur_err += this_err;
      //     cout <<  cur_node->name << "," << ct << "," << cur_err << endl;
      //     GetWrongCircuit(*answer, cur_node->name, ct);
      //   }
      // }
      //




      cout << "ESTIMATION ERROR RATE: " << cur_err << endl;

      auto end_time = chrono::system_clock::now();
      chrono::duration<double> diff_time = end_time - start_time;
      cout << endl << "Generation time is: " << setprecision(3) << diff_time.count() << " s" << endl;

      answer->WriteBlif("kaka.blif");
      if (spec_cir->inputs.size() > 30) {
        for (int i=0; i<10; i++) {
          Circuit *iii = answer->GetDuplicate("","","");
          RndSimEq(*spec_cir, *iii, 100000);
          delete iii;
        }
      }
      else  WCAEcheck(*spec_cir, *answer);


    }
    break;
    case 19: {
      auto start_time = chrono::system_clock::now();

      Circuit *spec_cir = new_spec_cir.GetDuplicate("","","");
      NodeVector gates;
      GetPossibleNodes(spec_cir, "", gates, 0);
      map<string, int>po_counts;
      map<string, int>pi_counts;
      map<string, double>pre_wcae;
      map<string, set<string>> effect_map;

      findPrimaryOutput(spec_cir, po_counts, effect_map);
      findPrimaryInput(spec_cir, pi_counts);


      vector<Circuit*> answers;
      Circuit *answer = spec_cir->GetDuplicate("","","");


      double cur_err = 0;
      double target_err = 1 * 0.01;

      long int origin_size = 0;
      NodeVector start_nodes;

      cout << "Please input target error rate " << endl;
      cin >> target_err;


      cout << "TAGRET ERROR RATE: " << target_err << endl;

      target_err *= 1;

      for (int i=0; i<spec_cir->all_nodes.size(); i++) {
        if (spec_cir->all_nodes[i]->inputs.size()==2) {
          origin_size++;
          Node* cur_node = spec_cir->all_nodes[i];
          if (cur_node->inputs[0]->is_input && cur_node->inputs[1]->is_input) {
            start_nodes.push_back(cur_node);
          }
        }
      }

      double cur_max = 0;

      int last_target = -1;
      Node* cur_node;
      NodeVector last_nodes;
      int k = 0;
      double real_min;

      while (cur_err<target_err && start_nodes.size()) {
        last_nodes.clear();
        // cout << cur_err << "," << cur_max << ","<< target_err << "," << start_nodes.size() << endl;
        cur_max = 1;
        real_min = 0;
        for (int i=0; i<start_nodes.size(); i++) {

          real_min = 1;
          int types = -1;
          cur_node = start_nodes[i];

          if (answer->all_nodes_map[cur_node->name]->inputs.size()!=2) continue;
          NodeType cur_type = cur_node->type;
          double this_err;
          if (cur_node->is_output) {
            this_err = 1 / pow(2, spec_cir->outputs.size() - (double)stoi(cur_node->name.substr(1,cur_node->name.size())) );
          } else this_err = 1 / pow(2, po_counts[cur_node->name]);
          // cout << cur_node->name << "," << this_err<<endl;
          cur_max = min(this_err, cur_max);
          real_min = min(this_err, real_min);
          if (this_err > target_err/10) continue;

          if (this_err < target_err) {
            pre_wcae[cur_node->name] = this_err;
            cur_err += this_err;
            if (pre_wcae.find(cur_node->inputs[0]->name) != pre_wcae.end()) cur_err -= pre_wcae[cur_node->inputs[0]->name];
            else cur_err -= pre_wcae[cur_node->inputs[1]->name];
            // if (pre_wcae.find(cur_node->inputs[1]->name) != pre_wcae.end())

            if (cur_node->inputs[0]->is_input && cur_node->inputs[1]->is_input) {
              if (cur_type != NODE_AND) continue;
              else types = 0;
            } else {
              if ( (cur_node->inputs[0]->is_input || cur_node->inputs[1]->is_input) && (cur_type == NODE_AND)) {
                types = 0;
              } else {
                // if (answer->all_nodes_map[cur_node->inputs[0]->name]->inputs.size() == 2 && answer->all_nodes_map[cur_node->inputs[1]->name]->inputs.size() == 2 ) continue;
                if (cur_type == NODE_XNOR || cur_type == NODE_XOR) { types = 0; }
                else {
                  if (returnPattern(cur_type, 0) == answer->all_nodes_map[cur_node->inputs[0]->name]->type - 1
                  || returnPattern(cur_type, 1) == answer->all_nodes_map[cur_node->inputs[1]->name]->type - 1 ) {
                  types = returnPattern(cur_type, -1);
                  } else {
                    types = 1 - returnPattern(cur_type, -1);
                  }
                }
              }
            }

            for (int j=0; j<cur_node->outputs.size(); j++) {
              last_nodes.push_back(cur_node->outputs[j]);
            }
            cout << cur_node->name << "," << types << ","<< cur_err << endl;
            GetWrongCircuit(*answer, cur_node->name, types);
          }
        }
        start_nodes.clear();
        start_nodes = last_nodes;
      }



      cout << "ESTIMATION ERROR RATE: " << cur_err << endl;

      auto end_time = chrono::system_clock::now();
      chrono::duration<double> diff_time = end_time - start_time;
      cout << endl << "Generation time is: " << setprecision(3) << diff_time.count() << " s" << endl;


      if (spec_cir->inputs.size() > 30) RndSimEq(*spec_cir, *answer, 50048);
      else WCAEcheck(*spec_cir, *answer);

      answer->WriteBlif(filename+"ap.blif");
    }
    break;
    case 20: {
      auto start_time = chrono::system_clock::now();

      Circuit *spec_cir = new_spec_cir.GetDuplicate("","","");
      NodeVector gates;
      GetPossibleNodes(spec_cir, "", gates, 0);
      map<string, int>po_counts;
      map<string, int>pi_counts;
      map<string, set<string>> effect_map;

      findPrimaryOutput(spec_cir, po_counts, effect_map);
      findPrimaryInput(spec_cir, pi_counts);


      vector<Circuit*> answers;
      Circuit *answer = spec_cir->GetDuplicate("","","");


      double cur_err = 0;
      double target_err = 0.1 * 0.01;

      long int origin_size = 0;
      NodeVector start_nodes;

      cout << "TAGRET ERROR RATE: " << target_err << endl;

      target_err *= 100;

      for (int i=0; i<spec_cir->all_nodes.size(); i++) {
        if (spec_cir->all_nodes[i]->inputs.size()==2) {
          origin_size++;
          Node* cur_node = spec_cir->all_nodes[i];
          if (cur_node->inputs[0]->is_input && cur_node->inputs[1]->is_input) {
            start_nodes.push_back(cur_node);
          }
        }
      }

      double cur_max = 0;

      int last_target = -1;
      Node* cur_node;
      NodeVector last_nodes;
      int k = 0;
      double real_min;

      while (cur_err<target_err && start_nodes.size()) {
        last_nodes.clear();
        // cout << cur_err << "," << cur_max << ","<< target_err << "," << start_nodes.size() << endl;
        cur_max = 1;
        real_min = 0;
        for (int i=0; i<start_nodes.size(); i++) {

          real_min = 1;
          int types = -1;
          cur_node = start_nodes[i];

          if (answer->all_nodes_map[cur_node->name]->inputs.size()!=2) continue;
          NodeType cur_type = cur_node->type;
          double this_err;
          if (cur_node->is_output) {
            this_err = 1 / pow(2, spec_cir->outputs.size() - (double)stoi(cur_node->name.substr(1,cur_node->name.size())) );
          } else this_err = 1 / pow(2, po_counts[cur_node->name]);
          // cout << cur_node->name << "," << this_err<<endl;
          cur_max = min(this_err, cur_max);
          real_min = min(this_err, real_min);
          if (this_err > target_err/100) continue;

          if (this_err + cur_err < target_err) {
            cur_err += this_err;
            if (cur_node->inputs[0]->is_input && cur_node->inputs[1]->is_input) {
              if (cur_type != NODE_AND) continue;
              else types = 0;
            } else {
              if ( (cur_node->inputs[0]->is_input || cur_node->inputs[1]->is_input) && (cur_type == NODE_AND)) {
                types = 0;
              } else {
                // if (answer->all_nodes_map[cur_node->inputs[0]->name]->inputs.size() == 2 && answer->all_nodes_map[cur_node->inputs[1]->name]->inputs.size() == 2 ) continue;
                if (cur_type == NODE_XNOR || cur_type == NODE_XOR) { types = 0; }
                else {
                  if (returnPattern(cur_type, 0) == answer->all_nodes_map[cur_node->inputs[0]->name]->type - 1
                  || returnPattern(cur_type, 1) == answer->all_nodes_map[cur_node->inputs[1]->name]->type - 1 ) {
                  types = returnPattern(cur_type, -1);
                  } else {
                    types = 1 - returnPattern(cur_type, -1);
                  }
                }
              }
            }

            for (int j=0; j<cur_node->outputs.size(); j++) {
              last_nodes.push_back(cur_node->outputs[j]);
            }
            // cout << cur_node->name << "," << types << ","<< cur_err << endl;
            GetWrongCircuit(*answer, cur_node->name, types);
          }
        }
        cout << cur_err << "," << cur_max << ","<< target_err << "," << start_nodes.size() << endl;

        start_nodes.clear();
        start_nodes = last_nodes;
      }



      cout << "ESTIMATION ERROR RATE: " << cur_err << endl;

      auto end_time = chrono::system_clock::now();
      chrono::duration<double> diff_time = end_time - start_time;
      cout << endl << "Generation time is: " << setprecision(3) << diff_time.count() << " s" << endl;


      // if (spec_cir->inputs.size() > 30) RndSimEq(*spec_cir, *answer, 50048);
      // else
      WCAEcheck(*spec_cir, *answer);


    }
    break;
  }
  return 0;
}

// NEURAL NETWORK PART
template<class T>
int readCSV(string filename,  vector<vector<T>> &inputs) {

  inputs.clear();

  ifstream file(filename);
  string line;
  vector<string> lines;
  vector<T> single_input;

  while (file.good()) {
    getline(file, line);
    if (line.length()>0) {
      lines.clear();
      lines = split(line, ',');
      single_input.clear();
      single_input.resize(lines.size(), 0);
      for (long int i=0; i<lines.size(); i++) {
        single_input[i] = stoi(line) / 255;

      }
      inputs.push_back(single_input);
    }
  }
  return 0;
}

int readCSV(string filename, vector<int> &labels) {

  labels.clear();

  ifstream file(filename);
  string line;
  int label;

  while (file.good()) {
    getline(file, line);
    if (line.length() == 1) {
      label = stoi(line);
      labels.push_back(label);
    }
  }
  return 0;
}

void formatConsole() {
  for (int i=0; i<30; i++) {
    cout << "=";
  }
  cout << endl;
}


int SetUpNN(string dataset_filename, string label_filename) {
  vector<vector<float>> inputs;
  vector<int> labels;
  readCSV(dataset_filename, inputs);
  readCSV(label_filename, labels);
  formatConsole();
  cout << "\t" <<"MODEL DETAIL" << "\t" << endl;
  cout << "TrainSet Size:\t" << inputs.size() << endl;
  cout << "Feature Size: \t" << inputs[0].size() << endl;
  cout << "Number of label: " << labels.size() << endl;
  formatConsole();
  Network network;
  // loat bias=1, int epoch=30, float lr=0.1, bool early_stopping=false
  network.init(1,300,0.01,false);
  // int size = 10;
  // inputs.resize(size, vector<float>(4.0));
  network.fit(inputs, labels);
  return 0;
}
