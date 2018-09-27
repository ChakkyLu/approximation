#include "global.h"
#include <cmath>

#include <bitset>

#include <iomanip>

#include <thread>

using namespace std;

#include "circuit.h"

using namespace nodecircuit;

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
        cout << progress_str << " -> number of minterm differences so far: " << num_diff << endl;
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

// function to be used for multi-threading
void DoSimEqTh2(Circuit *new_spec_cir, Circuit *new_impl_cir) {
  DoSimEq(*new_spec_cir, *new_impl_cir);
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
  cout << endl << "Simulation time is: " << setprecision(3) << diff_time.count() << " s" << endl;

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



int notfouts(NodeVector v1, Node* node){
  int key=1;
  for(int i=0; i<v1.size(); i++){
    if (v1[i] == node) return 0;
    if (v1[i]->level<=node->level) {
      long h1 = v1[i]->level;
      long h2 = node->level;
      if((h2-h1)<50) key=notfouts(v1[i]->outputs, node);
      else return 0;
    }
    if(key==0) return 0;
  }
  return key;
}

void debug_msg(string kk){
  cout << kk << endl;
}

int isGate(Node* n1, Node* n2){
  for (long i=0; i<n1->outputs.size(); i++){
    for (long j=0; j<n2->outputs.size(); j++){
      if (n2->outputs[j]->name == n1->outputs[i]->name) return 1;
    }
  }
  return 0;
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
    spec_cir->ResetLevels();
    spec_cir->LevelizeSortTopological(false);
    spec_cir->SetIndexes();
  }
  switch (mode) {
    case 0:
      for (long i=0; i<spec_cir->all_nodes.size(); i++){
        c1 = spec_cir->all_nodes[i];
        if(c1->inputs.size()==2 && !c1->is_input && !c1->is_output && c1->name!=first_node_name){
          second_nodes.push_back(c1);
        }
      }
      break;
    case 1:
      for (long i=0; i<spec_cir->all_nodes.size(); i++){
        c1 = spec_cir->all_nodes[i];
        if(first_node_name==""){
          if(c1->inputs.size()==2 && !c1->is_input && !c1->is_output && c1->name!=first_node_name){
            second_nodes.push_back(c1);
          }
        }
        if(c1->name!=first_node_name) continue;
        if(c1->name==first_node_name) first_node_name="";
      }
      break;
  }
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
  for (long i=0; i<spec_cir->all_nodes.size(); i++){
    n1 = spec_cir->all_nodes[i];
    if(notfouts(current_node->outputs, n1) and current_node->name != n1->name and current_node->inputs[0]->name != n1->name and n1->name != current_node->inputs[1]->name and !n1->is_input and !n1->is_output){
      possible_connections.push_back(n1);
    }
  }
  connections = possible_connections;
}

void GetWrongCircuit(Circuit &new_spec_cir, string target_node_name, int type){

  Node* target_node = new_spec_cir.all_nodes_map[target_node_name];

  if (target_node->inputs.size() == 2) {

    Node* cur_node = NULL;
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
  // new_spec_cir.Simplify();
  // new_spec_cir.SetIndexes();
}

int RndSimEq(Circuit *new_spec_cir, Circuit *new_impl_cir, int num_patterns) {

  if (new_spec_cir->inputs.size() == 0) {
    MSG("Cannot simulate circuit with no primary input!");
    return 0;
  }

  if (new_spec_cir->inputs.size() != new_impl_cir->inputs.size() || new_spec_cir->outputs.size() != new_impl_cir->outputs.size()) {
    MSG("Number of Pi/Po of two circuits are different!");
    return 0;
  }

  // new_spec_cir->Simplify();
  // new_spec_cir->ResetLevels();
  // new_spec_cir->LevelizeSortTopological(false);
  // new_spec_cir->SetIndexes();
  // //
  // // new_impl_cir->CreateCodedBlifValues();
  // // new_impl_cir->ConvertEquivalentTypeBlif();
  // // new_impl_cir->LevelizeSortTopological(false);
  // // new_impl_cir->RemoveBufNot();
  // new_impl_cir->Simplify();
  // new_impl_cir->ResetLevels();
  // new_impl_cir->LevelizeSortTopological(false);
  // new_impl_cir->SetIndexes();

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
      //cout << std::bitset<64>(new_rnd) << endl;
      // TODO: check it is different from previous simulation pattern(s)
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
  return num_diff;
}

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
  NodeVector change_nodes;
  double needtoadd;
  double count = 0;


  spec_cir = new_spec_cir.GetDuplicate("","","");
  target_node = spec_cir->all_nodes_map[target_node_name];
  change_nodes.clear();
  change_nodes.push_back(target_node);
  spec_cir->Simplify(change_nodes);
  spec_cir->ResetLevels();
  spec_cir->LevelizeSortTopological(false);
  spec_cir->SetIndexes();

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
     // po_diff_count[i][7] *= needtoadd;
     // total_diff_count[i][7] *= needtoadd;
    for (int j=0; j < 8; j++) {
      po_diff_count[i][j] *= needtoadd;
      total_diff_count[i][j] *= needtoadd;
    }
    delete impl_cirs[i];
  }
  delete spec_cir;
  return 0;
}

int DoSimEqMod(Circuit &spec_cir, NodeVector target_nodes, string spec_filename) {

  //ofstream outFile;
 // ofstream outText;
 // string filename = spec_filename;
 // size_t dotp = filename.find('.');
 // filename = filename.substr(0,dotp);
 // outFile.open(filename+".csv", ios::out | ios::app);
 // outText.open(filename+"v1.txt", ios::out | ios::app);


  const int num_changes = 15;
  Circuit* impl_cirs[num_changes];
  Circuit* new_spec_cir;
  std::thread sim_threads[num_changes];
  vector<bit64> po_diff_count[num_changes];
  vector<bit64> total_diff_count[num_changes];
  string target_node_name;
  int size1 = spec_cir.inputs.size();
  Node* target_node = new Node;
  int cur_num = 0;
  NodeVector change_nodes;
  double needtoadd;
  double count = 0;

  for(int m=0; m<target_nodes.size(); m++){

    target_node_name = target_nodes[m]->name;
    new_spec_cir = spec_cir.GetDuplicate("","","");
    target_node = new_spec_cir->all_nodes_map[target_node_name];
    change_nodes.clear();
    change_nodes.push_back(target_node);

    new_spec_cir->Simplify(change_nodes);
    new_spec_cir->ResetLevels();
    new_spec_cir->LevelizeSortTopological(false);
    new_spec_cir->SetIndexes();

    needtoadd = pow(2, (size1 - new_spec_cir->inputs.size()));

    cur_num=0;

    for (int i=0; i<num_changes+1; i++){
      if (i-1 == target_node->type)  continue;
      impl_cirs[cur_num] = new_spec_cir->GetDuplicate("","","");
      GetWrongCircuit(*impl_cirs[cur_num], target_node_name, i);
      sim_threads[cur_num] = std::thread(DoSimEqTh4, new_spec_cir, impl_cirs[cur_num], &po_diff_count[cur_num], &total_diff_count[cur_num]);
      cur_num++;
    }

    for (int t = 0; t < num_changes; t++)
      sim_threads[t].join();
    for (int t = 0; t < num_changes; t++) {
     // cout << endl << "Summary of final result:" << endl;
      //cout << "------target node----"<<target_node->name<<"---------- change " << impl_cirs[t]->all_nodes_map[target_node_name]->type << " ----------" << endl;
      //outText << target_node->inputs[0]->name << "," << target_node->inputs[1]->name  << "," << impl_cirs[t]->all_nodes_map[target_node_name]->type << "," << total_diff_count[t][7]*needtoadd << endl;
      cout << target_node->name << "," << target_node->type << "," << impl_cirs[t]->all_nodes_map[target_node_name]->type << "," << total_diff_count[t][7]*needtoadd << endl;
    }
  }

 // outText.close();
 // outFile.close();
  delete target_node;
  return 0;
}

int DoSimEqModRes(Circuit &new_spec_cir, NodeVector targets_nodes, string spec_filename, int mode) {
  ofstream outFile;
  ofstream outText;

  string filename = spec_filename;
  size_t dotp = filename.find('.');
  filename = filename.substr(0,dotp);

  if(mode==1 || mode==4){
    outFile.open(filename+"withRes.csv",  ios::out | ios::app);
    outText.open(filename+"withRes.txt", ios::out | ios::app);
  }
  if(mode==2){
    outFile.open(filename+"random.csv",  ios::out | ios::app);
    outText.open(filename+"random.txt",  ios::out | ios::app);
  }
  if(mode==3){
    outFile.open(filename+"random_analysis.csv",  ios::out | ios::app);
    outText.open(filename+"random_analysis.csv",  ios::out | ios::app);
  }


  const int num_changes = 15;
  Circuit* spec_cirs[num_changes];
  Circuit* impl_cirs[num_changes];
  std::thread sim_threads[num_changes];
  vector<bit64> po_diff_count[num_changes];
  vector<bit64> total_diff_count[num_changes];

  int numpatterns;
  int pattern_nums[5] = {1000,5000,10000,20000,50000};
  int numberofwrong[5] = {0, 0, 0, 0, 0};
  int threadnum = 0;
  std::thread rand2ex_thread[15];
  string nodename[15]={"","","","","","","","","","","","","","","",};
  long nodelevel[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  int nodechange[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  double inputsizechange[15]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  vector<bit64> total_node_count[15];
  Node* cur_node = NULL;
  int cur_change_num;
  Circuit* spec_cir;
  vector<int> real_sim_num;
  string target_node_name;
  Node* target_node = NULL;



  for(int m=0; m<targets_nodes.size(); m++){
    target_node_name = targets_nodes[m]->name;
    Node* target_node = new_spec_cir.all_nodes_map[target_node_name];
    NodeVector target_nodes;
    target_nodes.push_back(target_node);
    int size1 = new_spec_cir.inputs.size();

    if(mode==1 || mode==4) {
      new_spec_cir.Simplify(target_nodes);
      new_spec_cir.ResetLevels();
      new_spec_cir.LevelizeSortTopological(false);
      new_spec_cir.SetIndexes();
    }


    double needtoadd = 1;
    if(mode==1){
      outText<<"change node is "<<target_node->name<<" node level is "<<target_node->level<<" PIs "<< new_spec_cir.inputs.size()<<endl;
      needtoadd = pow(2, size1 - new_spec_cir.inputs.size());
    }
    // {
    //   cur_change_num = 0;
    //
    //   cout<<"target node is "<< target_node->name<<endl;
    //   spec_cirs[cur_change_num] = spec_cir->GetDuplicate("", "", "");
    //   impl_cirs[cur_change_num] = spec_cir->GetDuplicate("", "", "");
    //   cur_node = impl_cirs[cur_change_num]->all_nodes_map[target_node->name];
    //   cur_node->type = NODE_ZERO;
    //   if (cur_node->inputs[0]->outputs.size() == 1)
    //     cur_node->inputs[0]->outputs.clear();
    //   else {
    //     NodeVector::iterator it = cur_node->inputs[0]->outputs.begin();
    //     while (it != cur_node->inputs[0]->outputs.end() && *it != cur_node) {
    //       ++it;
    //     }
    //     if (it != cur_node->inputs[0]->outputs.end())
    //       cur_node->inputs[0]->outputs.erase(it);
    //   }
    //   if (cur_node->inputs[1]->outputs.size() == 1)
    //     cur_node->inputs[1]->outputs.clear();
    //   else {
    //     NodeVector::iterator it = cur_node->inputs[1]->outputs.begin();
    //     while (it != cur_node->inputs[1]->outputs.end() && *it != cur_node) {
    //       ++it;
    //     }
    //     if (it != cur_node->inputs[1]->outputs.end())
    //       cur_node->inputs[1]->outputs.erase(it);
    //   }
    //   cur_node->inputs.clear();
    //   impl_cirs[cur_change_num]->Simplify();
    //   impl_cirs[cur_change_num]->SetIndexes();
    //   cur_change_num++;
    //
    //   spec_cirs[cur_change_num] = spec_cir->GetDuplicate("", "", "");
    //   impl_cirs[cur_change_num] = spec_cir->GetDuplicate("", "", "");
    //   cur_node = impl_cirs[cur_change_num]->all_nodes_map[target_node->name];
    //   cur_node->type = NODE_ONE;
    //   if (cur_node->inputs[0]->outputs.size() == 1)
    //     cur_node->inputs[0]->outputs.clear();
    //   else {
    //     NodeVector::iterator it = cur_node->inputs[0]->outputs.begin();
    //     while (it != cur_node->inputs[0]->outputs.end() && *it != cur_node) {
    //       ++it;
    //     }
    //     if (it != cur_node->inputs[0]->outputs.end())
    //       cur_node->inputs[0]->outputs.erase(it);
    //   }
    //   if (cur_node->inputs[1]->outputs.size() == 1)
    //     cur_node->inputs[1]->outputs.clear();
    //   else {
    //     NodeVector::iterator it = cur_node->inputs[1]->outputs.begin();
    //     while (it != cur_node->inputs[1]->outputs.end() && *it != cur_node) {
    //       ++it;
    //     }
    //     if (it != cur_node->inputs[1]->outputs.end())
    //       cur_node->inputs[1]->outputs.erase(it);
    //   }
    //   cur_node->inputs.clear();
    //   impl_cirs[cur_change_num]->Simplify();
    //   impl_cirs[cur_change_num]->SetIndexes();
    //
    //   cur_change_num++;
    //
    //   spec_cirs[cur_change_num] = spec_cir->GetDuplicate("", "", "");
    //   impl_cirs[cur_change_num] = spec_cir->GetDuplicate("", "", "");
    //   cur_node = impl_cirs[cur_change_num]->all_nodes_map[target_node->name];
    //   cur_node->type = NODE_BUF;
    //   if (cur_node->inputs[1]->outputs.size() == 1)
    //     cur_node->inputs[1]->outputs.clear();
    //   else {
    //     NodeVector::iterator it = cur_node->inputs[1]->outputs.begin();
    //     while (it != cur_node->inputs[1]->outputs.end() && *it != cur_node) {
    //       ++it;
    //     }
    //     if (it != cur_node->inputs[1]->outputs.end())
    //       cur_node->inputs[1]->outputs.erase(it);
    //   }
    //   cur_node->inputs.resize(1);
    //   impl_cirs[cur_change_num]->Simplify();
    //   impl_cirs[cur_change_num]->SetIndexes();
    //
    //   cur_change_num++;
    //
    //   spec_cirs[cur_change_num] = spec_cir->GetDuplicate("", "", "");
    //   impl_cirs[cur_change_num] = spec_cir->GetDuplicate("", "", "");
    //   cur_node = impl_cirs[cur_change_num]->all_nodes_map[target_node->name];
    //   cur_node->type = NODE_BUF;
    //   if (cur_node->inputs[0]->outputs.size() == 1)
    //     cur_node->inputs[0]->outputs.clear();
    //   else {
    //     NodeVector::iterator it = cur_node->inputs[0]->outputs.begin();
    //     while (it != cur_node->inputs[0]->outputs.end() && *it != cur_node) {
    //       ++it;
    //     }
    //     if (it != cur_node->inputs[0]->outputs.end())
    //       cur_node->inputs[0]->outputs.erase(it);
    //   }
    //   cur_node->inputs[0] = cur_node->inputs[1];
    //   cur_node->inputs.resize(1);
    //   impl_cirs[cur_change_num]->Simplify();
    //   impl_cirs[cur_change_num]->SetIndexes();
    //
    //   cur_change_num++;
    //
    //   spec_cirs[cur_change_num] = spec_cir->GetDuplicate("", "", "");
    //   impl_cirs[cur_change_num] = spec_cir->GetDuplicate("", "", "");
    //   cur_node = impl_cirs[cur_change_num]->all_nodes_map[target_node->name];
    //   cur_node->type = NODE_NOT;
    //   if (cur_node->inputs[1]->outputs.size() == 1)
    //     cur_node->inputs[1]->outputs.clear();
    //   else {
    //     NodeVector::iterator it = cur_node->inputs[1]->outputs.begin();
    //     while (it != cur_node->inputs[1]->outputs.end() && *it != cur_node) {
    //       ++it;
    //     }
    //     if (it != cur_node->inputs[1]->outputs.end())
    //       cur_node->inputs[1]->outputs.erase(it);
    //   }
    //   cur_node->inputs.resize(1);
    //   impl_cirs[cur_change_num]->Simplify();
    //   impl_cirs[cur_change_num]->SetIndexes();
    //
    //   cur_change_num++;
    //
    //   spec_cirs[cur_change_num] = spec_cir->GetDuplicate("", "", "");
    //   impl_cirs[cur_change_num] = spec_cir->GetDuplicate("", "", "");
    //   cur_node = impl_cirs[cur_change_num]->all_nodes_map[target_node->name];
    //   cur_node->type = NODE_NOT;
    //   if (cur_node->inputs[0]->outputs.size() == 1)
    //     cur_node->inputs[0]->outputs.clear();
    //   else {
    //     NodeVector::iterator it = cur_node->inputs[0]->outputs.begin();
    //     while (it != cur_node->inputs[0]->outputs.end() && *it != cur_node) {
    //       ++it;
    //     }
    //     if (it != cur_node->inputs[0]->outputs.end())
    //       cur_node->inputs[0]->outputs.erase(it);
    //   }
    //   cur_node->inputs[0] = cur_node->inputs[1];
    //   cur_node->inputs.resize(1);
    //   impl_cirs[cur_change_num]->Simplify();
    //   impl_cirs[cur_change_num]->SetIndexes();
    //
    //   cur_change_num++;
    //
    //   for (int nn = NODE_AND; nn < NODE_BLIF; nn++) {
    //     if (nn == target_node->type)
    //       continue;
    //     spec_cirs[cur_change_num] = spec_cir->GetDuplicate("", "", "");
    //     impl_cirs[cur_change_num] = spec_cir->GetDuplicate("", "", "");
    //     cur_node = impl_cirs[cur_change_num]->all_nodes_map[target_node->name];
    //     cur_node->type = (NodeType)nn;
    //     cur_change_num++;
    //   }
    // }

    for (int i=0; i<num_changes; i++){
      spec_cirs[i] = new_spec_cir.GetDuplicate("","","");
      impl_cirs[i] = new_spec_cir.GetDuplicate("","","");
      GetWrongCircuit(*impl_cirs[i], target_node_name, i);
      // sim_threads[i] = std::thread(DoSimEqTh4, spec_cirs[i], impl_cirs[i], &po_diff_count[i], &total_diff_count[i]);
    }



    if (mode==1 || mode==4){ // FIRST RANDOM THEN Exhaustive
      for(int t=0; t<num_changes; t++){
        numpatterns = 10000;
        int h = RndSimEq(spec_cirs[t], impl_cirs[t], numpatterns);
        outText<<" change "<< t << " wrong cases "<< h << "in "<<numpatterns<<" 's input pattern cases"<<endl;
        if(h<=100 && threadnum<15) {
          total_node_count[threadnum] = total_diff_count[t];
          rand2ex_thread[threadnum] = std::thread(DoSimEqTh4, spec_cirs[t], impl_cirs[t], &po_diff_count[t], &total_node_count[threadnum]);
          nodename[threadnum] = target_node->name;
          nodelevel[threadnum] = target_node->level;
          nodechange[threadnum] = t;
          inputsizechange[threadnum] = needtoadd;
          threadnum ++;
        }
        if(threadnum==15 || m==targets_nodes.size()-1){
          cout<<"start thread"<<endl;
          for (int t=0; t<threadnum; t++){
            cout<<"node "<<nodename[t]<<" level "<<nodelevel[t]<<" change "<<nodechange[t]<<" has be sent to accurate simulation"<<endl;
            outText<<"node "<<nodename[t]<<" level "<<nodelevel[t]<<" change "<<nodechange[t]<<" has be sent to accurate simulation"<<endl;
          }
          for (int t=0; t<threadnum; t++){
            rand2ex_thread[t].join();
          }
          if(mode==1){
            for (int t = 0; t < threadnum; t++) {
              cout << endl << "Summary of final result:" << endl;
              cout << "-------"<<nodename[t]<<"-----"<<nodelevel[t] << "--------"<<nodechange[t]<<"----------" << endl;
              for (int i = 0; i < 8; i++){
                cout << total_node_count[t][i] << "\t";
              }
              cout<<"total differnces is "<< total_node_count[t][7]<<endl;
              outFile<< nodename[t] <<","<<nodelevel[t]<<","<<nodechange[t]<<","<<total_node_count[t][7]*inputsizechange[t]<<endl;
            }
          }
          threadnum = 0;
        }
      }
    }

    if(mode==2 ){
      for(int t=0; t<num_changes; t++){ //ALL RANDOM
        numpatterns = 50000;
        int h = RndSimEq(spec_cirs[t], impl_cirs[t], numpatterns);
        cout<<"node is: "<<target_node->name <<" node level is "<<target_node->level<<" change "<< t<< " wrong cases "<< h << "in "<<numpatterns<<" 's input pattern cases"<<endl;
        outText<<"node is: "<<target_node->name <<" node level is "<<target_node->level<<" change "<< t<< " wrong cases "<< h << "in "<<numpatterns<<" 's input pattern cases"<<endl;
        outFile<<target_node->name <<","<<target_node->level<<","<<t<<","<<h<<","<<numpatterns<<endl;
      }
    }

    if(mode==3){ //ALL RANDOM UNDER DIFFERENT PATTERNS
      for(int t=0; t<num_changes; t++){
        int h;
        cout<<"node is: "<<target_node->name <<" node level is "<<target_node->level<<" change "<< t;
        outText<<"node is: "<<target_node->name <<" node level is "<<target_node->level<<" change "<< t;
        outFile<<target_node->name <<","<<target_node->level<<","<<t<<",";
        for(int j=0; j<5; j++){
          h = RndSimEq(spec_cirs[t], impl_cirs[t], pattern_nums[j]);
          if(h==0) numberofwrong[j] += 1;
          cout<< " wrong cases "<< h << "in "<<pattern_nums[j]<<" 's input pattern cases"<<"\t";
          outText<< " wrong cases "<< h << "in "<<pattern_nums[j]<<" 's input pattern cases"<<"\t";
          outFile<<h<<",";
        }
        cout<<endl;
        outText<<endl;
        outFile<<endl;
      }
    }
  }
  if(mode==3){
    outFile<<""<<","<<""<<","<<""<<",";
    for(int j=0; j<5; j++){
      outFile<<numberofwrong[j]<<",";
    }
    outFile<<endl;
  }

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

int DoSimEqMulGate(Circuit &new_spec_cir, string spec_filename, NodeVector first_nodes, int mode){
  ofstream outFile;
  ofstream outFile2;
  string filename = spec_filename;
  size_t dotp = filename.find('.');
  filename = filename.substr(0,dotp);
  if(first_nodes.size()==1) {
    outFile.open("test_2G.csv",  ios::out | ios::app);
    outFile2.open("test_2G_0.csv", ios::out | ios::app);
  }

  else {
    outFile.open(filename+"_2G_type.csv",  ios::out | ios::app);
    outFile2.open(filename+"_2G_type_0.csv",  ios::out | ios::app);
  }

  const int num_changes = 15;

  Circuit* impl_cirs[num_changes];
  Circuit* origin_cir = new_spec_cir.GetDuplicate("", "", "");
  Circuit* impl_cir;
  Circuit* spec_cir;
  std::thread sim_threads[num_changes];
  vector<bit64> po_diff_count[num_changes];
  vector<bit64> total_diff_count[num_changes];

  Node* first_node = NULL;
  Node* second_node = NULL;
  string first_node_name;
  string second_node_name;

  double needtoadd = 1;
  double max_error = 0;
  int cur_num = 0;
  NodeVector second_nodes;
  NodeVector target_nodes;

  if (mode==0){
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

          for (int ii = 0; ii < num_changes+1; ii++){
            if (ii-1 == first_node->type) continue;
            impl_cir = spec_cir->GetDuplicate("","","");
            GetWrongCircuit(*impl_cir, first_node_name, ii);
            cur_num = 0;
            for (int jj = 0; jj < num_changes+1; jj++){
              if (jj-1 == second_node->type) continue;
              impl_cirs[cur_num] = impl_cir->GetDuplicate("","","");
              GetWrongCircuit(*impl_cirs[cur_num], second_node_name, jj);
              sim_threads[cur_num] = std::thread(DoSimEqTh4, spec_cir, impl_cirs[cur_num], &po_diff_count[cur_num], &total_diff_count[cur_num]);
              cur_num++;
            }
            for (int k = 0; k < 6; k++){
              sim_threads[k].join();
            }
            for (int t = 0; t < 6; t++) {
              cout << first_node_name << "," << first_node->type << "," << impl_cirs[t]->all_nodes_map[first_node_name]->type << "," << second_node_name << "," << second_node->type << "," << impl_cirs[t]->all_nodes_map[second_node_name]->type << "," << total_diff_count[t][7]*needtoadd << endl;
              if(total_diff_count[t][7] < max_error){
                outFile << first_node_name << "," << first_node->type << "," << impl_cirs[t]->all_nodes_map[first_node_name]->type << "," << second_node_name << "," << second_node->type << "," << impl_cirs[t]->all_nodes_map[second_node_name]->type << "," << total_diff_count[t][7]*needtoadd << endl;
              }
              cout << endl;
            }
            for (int k = 6; k < num_changes; k++){
              sim_threads[k].join();
            }
            for (int t = 6; t < num_changes; t++) {
              cout << first_node_name << "," << first_node->type << "," << impl_cirs[t]->all_nodes_map[first_node_name]->type << "," << second_node_name << "," << second_node->type << "," << impl_cirs[t]->all_nodes_map[second_node_name]->type << "," << total_diff_count[t][7]*needtoadd << endl;
              if(total_diff_count[t][7] < max_error){
                outFile << first_node_name << "," << first_node->type << "," << impl_cirs[t]->all_nodes_map[first_node_name]->type << "," << second_node_name << "," << second_node->type << "," << impl_cirs[t]->all_nodes_map[second_node_name]->type << "," << total_diff_count[t][7]*needtoadd << endl;
              }
              cout << endl;
            }
          }
        }
      }
    }
  }

  if (mode==1){
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
          for (int ii = 0; ii < num_changes+1; ii++){
            if (ii-1 == first_node->type) continue;
            impl_cir = spec_cir->GetDuplicate("","","");
            GetWrongCircuit(*impl_cir, first_node_name, ii);
            cur_num = 0;
            for (int jj = 0; jj < num_changes+1; jj++){
              if (jj-1 == second_node->type) continue;
              impl_cirs[cur_num] = impl_cir->GetDuplicate("","","");
              GetWrongCircuit(*impl_cirs[cur_num], second_node_name, jj);
              sim_threads[cur_num] = std::thread(DoSimEqTh4, spec_cir, impl_cirs[cur_num], &po_diff_count[cur_num], &total_diff_count[cur_num]);
              cur_num++;
            }
            for (int k = 0; k < num_changes; k++){
              sim_threads[k].join();
            }
            for (int k = 0; k < num_changes; k++){
              cout << first_node_name << "," << first_node->type << "," << impl_cirs[k]->all_nodes_map[first_node_name]->type << "," << second_node_name << "," << second_node->type << "," << impl_cirs[k]->all_nodes_map[second_node_name]->type << "," << total_diff_count[k][7] * needtoadd << endl;
              if(total_diff_count[k][7]>0){
                outFile << first_node_name << "," << first_node->type << "," << impl_cirs[k]->all_nodes_map[first_node_name]->type << "," << second_node_name << "," << second_node->type << "," << impl_cirs[k]->all_nodes_map[second_node_name]->type << "," << total_diff_count[k][7] * needtoadd;
                if(checkPattern(total_diff_count[k])==0){
                  for(int ii=0; ii<8; ii++){
                    outFile << "," << total_diff_count[k][ii]*needtoadd;
                  }
                }
                outFile << endl;
              }else{
                outFile << first_node_name << "," << first_node->type << "," << impl_cirs[k]->all_nodes_map[first_node_name]->type << "," << second_node_name << "," << second_node->type << "," << impl_cirs[k]->all_nodes_map[second_node_name]->type << "," << total_diff_count[k][7] * needtoadd;
              }
            }
          }
        }
      }
    }

  }
  return 0;
}

int DoknownSim(Circuit *new_spec_cir, vector<string> target_nodes, vector<int> types, int mode){

  int t1, t2, type ;
  string target_node_name ,nn1, nn2;
  if(mode<3) {
    target_node_name = target_nodes[0];
    Node* target_node = new_spec_cir->all_nodes_map[target_node_name];
    NodeVector target_nodes;
    type = types[0];
    target_nodes.push_back(target_node);
    new_spec_cir->Simplify(target_nodes);
    new_spec_cir->ResetLevels();
    new_spec_cir->LevelizeSortTopological(false);
    new_spec_cir->SetIndexes();
  }
  if(mode==3){
    nn1 = target_nodes[0];
    nn2 = target_nodes[1];
    Node* n1 = new_spec_cir->all_nodes_map[nn1];
    Node* n2 = new_spec_cir->all_nodes_map[nn2];
    t1 = types[0];
    t2 = types[1];
    NodeVector target_nodes;
    target_nodes.push_back(n1);
    target_nodes.push_back(n2);
    new_spec_cir->Simplify(target_nodes);
    new_spec_cir->ResetLevels();
    new_spec_cir->LevelizeSortTopological(false);
    new_spec_cir->SetIndexes();

  }

  // if(mode==0) new_spec_cir->Simplify();
  // else new_spec_cir->Simplify(target_nodes);
  // new_spec_cir->ResetLevels();
  // new_spec_cir->LevelizeSortTopological(false);
  // new_spec_cir->SetIndexes();

    Circuit* spec_cirs;
    Circuit* impl_cirs;

    Node* cur_node = NULL;

    spec_cirs = new_spec_cir->GetDuplicate("", "", "");
    impl_cirs = new_spec_cir->GetDuplicate("", "", "");
    int rnd;

    // GetWrongCircuit(*impl_cirs,target_node_name,type);

    switch(mode){
      case 0:
        GetWrongCircuit(*impl_cirs,target_node_name,type);
        DoSimEq(*spec_cirs,*impl_cirs);
        break;
      case 1:
        GetWrongCircuit(*impl_cirs,target_node_name,type);
        DoSimEq(*spec_cirs,*impl_cirs);
        break;
      case 2:
        GetWrongCircuit(*impl_cirs,target_node_name,type);
        return RndSimEq(spec_cirs,impl_cirs, 10000);
      case 3:
        // new_spec_cirs = &impl_cirs;
        GetWrongCircuit(*impl_cirs,nn1,t1);
        GetWrongCircuit(*impl_cirs,nn2,t2);
        // impl_cirs->Simplify();
        // impl_cirs->SetIndexes();
        DoSimEq(*spec_cirs,*impl_cirs);
        cout<<"need to be changed to "<< t1<<"    "<<t2<<endl;
        cout<<"now the type of node 1 is "<< impl_cirs->all_nodes_map[nn1]->type<<endl;
        cout<<"now the type of node 2 is "<< impl_cirs->all_nodes_map[nn2]->type<<endl;
        break;
    }
  return 0;
}

int DoSimEqMod(string spec_filename, int method) {
  Circuit new_spec_cir;
  string filename = spec_filename;
  size_t dotp = filename.find('.');
  filename = filename.substr(0,dotp);

  new_spec_cir.ReadBlif(spec_filename);

  if (new_spec_cir.inputs.size() == 0) {
    MSG("Cannot simulate circuit with no primary input!");
    return 0;
  }

  if (new_spec_cir.inputs.size() > 64) { // TODO: do we need to limit the inputs (because of high runtime?)
    MSG("Currently only simulating up to 64 inputs!");
    return 0;
  }

  if (method==1){
    const int num_changes = 15;
    string node_name;
    cin>>node_name;
    Node* n1 = new_spec_cir.all_nodes_map[node_name];
    NodeVector nv1;
    nv1.push_back(n1);
    vector<bit64> po_diff_count[num_changes];
    vector<bit64> total_diff_count[num_changes];
    DoSimEqMod(new_spec_cir, nv1, filename);
    delete n1;
  }
  // srand((unsigned)time( NULL));
  // // string node_name;
  // // int type = rand()%15;
  // int mode;
  //
  // Circuit* c1;
  // Circuit* c2;
  // if(method==0){
  //   string nn;
  //   c1 = new_spec_cir.GetDuplicate("","","");
  //   int seed = rand()%(c1->all_nodes.size() - 50 )+ 50;
  //   // Node* n1 = c1->all_nodes_map["n364"];
  //   Node* n1 = c1->all_nodes[seed];
  //   // type = 6;
  //   while (n1->inputs.size()!=2 || n1->is_input || n1->is_output){
  //     seed = rand()%(c1->all_nodes.size() - 50 )+ 50;
  //     n1 = c1->all_nodes[seed];
  //   }
  //   int t1 = rand()%15;
  //   while(n1->type==t1) t1 = rand()%15;
  //   nn = n1->name;
  //   vector<string> node_name;
  //   node_name.push_back(nn);
  //   std::vector<int> type;
  //   type.push_back(t1);
  //   cout<<"no simplify, exhaustive simulation"<<endl;
  //   DoknownSim(c1, node_name, type, 0);
  //   cout<<"after simplify, exhaustive simulation"<<endl;
  //   DoknownSim(c1, node_name, type, 1);
  //   cout<<"after simplify, random simulation"<<endl;
  //   DoknownSim(c1, node_name, type, 2);
  // }
  // if(method==1){
  //   c1 = new_spec_cir.GetDuplicate("","","");
  //   string node_name_1;
  //   string node_name_2;
  //   int n1_type;
  //   int n2_typ2;
  //   cin>>node_name_1;
  //   cin>>node_name_2;
  //   cin>>n1_type;
  //   cin>>n2_typ2;
  //   vector<string> node_name;
  //   node_name.push_back(node_name_1);
  //   node_name.push_back(node_name_2);
  //   std::vector<int> type;
  //   type.push_back(n1_type);
  //   type.push_back(n2_typ2);
  //   cout<<"first_node is "<<node_name_1<< " and its type is "<< c1->all_nodes_map[node_name_1]->type<< " "<<endl;
  //   cout<<"second_node is "<<node_name_2<< " and its type is "<< c1->all_nodes_map[node_name_2]->type<< " "<<endl;
  //   DoknownSim(c1, node_name, type, 3);
  //   return 0;
  // }
  // if(method==3){
  //   c1 = new_spec_cir.GetDuplicate("","","");
  //   int m = 4;
  //   string first_node_name;
  //   cin>>first_node_name;
  //   Node* n1 = c1->all_nodes_map[first_node_name];
  //   NodeVector nv1;
  //   GetPossibleNodes(c1, first_node_name, nv1, 0);
  //   string node_name_1 = first_node_name;
  //   string node_name_2;
  //   int n1_type;
  //   int n2_typ2;
  //   Node* n2;
  //   for(int i=0; i<nv1.size(); i++){
  //     vector<string> node_name;
  //     node_name_2 = nv1[i]->name;
  //     node_name.push_back(node_name_1);
  //     n2 = c1->all_nodes_map[node_name_2];
  //     node_name.push_back(node_name_2);
  //     for(int ii=0; ii<15; ii++){
  //       if(ii==n1->type) continue;
  //       for(int jj=0; jj<15; jj++){
  //         std::vector<int> type;
  //         if(jj==n2->type) continue;
  //         type.push_back(ii);
  //         type.push_back(jj);
  //         int k = DoknownSim(c1, node_name, type, 3);
  //         if(k==-1) cout<<node_name_1<<"  "<<node_name_2<<"   "<<ii <<"   "<<jj<<endl;
  //       }
  //     }
  //
  //   }
  //   // nv1.push_back(n1);
  //   // DoSimEqMulGate(c1, "test.blif", nv1,)
  //
  // }
  // if (method==2){
  //   ifstream inFile(filename+".csv", ios::in);
  //   ofstream outFile(filename+"info.csv", ios::app);
  //   if(!inFile){
  //     cout<<"Can't find this file!"<<endl;
  //     return 0;
  //   }
  //   string lineStr;
  //   Node* n1;
  //   c1 = new_spec_cir.GetDuplicate("","","");
  //   // long originalnum = nodenum(c1);
  //   // double all_input_pattern = pow(2, new_spec_cir.inputs.size());
  //   while (inFile>>lineStr)
  //   {
  //      c1 = new_spec_cir.GetDuplicate("","","");
  //      cout<<lineStr<<endl;
  //      stringstream ss;
  //      ss.str(lineStr);
  //      string str;
  //      string lineArray[4];
  //      int i=0;
  //      while (getline(ss, str, ',')){
  //        lineArray[i] = str;
  //        i++;
  //      }
  //      string node_name = lineArray[0];
  //      NodeVector nv1;
  //      n1 = c1->all_nodes_map[node_name];
  //      nv1.push_back(n1);
  //      c1->Simplify(nv1);
  //      c1->ResetLevels();
  //      c1->LevelizeSortTopological(false);
  //      c1->SetIndexes();
  //      outFile<<node_name<<","<<lineArray[2]<<","<<lineArray[3]<<","<<c1->outputs.size()<<endl;
  //      // long decreaseGate = originalnum - nodenum(c1);
  //      // if (stoi(lineArray[3])==0){
  //      //   DoknownSim(*c1, node_name, stoi(lineArray[2]), 3);
  //      //   DoSimEqLim(c1, spec_filename, node_name, stoi(lineArray[2]), 0, 5);
  //      // }
  //   }
  // }
  // if (method==4){
  //   ifstream inFile(filename+"G_0DIFF.csv", ios::in);
  //   ofstream outFile(filename+"G_0DIFF_info.csv", ios::app);
  //   if(!inFile){
  //     cout<<"Can't find this file!"<<endl;
  //     return 0;
  //   }
  //   string lineStr;
  //   Node* n1;
  //   Node* n2;
  //   c1 = new_spec_cir.GetDuplicate("","","");
  //   string str;
  //   string lineArray[7];
  //   string first_node_name;
  //   string second_node_name;
  //   while (inFile>>lineStr)
  //   {
  //      c1 = new_spec_cir.GetDuplicate("","","");
  //      stringstream ss;
  //      ss.str(lineStr);
  //      cout<<lineStr<<endl;
  //      int i=0;
  //      while (getline(ss, str, ',')){
  //        lineArray[i] = str;
  //        i++;
  //      }
  //      first_node_name = lineArray[0];
  //      second_node_name = lineArray[3];
  //      NodeVector nv1;
  //      n1 = c1->all_nodes_map[first_node_name];
  //      n2 = c1->all_nodes_map[first_node_name];
  //      nv1.push_back(n1);
  //      nv1.push_back(n2);
  //      c1->Simplify(nv1);
  //      c1->ResetLevels();
  //      c1->LevelizeSortTopological(false);
  //      c1->SetIndexes();
  //      n1 = c1->all_nodes_map[first_node_name];
  //      n2 = c1->all_nodes_map[second_node_name];
  //      outFile<<first_node_name<<","<<n1->level<<","<<second_node_name<<","<<n2->level<<","<<notfouts(n1->outputs, n2)<<","<<notfouts(n2->outputs, n1)
  //             <<","<<isGate(n1,n2)<<",";
  //      for(int i=0; i<n1->outputs.size(); i++){
  //        outFile<<n1->outputs[i]->name<<",";
  //      }
  //      outFile<<"----"<<",";
  //      for(int i=0; i<n2->outputs.size(); i++){
  //        outFile<<n2->outputs[i]->name<<",";
  //      }
  //      outFile<<endl;
  //   }
  // }
  // if (method==5){
  //   ifstream inFile(filename+"G_0pureAllKind.csv", ios::in);
  //   ofstream outFile(filename+"G_0_AllKind.csv", ios::app);
  //   if(!inFile){
  //     cout<<"Can't find this file!"<<endl;
  //     return 0;
  //   }
  //   string lineStr;
  //   Node* n1;
  //   Node* n2;
  //   c1 = new_spec_cir.GetDuplicate("","","");
  //   c2 = new_spec_cir.GetDuplicate("","","");
  //   string str;
  //   string lineArray[8];
  //   string first_node_name;
  //   string second_node_name;
  //   int firsttype;
  //   int secondtype;
  //   vector<bit64> po_diff_count;
  //   vector<bit64> total_diff_count;
  //   while (inFile>>lineStr)
  //   {
  //      c1 = new_spec_cir.GetDuplicate("","","");
  //      c2 = new_spec_cir.GetDuplicate("","","");
  //      stringstream ss;
  //      ss.str(lineStr);
  //      cout<<lineStr<<endl;
  //      int i=0;
  //      while (getline(ss, str, ',')){
  //        lineArray[i] = str;
  //        i++;
  //      }
  //      first_node_name = lineArray[0];
  //      second_node_name = lineArray[4];
  //      firsttype = stoi(lineArray[2]);
  //      secondtype = stoi(lineArray[6]);
  //      NodeVector nv1;
  //      n1 = c2->all_nodes_map[first_node_name];
  //      n2 = c2->all_nodes_map[second_node_name];
  //      nv1.push_back(n1);
  //      nv1.push_back(n2);
  //      c2->Simplify(nv1);
  //      c2->ResetLevels();
  //      c2->LevelizeSortTopological(false);
  //      c2->SetIndexes();
  //      NodeVector nv2;
  //      n1 = c1->all_nodes_map[first_node_name];
  //      n2 = c1->all_nodes_map[second_node_name];
  //      nv2.push_back(n1);
  //      nv2.push_back(n2);
  //      c1->Simplify(nv2);
  //      c1->ResetLevels();
  //      c1->LevelizeSortTopological(false);
  //      c1->SetIndexes();
  //      GetWrongCircuit(*c2, first_node_name, firsttype+1);
  //      GetWrongCircuit(*c2, second_node_name, secondtype+1);
  //      n1 = c2->all_nodes_map[first_node_name];
  //      n2 = c2->all_nodes_map[second_node_name];
  //      DoSimEq(*c1, *c2, po_diff_count, total_diff_count);
  //      outFile<<first_node_name<<","<<lineArray[1]<<","<<n1->type<<","<<second_node_name<<","<<lineArray[5]<<","<<n2->type<<","
  //      <<notfouts(n1->outputs, n2)<<","<<notfouts(n2->outputs, n1)<<","<<isGate(n1,n2)<<","<<total_diff_count[7]<<endl;
  //   }
  // }
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

int DoSimEqConnect(string spec_filename, int mode, int amount){
  Circuit new_spec_cir;
  string filename = spec_filename;
  size_t dotp = filename.find('.');
  ofstream outFile;
  ofstream outText;
  filename = filename.substr(0,dotp);
  // if(mode==1){
  //   outFile.open(filename+"conRes.csv", ios::out | ios::app);
  //   outText.open(filename+"conRes.txt", ios::out | ios::app);
  // }
  if(mode==0){
    // outFile.open(filename+"test.csv", ios::out | ios::app);
    // outText.open(filename+"test.txt", ios::out | ios::app);
    outFile.open(filename+"con.csv", ios::out | ios::app);
    outText.open(filename+"con.txt", ios::out | ios::app);
  }
  if(mode==2){
    outFile.open(filename+"ConRandom.csv", ios::out | ios::app);
    outText.open(filename+"ConRandom.txt", ios::out | ios::app);
  }
  // if(mode==3){
  //   outFile.open(filename+"conRadomAnaly.csv", ios::out | ios::app);
  //   outText.open(filename+"conRandomAnaly.txt", ios::out | ios::app);
  // }

  new_spec_cir.ReadBlif(spec_filename);
  if (new_spec_cir.inputs.size() == 0) {
    MSG("Cannot simulate circuit with no primary input!");
    return 0;
  }

  if (new_spec_cir.inputs.size() > 64) { // TODO: do we need to limit the inputs (because of high runtime?)
    MSG("Currently only simulating up to 64 inputs!");
    return 0;
  }


  // srand((unsigned)time( NULL));
  // int pattern_nums[5]={1000,5000,10000,20000,50000};

  if(mode==4){
    NodeVector select_nodes;
    Circuit* spec_cir = new_spec_cir.GetDuplicate("","","");
    NodeVector source_nodes;
    GetPossibleNodes(spec_cir, "", source_nodes, 0);
    RanSelectNodes(source_nodes, 100, select_nodes);
    for (int i=0; i<select_nodes.size(); i++){
      cout<<select_nodes[i]->name<<endl;
    }
    return 0;
  }

  if(amount==-1){
      int pattern_num = 10000;
      Circuit *this_spec_cir = new_spec_cir.GetDuplicate("","","");
      Circuit *pre_spec_cir = new_spec_cir.GetDuplicate("","","");
      NodeVector nodes;
      Node* this_node;
      Node* that_node;
      string this_node_name;
      string that_node_name;
      const int thread_num = 2;
      Circuit* spec_cirs[thread_num];
      Circuit* impl_cirs[thread_num];
      vector<bit64> po_diff_count[thread_num];
      vector<bit64> total_diff_count[thread_num];
      std::thread sim_threads[thread_num];
      GetPossibleNodes(pre_spec_cir, "", nodes, 0);
      NodeVector connections;
      NodeVector this_nodes;
      for (long i=0; i<nodes.size(); i++){
        this_node_name = nodes[i]->name;
        pre_spec_cir = new_spec_cir.GetDuplicate("","","");
        connections.clear();
        GetPossibleConnection(pre_spec_cir, this_node_name, connections);
        for (long j=0; j<connections.size(); j++){
          this_spec_cir = new_spec_cir.GetDuplicate("","","");
          that_node_name = connections[j]->name;
          this_node = this_spec_cir->all_nodes_map[this_node_name];
          that_node = this_spec_cir->all_nodes_map[that_node_name];
          this_nodes.clear();
          this_nodes.push_back(this_node);
          this_nodes.push_back(that_node);
          this_spec_cir->Simplify(this_nodes);
          this_spec_cir->SetIndexes();
          if (this_spec_cir->all_nodes_map.find(that_node_name) != this_spec_cir->all_nodes_map.end() and this_spec_cir->all_nodes_map.find(this_node_name) != this_spec_cir->all_nodes_map.end())  {
            spec_cirs[0] = this_spec_cir->GetDuplicate("","","");
            impl_cirs[0] = this_spec_cir->GetDuplicate("","","");
            GetWrongConnect(*impl_cirs[0], this_node_name, that_node_name, 0);
            if (mode==0) sim_threads[0] = std::thread(DoSimEqTh4, spec_cirs[0], impl_cirs[0], &po_diff_count[0], &total_diff_count[0]);
            spec_cirs[1] = this_spec_cir->GetDuplicate("","","");
            impl_cirs[1] = this_spec_cir->GetDuplicate("","","");
            GetWrongConnect(*impl_cirs[1], this_node_name, that_node_name, 1);
            if (mode==0) sim_threads[1] = std::thread(DoSimEqTh4, spec_cirs[1], impl_cirs[1], &po_diff_count[1], &total_diff_count[1]);
            if (mode==0){
              for(int k=0; k<2; k++){
                sim_threads[k].join();
              }
              for(int t=0; t<2; t++){
                double max_error = ceil(pow(2, impl_cirs[t]->inputs.size()) * 0.01);
                double needtoadd = pow(2, new_spec_cir.inputs.size() - this_spec_cir->inputs.size());
                if (total_diff_count[t][7] <= max_error){
                  outText << "----- changed node " << this_node->name << " ----"<< "connection " << that_node->name<<"---"<<"which input-----"<<t<<"----";
                  outText <<"wrong cases "<< total_diff_count[t][7]*needtoadd<<endl;
                  outFile << this_node->name<<","<<this_node->level<<","<< that_node->name<<","<<that_node->level<<","<<notfouts(that_node->outputs, this_node)<<","<<t<<","<<total_diff_count[t][7]*needtoadd <<endl;
                  for (int ii=0; ii<impl_cirs[t]->outputs.size(); ii++){
                    outText<<impl_cirs[t]->outputs[ii]->name<<"--";
                    cout<<endl;
                    outText<<endl;
                  }
                }
              }
            }
            if (mode==2){
              int h;
              for (int t=0; t<2; t++){
                cout << "------change node----- " << this_node_name << " ----connection name-----" << that_node_name << "----- "<<endl;
                for(int t=0; t<2; t++){
                  h = RndSimEq(spec_cirs[t], impl_cirs[t], pattern_num);
                  cout << "------change node----- " << this_node_name << " ----connection name-----" << that_node_name << "-----which input---- " << t << "---wrong case----" << h << endl;
                  outText << "------change node----- " << this_node_name << " ----connection name-----" << that_node_name << "-----which input---- " << t << "---wrong case----" << h << endl;
                  outFile << this_node_name << "," << that_node_name << "," << t << "," << notfouts(that_node->outputs, this_node) <<","<< h <<endl;
                }
              }
            }
          }
        }
      }
      return 0;
    }


  else if(amount>0){
    const int threads_num = 16;
    Circuit* spec_cirs[threads_num];
    Circuit* impl_cirs[threads_num];
    vector<bit64> po_diff_count[threads_num];
    vector<bit64> total_diff_count[threads_num];
    NodeVector outer_nodes;
    NodeVector inner_nodes;
    std::thread sim_threads[threads_num];
    int inner_amount = 20;
    Node* outer_node = NULL;
    Node* inner_node = NULL;
    Circuit* this_spec_cir;
    string outer_node_name;
    string inner_node_name;
    int cur_thread_num = 0;
    int pattern_num = 10000;
    NodeVector select_inner_nodes;
    NodeVector source_nodes;


    Circuit* spec_cir = new_spec_cir.GetDuplicate("","","");
    NodeVector total_nodes;
    GetPossibleNodes(spec_cir, "", total_nodes, 0);
    NodeVector select_outer_nodes;
    RanSelectNodes(total_nodes, amount, select_outer_nodes);
    for (int i=0; i<select_outer_nodes.size(); i++){
      outer_node_name = select_outer_nodes[i]->name;
      source_nodes.clear();
      select_inner_nodes.clear();
      inner_amount = 20;
      GetPossibleConnection(spec_cir, outer_node_name, source_nodes);
      if(source_nodes.size()<inner_amount) inner_amount = source_nodes.size();
      RanSelectNodes(source_nodes, inner_amount, select_inner_nodes);
      for (int j=0; j<select_inner_nodes.size(); j++){
        inner_node_name = select_inner_nodes[j]->name;
        this_spec_cir = new_spec_cir.GetDuplicate("","","");
        outer_node = this_spec_cir->all_nodes_map[outer_node_name];
        inner_node = this_spec_cir->all_nodes_map[inner_node_name];
        NodeVector this_nodes;
        this_nodes.push_back(outer_node);
        this_nodes.push_back(inner_node);
        this_spec_cir->Simplify(this_nodes);
        this_spec_cir->SetIndexes();
        if (this_spec_cir->all_nodes_map.find(outer_node_name) != this_spec_cir->all_nodes_map.end() and this_spec_cir->all_nodes_map.find(inner_node_name) != this_spec_cir->all_nodes_map.end())  {
          spec_cirs[cur_thread_num] = this_spec_cir->GetDuplicate("","","");
          impl_cirs[cur_thread_num] = this_spec_cir->GetDuplicate("","","");
          GetWrongConnect(*impl_cirs[cur_thread_num], outer_node_name, inner_node_name, 0);
          if(mode==0) sim_threads[cur_thread_num] = std::thread(DoSimEqTh4, spec_cirs[cur_thread_num], impl_cirs[cur_thread_num], &po_diff_count[cur_thread_num], &total_diff_count[cur_thread_num]);
          outer_nodes.push_back(outer_node);
          inner_nodes.push_back(inner_node);
          cur_thread_num++;
          spec_cirs[cur_thread_num] = this_spec_cir->GetDuplicate("","","");
          impl_cirs[cur_thread_num] = this_spec_cir->GetDuplicate("","","");
          GetWrongConnect(*impl_cirs[cur_thread_num], outer_node_name, inner_node_name, 1);
          if(mode==0) sim_threads[cur_thread_num] = std::thread(DoSimEqTh4, spec_cirs[cur_thread_num], impl_cirs[cur_thread_num], &po_diff_count[cur_thread_num], &total_diff_count[cur_thread_num]);
          cur_thread_num++;
        }
        if (mode==2){
          int h;
          cout << "------change node----- "<< outer_node_name << " ----connection name-----"<<inner_node_name<<"----- "<<endl;
          for(int t=0; t<2; t++){
            h=RndSimEq(spec_cirs[t], impl_cirs[t], pattern_num);
            cout << "------change node----- " << outer_node_name << " ----connection name-----" << inner_node_name << "-----which input---- " << t << "---wrong case----" << h << endl;
            outText << "------change node----- " << outer_node_name << " ----connection name-----" << inner_node_name << "-----which input---- " << t << "---wrong case----" << h << endl;
            outFile << outer_node_name << "," << inner_node_name << "," << t << "," << notfouts(inner_node->outputs, outer_node) <<","<< h <<endl;
          }
          cur_thread_num = 0;
        }
        if (mode == 0 and cur_thread_num == threads_num){
          for (int t=0; t<cur_thread_num; t++){
            sim_threads[t].join();
          }
          for (int t=0; t<cur_thread_num; t++){
            double needtoadd = pow(2, new_spec_cir.inputs.size() - spec_cirs[t]->inputs.size());
            double max_error = ceil(pow(2, spec_cirs[t]->inputs.size()) * 0.01);
            // if (total_diff_count[t][7] <= max_error){
            //   outFile << outer_nodes[floor(t/2)] << "," << inner_nodes[floor(t/2)] << "," << t%2 << "," << notfouts(inner_nodes[floor(t/2)]->outputs, outer_nodes[floor(t/2)]) <<","<< total_diff_count[t][7] <<endl;
            // }
            cout << " change node " << outer_nodes[floor(t/2)]->name << " connection name " << inner_nodes[floor(t/2)]->name << " which input " << t%2 << " wrong case " << total_diff_count[t][7]*needtoadd << endl;
            outFile << outer_nodes[floor(t/2)]->name << "," << inner_nodes[floor(t/2)]->name << "," << t%2 << "," << notfouts(inner_nodes[floor(t/2)]->outputs, outer_nodes[floor(t/2)]) <<","<< total_diff_count[t][7]*needtoadd <<endl;
          }
          cur_thread_num = 0;
          inner_nodes.clear();
          outer_nodes.clear();
        }
      }
    }
    return 0;
  }

  outFile.close();
  outText.close();
  return 0;
}

int DosimEqMulCon(Circuit &new_spec_cir, string spec_filename, NodeVector gate_nodes, int mode){
  ofstream outFile;
  ofstream outFile2;
  string filename = spec_filename;
  size_t dotp = filename.find('.');
  filename = filename.substr(0,dotp);
  if(gate_nodes.size()==1) {
    outFile.open("test_2C.csv",  ios::out | ios::app);
    outFile2.open("test_2C_0.csv", ios::out | ios::app);
  }

  else {
    outFile.open(filename+"_2C_type.csv",  ios::out | ios::app);
    outFile2.open(filename+"_2C_type_0.csv",  ios::out | ios::app);
  }

  // const int num_changes = 15;

  // Circuit* impl_cirs[num_changes];
  Circuit* origin_cir = new_spec_cir.GetDuplicate("", "", "");
  Circuit* impl_cir;
  Circuit* spec_cir;
  // std::thread sim_threads[num_changes];
  vector<bit64> po_diff_count;
  vector<bit64> total_diff_count;

  Node* first_node = new Node;
  Node* second_node = new Node;
  Node* gate_node = new Node;
  string first_node_name;
  string second_node_name;
  string gate_node_name;

  double needtoadd = 1;
  double max_error = 0;
  int cur_num = 0;
  NodeVector connection_nodes;
  NodeVector target_nodes;
  NodeVector tmp_nodes;
  int amount = 100;


  for (int i=0; i<gate_nodes.size(); i++){
    gate_node_name = gate_nodes[i]->name;
    connection_nodes.clear();
    tmp_nodes.clear();
    GetPossibleConnection(origin_cir, gate_node_name, tmp_nodes);
    RanSelectNodes(tmp_nodes, amount, connection_nodes);
    for (int j=0; j<connection_nodes.size(); j++){
      first_node_name = connection_nodes[j]->name;
      for (int jj=j+1; jj<connection_nodes.size(); jj++){
        second_node_name = connection_nodes[jj]->name;
        if (first_node_name == second_node_name) continue;
        spec_cir = new_spec_cir.GetDuplicate("","","");
        gate_node = spec_cir->all_nodes_map[gate_node_name];
        first_node = spec_cir->all_nodes_map[first_node_name];
        second_node = spec_cir->all_nodes_map[second_node_name];
        target_nodes.clear();
        target_nodes.push_back(gate_node);
        target_nodes.push_back(first_node);
        target_nodes.push_back(second_node);
        spec_cir->Simplify(target_nodes);
        spec_cir->SetIndexes();
        if (spec_cir->all_nodes_map.find(gate_node_name) != spec_cir->all_nodes_map.end() && spec_cir->all_nodes_map.find(first_node_name) != spec_cir->all_nodes_map.end() && spec_cir->all_nodes_map.find(second_node_name) != spec_cir->all_nodes_map.end()){
          impl_cir = spec_cir->GetDuplicate("","","");
          GetWrongConnect(*impl_cir, gate_node_name, first_node_name, 0);
          GetWrongConnect(*impl_cir, gate_node_name, second_node_name, 1);
          DoSimEq(*spec_cir, *impl_cir, po_diff_count, total_diff_count);
          cout << gate_node_name << "," << first_node_name << "," << second_node_name << "," << total_diff_count[7]*pow(2, new_spec_cir.inputs.size() - spec_cir->inputs.size()) << endl;
          if (total_diff_count[7]==0) outFile2 << gate_node_name << "," << first_node_name << "," << second_node_name << "," << 0 << endl;
          else outFile << gate_node_name << "," << first_node_name << "," << second_node_name << "," << total_diff_count[7]*pow(2, new_spec_cir.inputs.size() - spec_cir->inputs.size()) << endl;
        }
        delete spec_cir;
        // delete first_node;
        // delete second_node;
      }
    }
  }
  // delete origin_cir;
  return 0;
}

int DoMixedChange(Circuit &new_spec_cir, string spec_filename, int amount, int mode){
  ofstream outFile;
  string filename = spec_filename;
  size_t dotp = filename.find('.');
  filename = filename.substr(0,dotp);
  outFile.open(filename+"mixed_4test.csv",  ios::out | ios::app);

  srand((unsigned)time( NULL));

  const int num_changes = 15;

  if (mode==0){
    NodeVector gate_nodes;
    Node* gate_node = NULL;
    Circuit* ex_impl_cir;
    Circuit* ex_spec_cir;
    Circuit* ex_cir;
    Circuit* impl_cir;
    Circuit* spec_cir = new_spec_cir.GetDuplicate("","","");
    vector<bit64> po_diff_count;
    vector<bit64> total_diff_count;
    NodeVector simplify_nodes;
    string gate_node_name;
    string outer_node_name;
    string inner_node_name;
    Node* outer_node = NULL;
    Node* inner_node = NULL;
    NodeVector outer_nodes;
    NodeVector inner_nodes;
    NodeVector beitai_nodes;
    NodeVector all_nodes;
    GetPossibleNodes(spec_cir, "", all_nodes, 0);
    RanSelectNodes(all_nodes, 100, outer_nodes);

    for (int i=0; i<outer_nodes.size(); i++){
      impl_cir = new_spec_cir.GetDuplicate("","","");
      outer_node_name = outer_nodes[i]->name;
      GetPossibleConnection(impl_cir, outer_node_name, beitai_nodes);
      RanSelectNodes(beitai_nodes, 100, inner_nodes);
      for (int j=0; j<inner_nodes.size(); j++){
        inner_node_name = inner_nodes[j]->name;
        outer_node = impl_cir->all_nodes_map[outer_node_name];
        inner_node = impl_cir->all_nodes_map[inner_node_name];
        simplify_nodes.clear();
        simplify_nodes.push_back(outer_node);
        simplify_nodes.push_back(inner_node);
        impl_cir->Simplify(simplify_nodes);
        impl_cir->ResetLevels();
        impl_cir->LevelizeSortTopological(false);
        impl_cir->SetIndexes();
        if (impl_cir->all_nodes_map.find(outer_node_name) != impl_cir->all_nodes_map.end() and impl_cir->all_nodes_map.find(inner_node_name) != impl_cir->all_nodes_map.end()) {
          ex_impl_cir = impl_cir->GetDuplicate("","","");
          GetWrongConnect(*ex_impl_cir, outer_node_name, inner_node_name, rand()%1);
          GetPossibleNodes(impl_cir, "", beitai_nodes, 0);
          if (beitai_nodes.size()<100) RanSelectNodes(beitai_nodes, 100, gate_nodes);
          else RanSelectNodes(beitai_nodes, beitai_nodes.size(), gate_nodes);
          for (int k=0; k<gate_nodes.size(); k++){
            int type_change = rand()%16;
            gate_node_name = gate_nodes[k]->name;
            gate_node = ex_impl_cir->all_nodes_map[gate_node_name];
            while (type_change-1 == gate_node->type) type_change = rand()%16;
            GetWrongCircuit(*ex_impl_cir, gate_node_name, type_change);
            DoSimEq(*impl_cir, *ex_impl_cir, po_diff_count, total_diff_count);
            cout << "connection from " << inner_node_name <<" to gate " << outer_node_name << " gate change " << gate_node_name << " to type " << ex_impl_cir->all_nodes_map[gate_node_name]->type << " wrong is " <<total_diff_count[7] << endl;
            outFile << inner_node_name <<"," << outer_node_name << "," << gate_node_name << "," << ex_impl_cir->all_nodes_map[gate_node_name]->type << "," <<total_diff_count[7] << endl;
          }
        }
      }
    }
  }
  if (mode==1) {
    Circuit* origin_cir = new_spec_cir.GetDuplicate("","","");
    Circuit* spec_cir;
    Circuit* impl_cirs[num_changes];
    NodeVector gate_nodes;
    NodeVector connection_nodes;
    NodeVector simplify_nodes;
    std::vector<bit64> po_diff_count[num_changes];
    std::vector<bit64> total_diff_count[num_changes];
    std::thread sim_threads[num_changes];
    Node* gate_node = new Node;
    Node* connection_node = new Node;
    string gate_node_name;
    string connection_node_name;
    NodeVector tmp_nodes;
    int cur_num = 0;
    double needtoadd =  1;
    tmp_nodes.clear();
    GetPossibleNodes(origin_cir, "", tmp_nodes, 0);
    RanSelectNodes(tmp_nodes, amount, gate_nodes);
    for (int i=0; i<gate_nodes.size(); i++){
      gate_node_name = gate_nodes[i]->name;
      tmp_nodes.clear();
      connection_nodes.clear();
      GetPossibleConnection(origin_cir, gate_node_name, tmp_nodes);
      RanSelectNodes(tmp_nodes, amount, connection_nodes);
      for (int j=0; j<connection_nodes.size(); j++){
        connection_node_name = connection_nodes[j]->name;
        spec_cir = origin_cir->GetDuplicate("", "", "");
        gate_node = spec_cir->all_nodes_map[gate_node_name];
        connection_node = spec_cir->all_nodes_map[connection_node_name];
        simplify_nodes.clear();
        simplify_nodes.push_back(gate_node);
        simplify_nodes.push_back(connection_node);
        spec_cir->Simplify();
        spec_cir->SetIndexes();
        needtoadd = pow(2, new_spec_cir.inputs.size() - spec_cir->inputs.size());
        for (int mm=0; mm<2; mm++){
          cur_num = 0;
          for (int jj = 0; jj < num_changes+1; jj++){
            if (jj-1 == gate_node->type) continue;
            impl_cirs[cur_num] = spec_cir->GetDuplicate("","","");
            GetWrongCircuit(*impl_cirs[cur_num], gate_node_name, jj);
            GetWrongConnect(*impl_cirs[cur_num], gate_node_name, connection_node_name, mm);
            sim_threads[cur_num] = std::thread(DoSimEqTh4, spec_cir, impl_cirs[cur_num], &po_diff_count[cur_num], &total_diff_count[cur_num]);
            cur_num++;
          }
          for (int k = 0; k < num_changes; k++){
            sim_threads[k].join();
          }
          for (int k = 0; k < num_changes; k++){
            cout << gate_node_name << "," << gate_node->type << "," << impl_cirs[k]->all_nodes_map[gate_node_name]->type << "," << connection_node_name <<"," << total_diff_count[k][7] * needtoadd << endl;
            outFile << gate_node_name << "," << gate_node->type << "," << impl_cirs[k]->all_nodes_map[gate_node_name]->type << "," << connection_node_name <<"," << total_diff_count[k][7] * needtoadd;
            // if(checkPattern(total_diff_count[k])==0){
            //   for(int ii=0; ii<8; ii++){
            //     outFile << "," << total_diff_count[k][ii]*needtoadd;
            //   }
            // }
            outFile << endl;
          }
        }
        delete gate_node;
        delete connection_node;
      }
    }
  }
  return 0;
}

int DoSimEqLim(string spec_filename, int mode ,int amount){
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
      DoSimEqMulGate(new_spec_cir, spec_filename, first_nodes, 1);
      break;
    case 5:
      DoMixedChange(new_spec_cir, spec_filename, amount, 1);
      break;
    case 6:
      DosimEqMulCon(new_spec_cir, spec_filename, first_nodes, 0);
    case -1:
      new_spec_cir.WriteBlif(spec_filename);
      break;
    default:
      break;
  }
  // delete cur_node;
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

/*
dosimeqmodtwo: do simeq on gates to two circuits
*/
int DoSimEqModTwo(string spec_filename_1, string spec_filename_2, int mode, int amount){
  Circuit new_spec_cir_1, new_spec_cir_2;
  new_spec_cir_1.ReadBlif(spec_filename_1);
  new_spec_cir_2.ReadBlif(spec_filename_2);


  if (new_spec_cir_1.inputs.size() == 0 || new_spec_cir_2.inputs.size() == 0 ) {
    MSG("Cannot simulate circuit with no primary input!");
    return 0;
  }

  if (new_spec_cir_1.inputs.size() > 64 || new_spec_cir_2.inputs.size() > 64) { // TODO: do we need to limit the inputs (because of high runtime?)
    MSG("Currently only simulating up to 64 inputs!");
    return 0;
  }


  NodeVector target_nodes;
  Circuit* spec_cir;

  if (amount == -1){
    if(new_spec_cir_1.all_nodes.size() < new_spec_cir_2.all_nodes.size()){
      spec_cir = new_spec_cir_1.GetDuplicate("","","");
      GetPossibleNodes(spec_cir, "", target_nodes, 0);
      // for (int i=0; i<spec_cir->all_nodes.size(); i++){
      //   cur_node = spec_cir->all_nodes[i];
      //   if (cur_node->inputs.size()==2 && !cur_node->is_input && !cur_node->is_output){
      //     for (int j=0; j<impl_cir->all_nodes.size(); j++){
      //       that_node = impl_cir->all_nodes[j];
      //       if (that_node->inputs.size()==2 && !that_node->is_input && !that->is_output && issamenode(cur_node, that_node)){
      //         target_node_names_1.push_back(cur_node->name);
      //         target_node_names_2.push_back(that_node->name);
      //         break;
      //       }
      //       if (abs(that_node_name->level-cur_node->level)>2) continue;
      //     }
      //   }
      // }
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
      // target_nodes
    }else{
      spec_cir = new_spec_cir_2.GetDuplicate("","","");
      NodeVector all_nodes;
      GetPossibleNodes(spec_cir, "", all_nodes, 0);
      RanSelectNodes(all_nodes, amount, target_nodes);
    }
  }
  delete spec_cir;

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
    if (total_diff_count[7]>=pow(2, new_spec_cir_1.inputs.size()) * 0.01 ) {
      return 0;
    }
    outText << filename << "," << total_diff_count[7] << ",";
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
    int b;
    int count=0;

    for (int m=0; m<target_nodes.size(); m++){
      target_node_name = target_nodes[m]->name;
      spec_cir_1 = new_spec_cir_1.GetDuplicate("","","");
      spec_cir_2 = new_spec_cir_2.GetDuplicate("","","");
      if (spec_cir_1->all_nodes_map.find(target_node_name) == spec_cir_1->all_nodes_map.end() || spec_cir_2->all_nodes_map.find(target_node_name) == spec_cir_2->all_nodes_map.end()) continue;
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
          // long int mm = 0;
          // long int nn = 0;
          // for (int j = 0; j < 8; j++){
          //   if (j==0){
          //     mm = total_diff_count_1[t][0];
          //     nn = total_diff_count_2[t][0];
          //   }else{
          //     mm = total_diff_count_1[t][j] - total_diff_count_1[t][j-1];
          //     nn = total_diff_count_2[t][j] - total_diff_count_2[t][j-1] ;
          //   }
          //   outDetail << "\t" << mm << "\t" << nn;
          // }
          // outDetail << endl;
        }
      }
      delete spec_cir_1;
      delete spec_cir_2;
    }
    cout << "all have " << count <<endl;
    outText << count << endl;
    outText.close();
    outDetail.close();
  }

  if (mode==1) {
    ofstream outText;
    ofstream outDetail;
    string filename = spec_filename_2;
    size_t dotp = filename.find('.');
    filename = filename.substr(0,dotp);
    outText.open(filename+"condiff.csv", ios::out | ios::app);
    outDetail.open(filename+"condetail.csv", ios::out | ios::app);
    vector<bit64> po_diff_count;
    vector<bit64> total_diff_count;
    DoSimEq(new_spec_cir_1, new_spec_cir_2, po_diff_count, total_diff_count);
    outText << filename << "," << total_diff_count[7] << ",";
    outDetail << spec_filename_2 << "\t" << total_diff_count[7] << endl;
    for (int i=0; i< po_diff_count.size(); i++) {
      if (po_diff_count[i]!=0)
        outDetail << "port " << i << "\t" << po_diff_count[i] << endl;
    }
    outDetail << "detail" << endl;

    string target_node_name, connection_name;
    NodeVector connections;
    Circuit *spec_cir_1;
    Circuit *spec_cir_2;
    int doconnection = 20;
    int count = 0;
    std::thread sim_threads[doconnection];
    vector<bit64> total_diff_count_1[doconnection];
    vector<bit64> total_diff_count_2[doconnection];
    vector<string> conName;

    for (int m=0; m<target_nodes.size(); m++){
      target_node_name = target_nodes[m]->name;
      debug_msg(target_node_name);
      connections.clear();
      spec_cir_1 = new_spec_cir_1.GetDuplicate("","","");
      GetPossibleConnection(spec_cir_1, target_node_name, connections);
      if (spec_cir_1->all_nodes_map.find(target_node_name) == spec_cir_1->all_nodes_map.end()) continue;
      int nnsize = connections.size();
      int kksize = min(nnsize, doconnection);
      conName.clear();
      for (int n=0; n<kksize; n++) {
        spec_cir_1 = new_spec_cir_1.GetDuplicate("","","");
        spec_cir_2 = new_spec_cir_2.GetDuplicate("","","");
        connection_name = connections[n]->name;
        sim_threads[n] = std::thread(doConTwoTh4, spec_cir_1, spec_cir_2, &total_diff_count_1[n], &total_diff_count_2[n], target_node_name, connection_name);
        conName.push_back(connection_name);
      }
      for (int k=0; k<kksize; k++) {
        sim_threads[k].join();
      }
      for (int k=0; k<kksize; k++) {
        if (total_diff_count_1[k][7]!=total_diff_count_2[k][7]) {
          count ++;
          outDetail << target_node_name << "\t" << conName[k] << "\t" << total_diff_count_1[k][7] << "\t" << total_diff_count_2[k][7] << endl;
        }
        cout << target_node_name << "," << conName[k] << ", " << total_diff_count_1[k][7] << "," <<  total_diff_count_2[k][7]  << endl;

      }
      delete spec_cir_1;
      delete spec_cir_2;
        // DoConnectionTwo(*spec_cir_1, *spec_cir_2, total_diff_count, target_node_name, connection_name)
    }
    outText << count << endl;
    cout << "All have " << count << endl;
  }
  return 0;
}
