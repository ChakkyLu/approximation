#include "global.h"
#include <bitset>
#include <ctime>

#include <iomanip>

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

extern bit64 MyRand64();

class PatternType {
public:
  PatternType() { high_sa = -1.0; low_sa = -1.0; }
  PatternType(const PatternType &cp) { pis = cp.pis; ff_cur = cp.ff_cur; ff_nxt = cp.ff_nxt; high_sa = cp.high_sa; low_sa = cp.low_sa; all_sas = cp.all_sas; }
  ~PatternType() { pis.clear(); ff_cur.clear(); ff_nxt.clear(); all_sas.clear(); }

  ValVector      pis;
  ValVector      ff_cur;
  ValVector      ff_nxt;
  double         high_sa;
  double         low_sa;
  vector<double> all_sas;
};

typedef vector<PatternType> PatternVector;

typedef vector<PatternVector> PatternVector2D;

int SimActivity(Circuit &new_spec_cir, int min_act_level) {

  auto start_time = chrono::system_clock::now();

  MSG("Simulating...");
  std::srand(std::time(0));

  Val64Vector pi_vect;
  Val64Vector po_vect_spec;
  pi_vect.resize(new_spec_cir.inputs.size());
  po_vect_spec.resize(new_spec_cir.outputs.size());

  int pi_index = 0;
  while (pi_index < new_spec_cir.inputs.size()) {
    bit64 new_rnd = MyRand64();
    //cout << std::bitset<64>(new_rnd) << endl;
    // TODO: check it is different from previous simulation pattern(s)
    pi_vect[pi_index] = new_rnd;
    pi_index++;
  }

  Val64Vector all_vals_spec;
  new_spec_cir.Simulate(pi_vect, po_vect_spec, all_vals_spec);


  MSG("Counting...");
  Val64Vector sig_activities;
  sig_activities.resize(63, 0);
  bit64 sig_count = 0;
  for (int i = 0; i < new_spec_cir.all_nodes.size(); i++) {
    Node* cur_node = new_spec_cir.all_nodes[i];
    if (!cur_node->is_input && !cur_node->is_output && cur_node->index < all_vals_spec.size() && cur_node->type != NODE_ZERO && cur_node->type != NODE_ONE) {
      sig_count++;
      bit64 diff_pattern = all_vals_spec[cur_node->index];
      diff_pattern = diff_pattern ^ (diff_pattern >> 1);
      for (int k = 0; k < 63; k++) {
        if (diff_pattern%2 == 1) {
          sig_activities[k]++;
        }
        diff_pattern >>= 1;
      }
    }
  }


  auto end_time = chrono::system_clock::now();
  chrono::duration<double> diff_time = end_time - start_time;

  if (min_act_level < 0) {
    cout << "total signals: " << sig_count << endl;
    for (int k = 0; k < 63; k++) {
      cout << sig_activities[k] << " -> " << setprecision(2) << (double) sig_activities[k] / sig_count << endl;
    }
  }
  else {
    int index = 0;
    do {
      double activity = (double) sig_activities[index] / sig_count;
      cout << sig_activities[index] << " / " << sig_count << " -> " << setprecision(2) << activity << endl;
      if (activity >= (double)min_act_level/100) {
        MSG("Found the patterns! dumping to pattern files...");
        ofstream pat1file("pattern1");
        ofstream pat2file("pattern2");
        if (pat1file.good() && pat2file.good()) {
          for (int pi_k = 0; pi_k < pi_vect.size(); pi_k++) {
            bit64 val = pi_vect[pi_k];
            val >>= index;
            pat1file << val%2 ? '1' : '0';
            val >>= 1;
            pat2file << val%2 ? '1' : '0';
          }
          pat1file << endl;
          pat2file << endl;
          pat1file.close();
          pat2file.close();
        }
        else {
          ERR("Error writing to pattern files!");
        }
        index = 64;
      }
      else {
        index++;
      }
    } while (index < 63);
  }
  cout << endl << "Simulation and activity count time is: " << setprecision(3) << diff_time.count() << " s" << endl;

  return 0;
}

int SimActivity(string spec_filename) {
  Circuit new_spec_cir;
  new_spec_cir.ReadBlif(spec_filename);

  if (new_spec_cir.inputs.size() == 0) {
    MSG("Cannot simulate circuit with no primary input!");
    return 0;
  }

  return SimActivity(new_spec_cir, -1);
}

Circuit* GetDiffCircuit2(Circuit &org_circuit, NodeVector& diff_h_nodes, NodeVector& diff_l_nodes, bool simplify = true) {
  // create a new circuit
  Circuit *diff_cir = org_circuit.GetDuplicate("", "", "");
  int num_diff_nodes = diff_h_nodes.size()+diff_l_nodes.size();

  string name_prefix = "cp_";
  Circuit* cp_spec_cir = org_circuit.GetDuplicate(name_prefix, name_prefix, name_prefix);

  diff_cir->inputs.reserve(diff_cir->inputs.size()*2);
  diff_cir->inputs.insert(diff_cir->inputs.end(), cp_spec_cir->inputs.begin(), cp_spec_cir->inputs.end());
  diff_cir->all_nodes.reserve(diff_cir->all_nodes.size()*2);
  diff_cir->all_nodes.insert(diff_cir->all_nodes.end(), cp_spec_cir->all_nodes.begin(), cp_spec_cir->all_nodes.end());
  for (map<string, Node *>::iterator it = cp_spec_cir->all_nodes_map.begin(); it != cp_spec_cir->all_nodes_map.end(); ++it)
    diff_cir->all_nodes_map.insert(*it);

  // remove the previous outputs
  for (int k = 0; k < diff_cir->outputs.size(); k++) {
    diff_cir->outputs[k]->is_output = false;
    cp_spec_cir->outputs[k]->is_output = false;
  }
  diff_cir->outputs.clear();
  diff_cir->outputs.reserve(num_diff_nodes);
  diff_cir->all_nodes.reserve(diff_cir->all_nodes.size()+num_diff_nodes);

  int cur_node_index = 0;
  while (cur_node_index < diff_h_nodes.size()) {
    string node_name = diff_h_nodes[cur_node_index]->name;
    Node* cur_node = diff_cir->all_nodes_map[node_name];
    Node* cur_node2 = diff_cir->all_nodes_map[name_prefix+node_name];
    Node* new_out_node = new Node(NODE_XOR);
    new_out_node->name = "diff_out_"+node_name;
    new_out_node->is_output = true;
    new_out_node->inputs.push_back(cur_node);
    new_out_node->inputs.push_back(cur_node2);
    cur_node->outputs.push_back(new_out_node);
    cur_node2->outputs.push_back(new_out_node);
    diff_cir->outputs.push_back(new_out_node);
    diff_cir->all_nodes.push_back(new_out_node);
    diff_cir->all_nodes_map[new_out_node->name] = new_out_node;
    cur_node_index++;
  }
  cur_node_index = 0;
  while (cur_node_index < diff_l_nodes.size()) {
    string node_name = diff_l_nodes[cur_node_index]->name;
    Node* cur_node = diff_cir->all_nodes_map[node_name];
    Node* cur_node2 = diff_cir->all_nodes_map[name_prefix+node_name];
    Node* new_out_node = new Node(NODE_XOR);
    new_out_node->name = "diff_out_"+node_name;
    new_out_node->is_output = true;
    new_out_node->inputs.push_back(cur_node);
    new_out_node->inputs.push_back(cur_node2);
    cur_node->outputs.push_back(new_out_node);
    cur_node2->outputs.push_back(new_out_node);
    diff_cir->outputs.push_back(new_out_node);
    diff_cir->all_nodes.push_back(new_out_node);
    diff_cir->all_nodes_map[new_out_node->name] = new_out_node;
    cur_node_index++;
  }

  if (simplify)
    diff_cir->Simplify();

  return diff_cir;
}

Circuit* GetDiffCircuit(Circuit &org_circuit, NodeVector& diff_h_nodes, NodeVector& diff_l_nodes, ValVector &allsig_vals, bool simplify = true) {
  // create a new circuit
  Circuit *diff_cir = org_circuit.GetDuplicate("", "", "");
  int num_diff_nodes = diff_h_nodes.size()+diff_l_nodes.size();

  NodeVector ffs_inputs;
  ffs_inputs.reserve(diff_cir->ffs.size());
  // remove FFs and make the circuit combinational
  for (int i = 0; i < diff_cir->ffs.size(); i++) {
    Node* dff = diff_cir->ffs[i];
    Node* dff_in = dff->inputs[0];
    for (NodeVector::iterator it = dff_in->outputs.begin(); it != dff_in->outputs.end(); ++it) {
      if (*it == dff) {
        dff_in->outputs.erase(it);
        break;
      }
    }
    dff->inputs.clear();
    dff->type = NODE_UNKNOWN;
    dff->is_output = false;
    dff->is_input = true;
    ffs_inputs.push_back(dff_in);
  }
  if (diff_cir->ffs.size() > 0) {
    diff_cir->inputs.reserve(diff_cir->inputs.size() + diff_cir->ffs.size());
    diff_cir->inputs.insert(diff_cir->inputs.end(), diff_cir->ffs.begin(), diff_cir->ffs.end());
    diff_cir->ffs.clear();
  }

  // remove the previous outputs
  for (int k = 0; k < diff_cir->outputs.size(); k++) {
    diff_cir->outputs[k]->is_output = false;
  }
  diff_cir->outputs.clear();
  diff_cir->outputs.reserve(num_diff_nodes+org_circuit.ffs.size());
  diff_cir->all_nodes.reserve(diff_cir->all_nodes.size()+num_diff_nodes);

  int cur_node_index = 0;
  while (cur_node_index < diff_h_nodes.size()) {
    string node_name = diff_h_nodes[cur_node_index]->name;
    Node* cur_node = diff_cir->all_nodes_map[node_name];
    Node* new_node = NULL;
    if (allsig_vals[diff_h_nodes[cur_node_index]->index])
      new_node = new Node(NODE_NOT);
    else
      new_node = new Node(NODE_BUF);
    new_node->name = "diff_out_"+cur_node->name;
    new_node->is_output = true;
    new_node->inputs.push_back(cur_node);
    cur_node->outputs.push_back(new_node);
    diff_cir->outputs.push_back(new_node);
    diff_cir->all_nodes.push_back(new_node);
    diff_cir->all_nodes_map[new_node->name] = new_node;
    cur_node_index++;
  }
  cur_node_index = 0;
  while (cur_node_index < diff_l_nodes.size()) {
    string node_name = diff_l_nodes[cur_node_index]->name;
    Node* cur_node = diff_cir->all_nodes_map[node_name];
    Node* new_node;
    if (allsig_vals[diff_l_nodes[cur_node_index]->index])
      new_node = new Node(NODE_NOT);
    else
      new_node = new Node(NODE_BUF);
    new_node->name = "diff_out_"+cur_node->name;
    new_node->is_output = true;
    new_node->inputs.push_back(cur_node);
    cur_node->outputs.push_back(new_node);
    diff_cir->outputs.push_back(new_node);
    diff_cir->all_nodes.push_back(new_node);
    diff_cir->all_nodes_map[new_node->name] = new_node;
    cur_node_index++;
  }

  diff_cir->outputs.insert(diff_cir->outputs.end(), ffs_inputs.begin(), ffs_inputs.end());
  // TODO: ff inputs are assumed not to be in the output list! otherwise trouble?
  for (int i = 0; i < ffs_inputs.size(); i++)
    ffs_inputs[i]->is_output = true;

  diff_cir->WriteBlif("diff_pre_simp.blif");
  if (simplify) {
    diff_cir->Simplify();
  }
  return diff_cir;
}

Circuit* GetDiffCircuit(Circuit &org_circuit, NodeVector& diff_h_nodes, NodeVector& diff_l_nodes, bool simplify = true) {
  // create a new circuit
  Circuit *diff_cir = org_circuit.GetDuplicate("", "", "");
  Circuit *init_cir = org_circuit.GetDuplicate("cp_", "cp_", "cp_");
  int num_diff_nodes = diff_h_nodes.size()+diff_l_nodes.size();

  diff_cir->inputs.reserve(2*org_circuit.inputs.size()+org_circuit.ffs.size());
  diff_cir->outputs.reserve(org_circuit.ffs.size()+num_diff_nodes);
  diff_cir->all_nodes.reserve(2*org_circuit.all_nodes.size()+num_diff_nodes);

  diff_cir->inputs.insert(diff_cir->inputs.end(), init_cir->inputs.begin(), init_cir->inputs.end());
  diff_cir->all_nodes.insert(diff_cir->all_nodes.end(), init_cir->all_nodes.begin(), init_cir->all_nodes.end());
  diff_cir->all_nodes_map.merge(init_cir->all_nodes_map);

  NodeVector ffs_inputs;
  ffs_inputs.reserve(org_circuit.ffs.size());
  // remove FFs and make the circuit combinational
  for (int i = 0; i < org_circuit.ffs.size(); i++) {
    Node* dff2 = init_cir->ffs[i];
    Node* dff2_in = dff2->inputs[0];
    for (NodeVector::iterator it = dff2_in->outputs.begin(); it != dff2_in->outputs.end(); ++it) {
      if (*it == dff2) {
        dff2_in->outputs.erase(it);
        break;
      }
    }
    dff2->inputs.clear();
    dff2->type = NODE_UNKNOWN;
    dff2->is_output = false;
    dff2->is_input = true;

    diff_cir->inputs.push_back(dff2);

    Node* dff = diff_cir->ffs[i];
    Node* dff_in = dff->inputs[0];
    for (NodeVector::iterator it = dff_in->outputs.begin(); it != dff_in->outputs.end(); ++it) {
      if (*it == dff) {
        dff_in->outputs.erase(it);
        break;
      }
    }
    dff->inputs[0] = dff2_in;
    dff->type = NODE_BUF;
    dff->is_output = false;
    dff->is_input = false;

    dff_in->is_output = true;
    ffs_inputs.push_back(dff_in);
  }
  diff_cir->ffs.clear();
  init_cir->ffs.clear();


  // remove the previous outputs
  for (int k = 0; k < org_circuit.outputs.size(); k++) {
    diff_cir->outputs[k]->is_output = false;
    init_cir->outputs[k]->is_output = false;
  }
  diff_cir->outputs.clear();
  init_cir->outputs.clear();

  int cur_node_index = 0;
  while (cur_node_index < diff_h_nodes.size()) {
    string node_name = diff_h_nodes[cur_node_index]->name;
    Node* cur_node = diff_cir->all_nodes_map[node_name];
    Node* cur_node2 = diff_cir->all_nodes_map["cp_"+node_name];
    Node* new_node = new Node(NODE_XOR);
    new_node->name = "diff_out_"+cur_node->name;
    new_node->is_output = true;
    new_node->inputs.push_back(cur_node);
    new_node->inputs.push_back(cur_node2);
    cur_node->outputs.push_back(new_node);
    cur_node2->outputs.push_back(new_node);
    diff_cir->outputs.push_back(new_node);
    diff_cir->all_nodes.push_back(new_node);
    diff_cir->all_nodes_map[new_node->name] = new_node;
    cur_node_index++;
  }
  cur_node_index = 0;
  while (cur_node_index < diff_l_nodes.size()) {
    string node_name = diff_l_nodes[cur_node_index]->name;
    Node* cur_node = diff_cir->all_nodes_map[node_name];
    Node* cur_node2 = diff_cir->all_nodes_map["cp_"+node_name];
    Node* new_node = new Node(NODE_XOR);
    new_node->name = "diff_out_"+cur_node->name;
    new_node->is_output = true;
    new_node->inputs.push_back(cur_node);
    new_node->inputs.push_back(cur_node2);
    cur_node->outputs.push_back(new_node);
    cur_node2->outputs.push_back(new_node);
    diff_cir->outputs.push_back(new_node);
    diff_cir->all_nodes.push_back(new_node);
    diff_cir->all_nodes_map[new_node->name] = new_node;
    cur_node_index++;
  }

  diff_cir->outputs.insert(diff_cir->outputs.end(), ffs_inputs.begin(), ffs_inputs.end());
  // TODO: ff inputs are assumed not to be in the output list! otherwise trouble?

  diff_cir->WriteBlif("diff_pre_simp.blif");
  if (simplify) {
    diff_cir->ResetLevels();
    diff_cir->LevelizeSortTopological(false);
    diff_cir->RemoveBufNot();
    diff_cir->Simplify();
  }
  return diff_cir;
}

int AddZeroOneOutputs(Circuit* org_circuit, int num_zero, int num_one) {
  if (org_circuit == NULL || num_zero < 0 || num_one < 0)
    return 0;
  if (num_zero == 0 && num_one == 0)
    return 0;
  int org_output_size = org_circuit->outputs.size();
  org_circuit->outputs.resize(org_output_size+num_zero+num_one, NULL);

  int cur_index = org_circuit->outputs.size() - 1;
  for (int i = num_one; i > 0; ) {
    i--;
    Node* one_node = new Node(NODE_ONE);
    one_node->is_output = true;
    one_node->name = "diff_out_one_"+std::to_string(i);
    org_circuit->outputs[cur_index] = one_node;
    org_circuit->all_nodes.push_back(one_node);
    org_circuit->all_nodes_map[one_node->name] = one_node;
    cur_index--;
  }
  if (num_zero > 0) {
    for (int j = org_output_size; j > 0;) {
      j--;
      org_circuit->outputs[cur_index] = org_circuit->outputs[j];
      cur_index--;
    }
    for (int k = num_zero; k > 0; ) {
      k--;
      Node* zero_node = new Node(NODE_ZERO);
      zero_node->is_output = true;
      zero_node->name = "diff_out_zero_"+std::to_string(k);
      org_circuit->outputs[cur_index] = zero_node;
      org_circuit->all_nodes.push_back(zero_node);
      org_circuit->all_nodes_map[zero_node->name] = zero_node;
      cur_index--;
    }
  }

  return 0;
}

Circuit* CreateTopSatCircuit(Circuit* diff_cir, int num_cur_high_outputs, int num_cur_low_outputs, int num_final_high_outputs, int num_final_low_outputs, int target_high_sa, int target_low_sa) {
  if (num_final_high_outputs == 0 && num_final_low_outputs == 0)
    return NULL;
  if (num_final_high_outputs != 0 && num_final_high_outputs != 1024 && num_final_high_outputs != 4096 && num_final_high_outputs < 1024 && num_final_high_outputs%128 != 0)
    return NULL;
  if (num_final_low_outputs != 0 && num_final_low_outputs != 1024 && num_final_low_outputs != 4096 && num_final_low_outputs < 1024 && num_final_low_outputs%128 != 0)
    return NULL;

  int num_ffs = diff_cir->outputs.size() - num_cur_high_outputs - num_cur_low_outputs;
  //if (diff_cir->outputs.size() != num_cur_high_outputs + num_cur_low_outputs)
  //  return NULL;

  if (target_high_sa < 0 || target_high_sa > 100)
    return NULL;
  if (target_low_sa < 0 || target_low_sa > 100)
    return NULL;

  int num_rem = num_final_high_outputs+num_final_low_outputs-diff_cir->outputs.size()+num_ffs;
  if (num_rem < 0)
    return NULL;

  if (num_final_high_outputs < num_cur_high_outputs)
    return NULL;
  if (num_final_low_outputs < num_cur_low_outputs)
    return NULL;

  if (num_rem > 0) {
    int num_rem_high = num_final_high_outputs-num_cur_high_outputs;
    int num_one_high = 0;//(num_rem_high*target_high_sa)/100;
    int num_zero_high = num_rem_high-num_one_high;

    int num_rem_low = num_final_low_outputs-num_cur_low_outputs;
    int num_one_low = 0;//(num_rem_low*target_low_sa)/100;
    int num_zero_low = num_rem_low-num_one_low;

    diff_cir->outputs.reserve(num_final_high_outputs+num_final_low_outputs+num_ffs);
    if (num_final_high_outputs == 0)
      AddZeroOneOutputs(diff_cir, num_zero_low, num_one_low);
    else if (num_final_low_outputs == 0)
      AddZeroOneOutputs(diff_cir, num_zero_high, num_one_high);
    else {
      NodeVector temp_high(num_cur_high_outputs);
      NodeVector temp_low(num_cur_low_outputs);
      temp_high.insert(temp_high.end(), diff_cir->outputs.begin(), diff_cir->outputs.begin()+num_cur_high_outputs);
      temp_low.insert(temp_low.end(), diff_cir->outputs.begin()+num_cur_high_outputs, diff_cir->outputs.end());

      diff_cir->outputs = temp_low;
      AddZeroOneOutputs(diff_cir, num_zero_low, num_one_low);
      temp_low = diff_cir->outputs;
      diff_cir->outputs = temp_high;
      AddZeroOneOutputs(diff_cir, num_zero_high, num_one_high);
      //temp_high = diff_cir->outputs;
      diff_cir->outputs.insert(diff_cir->outputs.end(), temp_low.begin(), temp_low.end());
    }
  }

  diff_cir->WriteBlif("diff_act.blif");

  // write the top level model connecting the diff circuit and sorter
  ofstream top_cir_file("top.blif");
  if (top_cir_file.good()) {
    top_cir_file << ".model top" << endl;
    top_cir_file << ".inputs ";
    for (int i = 0; i < diff_cir->inputs.size(); i++)
      top_cir_file << "\\" << endl << diff_cir->inputs[i]->name;
    top_cir_file << endl << ".outputs ";
    if (num_final_high_outputs > 0) {
      for (int i = 0; i < 10; i++)
        top_cir_file << "\\" << endl << "yh0" << std::to_string(i);
      for (int i = 10; i < num_final_high_outputs; i++)
        top_cir_file << "\\" << endl << "yh" << std::to_string(i);
    }
    if (num_final_low_outputs > 0) {
      for (int i = 0; i < 10; i++)
        top_cir_file << "\\" << endl << "yl0" << std::to_string(i);
      for (int i = 10; i < num_final_low_outputs; i++)
        top_cir_file << "\\" << endl << "yl" << std::to_string(i);
    }
    for (int ff = num_final_high_outputs+num_final_low_outputs; ff < diff_cir->outputs.size(); ff++) {
      top_cir_file << "\\" << endl << diff_cir->outputs[ff]->name;
    }
    top_cir_file << endl << endl;
    top_cir_file << ".subckt " << diff_cir->name;
    for (int i = 0; i < diff_cir->inputs.size(); i++)
      top_cir_file << "\\" << endl << diff_cir->inputs[i]->name << "=" << diff_cir->inputs[i]->name;
    if (num_final_high_outputs > 0) {
      for (int i = 0; i < 10; i++)
        top_cir_file << "\\" << endl << diff_cir->outputs[i]->name << "=" << "xh0" << std::to_string(i);
      for (int i = 10; i < num_final_high_outputs; i++)
        top_cir_file << "\\" << endl << diff_cir->outputs[i]->name << "=" << "xh" << std::to_string(i);
    }
    if (num_final_low_outputs > 0) {
      for (int i = 0; i < 10; i++)
        top_cir_file << "\\" << endl << diff_cir->outputs[i + num_final_high_outputs]->name << "=" << "xl0" << std::to_string(i);
      for (int i = 10; i < num_final_low_outputs; i++)
        top_cir_file << "\\" << endl << diff_cir->outputs[i + num_final_high_outputs]->name << "=" << "xl" << std::to_string(i);
    }
    for (int ff = num_final_high_outputs+num_final_low_outputs; ff < diff_cir->outputs.size(); ff++) {
      top_cir_file << "\\" << endl << diff_cir->outputs[ff]->name << "=" << diff_cir->outputs[ff]->name;
    }
    if (num_final_high_outputs > 0) {
      top_cir_file << endl << ".subckt Sorter" << num_final_high_outputs;
      for (int i = 0; i < 10; i++)
        top_cir_file << "\\" << endl << "x0" << std::to_string(i) << "=" << "xh0" << std::to_string(i);
      for (int i = 10; i < num_final_high_outputs; i++)
        top_cir_file << "\\" << endl << "x" << std::to_string(i) << "=" << "xh" << std::to_string(i);
      for (int i = 0; i < 10; i++)
        top_cir_file << "\\" << endl << "y0" << std::to_string(i) << "=" << "yh0" << std::to_string(i);
      for (int i = 10; i < num_final_high_outputs; i++)
        top_cir_file << "\\" << endl << "y" << std::to_string(i) << "=" << "yh" << std::to_string(i);
      top_cir_file << endl;
    }
    if (num_final_low_outputs > 0) {
      top_cir_file << endl << ".subckt Sorter" << num_final_low_outputs;
      for (int i = 0; i < 10; i++)
        top_cir_file << "\\" << endl << "x0" << std::to_string(i) << "=" << "xl0" << std::to_string(i);
      for (int i = 10; i < num_final_low_outputs; i++)
        top_cir_file << "\\" << endl << "x" << std::to_string(i) << "=" << "xl" << std::to_string(i);
      for (int i = 0; i < 10; i++)
        top_cir_file << "\\" << endl << "y0" << std::to_string(i) << "=" << "yl0" << std::to_string(i);
      for (int i = 10; i < num_final_low_outputs; i++)
        top_cir_file << "\\" << endl << "y" << std::to_string(i) << "=" << "yl" << std::to_string(i);
      top_cir_file << endl;
    }
    top_cir_file << ".end" << endl << endl;
    top_cir_file.close();

    string cmd_str = "cat top.blif diff_act.blif";
    if (num_final_high_outputs == 1024 || num_final_low_outputs == 1024)
      cmd_str += " sorter1k.blif";
    if (num_final_high_outputs == 4096 || num_final_low_outputs == 4096)
      cmd_str += " sorter4k.blif";
    if (num_final_high_outputs < 1024 && num_final_high_outputs > 0)
      cmd_str += " sorter" + std::to_string(num_final_high_outputs) + ".blif";
    if (num_final_low_outputs < 1024 && num_final_low_outputs > 0 && num_final_low_outputs != num_final_high_outputs)
      cmd_str += " sorter" + std::to_string(num_final_low_outputs) + ".blif";
    cmd_str += " > top_all.blif";
    system(cmd_str.c_str());

    ABC::Abc_Frame_t *pAbc;
    pAbc = ABC::Abc_FrameGetGlobalFrame();
    //string command = "read top_all.blif; strash ; write temp.bench ; read temp.bench ; topo ; write new_top_all.blif";
    string command = "read top_all.blif; strash ; write new_top_all.blif";
    if (ABC::Cmd_CommandExecute(pAbc, command.c_str())) {
      ERR(string("ABC Command execution error: ") + command);
      return NULL;
    }

    Circuit *top_cir = new Circuit;
    top_cir->ReadBlif("new_top_all.blif", false, false, false);
    //top_cir->LevelizeSortTopological(false);
    top_cir->RemoveBufNot();
    top_cir->Simplify();
    //top_cir->LevelizeSortTopological(false);
    //top_cir->SetIndexes();

    return top_cir;
  }

  return NULL;
}


int SatSimIncActivity(Circuit &new_spec_cir, ValVector& allsig_vals, int min_high_act_level, int min_low_act_level, int max_high_active_outputs = 1024, int max_low_active_outputs = 1024, int counter = 0) {
  const int num_sat_solutions = 3;
  ValVector allsig_vect[num_sat_solutions];

  bool do_high_sa = (max_high_active_outputs > 0);
  bool do_low_sa = (max_low_active_outputs > 0);

  if (!do_high_sa && !do_low_sa) {
    return 0;
  }
  if (do_high_sa && min_high_act_level <= min_low_act_level) {
    ERR("high activity should be greater than low activity!");
    return 0;
  }

  if (new_spec_cir.all_nodes.size()-new_spec_cir.inputs.size()-new_spec_cir.outputs.size() < max_high_active_outputs+max_low_active_outputs) {
    ERR("High/Low activity only works for circuits with more than 2k nodes!");
    return 0;
  }

  cout << "preparing files...." << endl;

  NodeVector sa_high_nodes, sa_low_nodes;

  int index = new_spec_cir.all_nodes.size() - max_high_active_outputs - max_low_active_outputs;
  if (do_low_sa) {
    sa_low_nodes.reserve(max_low_active_outputs);
    while (index < new_spec_cir.all_nodes.size() - max_high_active_outputs) {
      sa_low_nodes.push_back(new_spec_cir.all_nodes[index]);
      index++;
    }
  }
  if (do_high_sa) {
    sa_high_nodes.reserve(max_high_active_outputs);
    while (index < new_spec_cir.all_nodes.size()) {
      sa_high_nodes.push_back(new_spec_cir.all_nodes[index]);
      index++;
    }
  }

  // create a new circuit
  Circuit *diff_cir = GetDiffCircuit(new_spec_cir, sa_high_nodes, sa_low_nodes, allsig_vals);
  Circuit *top_cir = CreateTopSatCircuit(diff_cir, sa_high_nodes.size(), sa_low_nodes.size(), max_high_active_outputs, max_low_active_outputs, min_high_act_level, min_low_act_level);

  if (top_cir) {
    cout << "preparing for SAT..." <<endl;

    auto start_time = chrono::system_clock::now();

    GInterface sat_act_gen;

    Cnf top_cnf;
    top_cnf.Convert2Cnf(top_cir, NULL, NULL, false); // do not optimize as it is very time consuming, and useless!

    int h_low_i_index = (max_high_active_outputs*(100-min_high_act_level)*9)/1000;
    int h_high_i_index = (max_high_active_outputs*(100-min_high_act_level)*11)/1000;
    int l_low_i_index = (max_low_active_outputs*(100-min_low_act_level)*9)/1000;
    int l_high_i_index = (max_low_active_outputs*(100-min_low_act_level)*12)/1000;
    int min_low_margin = max_low_active_outputs/64;
    int min_high_margin = max_high_active_outputs/64;
    if (l_low_i_index < min_low_margin)
      l_low_i_index = min_low_margin;
    if (l_high_i_index > max_low_active_outputs-min_low_margin)
      l_high_i_index = max_low_active_outputs-min_low_margin;
    if (h_low_i_index < min_high_margin)
      h_low_i_index = min_high_margin;
    if (h_high_i_index > max_high_active_outputs-min_high_margin)
      h_high_i_index = max_high_active_outputs-min_high_margin;

    int cur_i = 0;
    if (do_high_sa) {
      while (cur_i < 10) {
        LiteralType lit_i = top_cnf.literal_str_map["yh0" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(-lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
      while (cur_i < h_low_i_index) {
        LiteralType lit_i = top_cnf.literal_str_map["yh" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(-lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
      cur_i = h_high_i_index;
      while (cur_i < max_high_active_outputs) {
        LiteralType lit_i = top_cnf.literal_str_map["yh" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
    }
    if (do_low_sa) {
      cur_i = 0;
      while (cur_i < 10) {
        LiteralType lit_i = top_cnf.literal_str_map["yl0" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(-lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
      while (cur_i < l_low_i_index) {
        LiteralType lit_i = top_cnf.literal_str_map["yl" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(-lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
      cur_i = l_high_i_index;
      while (cur_i < max_low_active_outputs) {
        LiteralType lit_i = top_cnf.literal_str_map["yl" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
    }


    top_cnf.WriteDimacs("temp.cnf");
    sat_act_gen.AddCnf(top_cnf.clauses, 0);

    for (int sol_num = 0; sol_num < num_sat_solutions; sol_num++) {
      if (do_high_sa)
        cout << "High activity level range is: " << 100-((100*h_high_i_index)/max_high_active_outputs) << "~" << 100-((100*h_low_i_index)/max_high_active_outputs) << " %" << endl;
      if (do_low_sa)
        cout << "Low activity level range is: " << 100-((100*l_high_i_index)/max_low_active_outputs) << "~" << 100-((100*l_low_i_index)/max_low_active_outputs) << " %" << endl;


      bool sat_res = sat_act_gen.GetAnswer();
      if (!sat_res) {
        if (sol_num == 0) {
          MSG("No Answer!");
          return -1;
        }
        MSG("No more answers!");
        break;
      }
      MSG("Found the patterns! dumping to pattern files...");

      if (do_high_sa) {
        ClauseVector blocks_clauses;
        int cur_i = h_high_i_index;
        bool sat_res_var = false;
        do {
          LiteralType lit_i = top_cnf.literal_str_map["yh" + std::to_string(cur_i)];
          sat_res_var = sat_act_gen.ModelNum(lit_i) == l__True;
          LiteralVector clause_i;
          clause_i.push_back(lit_i);
          blocks_clauses.push_back(clause_i);
          cur_i--;
        } while (sat_res_var);
        sat_act_gen.AddCnf(blocks_clauses, top_cnf.clauses.size());
        top_cnf.clauses.insert(top_cnf.clauses.end(), blocks_clauses.begin(), blocks_clauses.end());
        h_high_i_index = cur_i;
      }
      if (do_low_sa) {
        ClauseVector blocks_clauses;
        int cur_i = l_low_i_index;
        bool sat_res_var = false;
        do {
          LiteralType lit_i = top_cnf.literal_str_map["yl" + std::to_string(cur_i)];
          sat_res_var = sat_act_gen.ModelNum(lit_i) == l__True;
          LiteralVector clause_i;
          clause_i.push_back(-lit_i);
          blocks_clauses.push_back(clause_i);
          cur_i++;
        } while (!sat_res_var);
        sat_act_gen.AddCnf(blocks_clauses, top_cnf.clauses.size());
        top_cnf.clauses.insert(top_cnf.clauses.end(), blocks_clauses.begin(), blocks_clauses.end());
        l_low_i_index = cur_i;
      }

      ValVector pi_vect(new_spec_cir.inputs.size(), false);
      ValVector po_vect(new_spec_cir.outputs.size(), false);

      ofstream patfile("patterns.csv", ofstream::out|ofstream::app);
      if (patfile.good()) {
        patfile << counter << "," << sol_num << ",";
        for (int pi_k = 0; pi_k < new_spec_cir.inputs.size(); pi_k++) {
          LiteralType lit_block = top_cnf.GetLiteral(new_spec_cir.inputs[pi_k]->name);
          bool sat_res_var1 = sat_act_gen.ModelNum(lit_block) == l__True;
          patfile << sat_res_var1 ? '1' : '0';
          if (lit_block <= 0)
            patfile << '?';
          pi_vect[pi_k] = sat_res_var1;
        }
        patfile << endl;
        patfile.close();
      }
      else {
        ERR("Error opening the pattern file to write the results!");
        return -1;
      }

      allsig_vect[sol_num].resize(new_spec_cir.all_nodes.size(), false);
      new_spec_cir.Simulate(pi_vect, po_vect, allsig_vect[sol_num]);


      int num_diff_nodes = 0;
      bit64 diff_count = 0;
      for (int i = 0; i < allsig_vect[sol_num].size(); i++) {
        Node *cur_node = new_spec_cir.all_nodes[i];
        if (!cur_node->is_input && cur_node->type != NODE_ZERO && cur_node->type != NODE_ONE) {
          num_diff_nodes++;
          if (allsig_vect[sol_num][i] != allsig_vals[i])
            diff_count++;
        }
      }
      bit64 diff_count_h = 0;
      for (int i = 0; i < sa_high_nodes.size(); i++) {
        index = sa_high_nodes[i]->index;
        if (allsig_vect[sol_num][index] != allsig_vals[index])
          diff_count_h++;
      }
      bit64 diff_count_l = 0;
      for (int i = 0; i < sa_low_nodes.size(); i++) {
        index = sa_low_nodes[i]->index;
        if (allsig_vect[sol_num][index] != allsig_vals[index])
          diff_count_l++;
      }
      cout << "Activity rate is: " << setprecision(2) << (double) diff_count / num_diff_nodes << " (" << num_diff_nodes << ") :::";
      if (do_high_sa)
        cout << " High=" << (double) diff_count_h / max_high_active_outputs;
      if (do_low_sa)
        cout << " Low=" << (double) diff_count_l / max_low_active_outputs;
      cout << endl;
    }

    auto end_time = chrono::system_clock::now();
    chrono::duration<double> diff_time = end_time - start_time;
    cout << endl << "Pattern creation by SAT time is: " << setprecision(3) << diff_time.count() << " s" << endl;
  }
  else {
    ERR("Cannot write top level module!");
  }

  delete top_cir;
  delete diff_cir;

  if (counter > 0) {
    int index = num_sat_solutions-1;
    while (allsig_vect[index].size() == 0)
      index--;
    while (index > 0) {
      int ret = SatSimIncActivity(new_spec_cir, allsig_vect[index], min_high_act_level, min_low_act_level, max_high_active_outputs, max_low_active_outputs, counter - 1);
      if (ret >= 0)
        return 0;
      index--;
    }
    return -1;
  }
  return 0;
}

int SatSatIncActivity(Circuit &new_spec_cir, int min_high_act_level, int min_low_act_level, int max_high_active_outputs = 1024, int max_low_active_outputs = 1024) {
  if (max_high_active_outputs == 0 && max_low_active_outputs == 0) {
    return 0;
  }
  if (min_high_act_level <= min_low_act_level) {
    ERR("high activity should be greater than low activity!");
    return 0;
  }

  if (new_spec_cir.all_nodes.size()-new_spec_cir.inputs.size()-new_spec_cir.outputs.size() < max_high_active_outputs+max_low_active_outputs) {
    ERR("High/Low activity only works for circuits with more than 2k nodes!");
    return 0;
  }

  bool do_high_sa = (max_high_active_outputs > 0);
  bool do_low_sa = (max_low_active_outputs > 0);

  cout << "preparing files...." << endl;

  NodeVector sa_high_nodes, sa_low_nodes;

  int index = new_spec_cir.all_nodes.size() - max_high_active_outputs - max_low_active_outputs;
  if (do_low_sa) {
    sa_low_nodes.reserve(max_low_active_outputs);
    while (index < new_spec_cir.all_nodes.size() - max_high_active_outputs) {
      sa_low_nodes.push_back(new_spec_cir.all_nodes[index]);
      index++;
    }
  }
  if (do_high_sa) {
    sa_high_nodes.reserve(max_high_active_outputs);
    while (index < new_spec_cir.all_nodes.size()) {
      sa_high_nodes.push_back(new_spec_cir.all_nodes[index]);
      index++;
    }
  }

  // create a new circuit
  Circuit *diff_cir = GetDiffCircuit2(new_spec_cir, sa_high_nodes, sa_low_nodes);
  Circuit *top_cir = CreateTopSatCircuit(diff_cir, sa_high_nodes.size(), sa_low_nodes.size(), max_high_active_outputs, max_low_active_outputs, min_high_act_level, min_low_act_level);

  if (top_cir) {
    cout << "preparing for SAT..." <<endl;

    auto start_time = chrono::system_clock::now();

    GInterface sat_act_gen;

    Cnf top_cnf;
    top_cnf.Convert2Cnf(top_cir, NULL, NULL, false); // do not optimize as it is very time consuming, and useless!

    int h_low_i_index = (max_high_active_outputs*(100-min_high_act_level)*9)/1000;
    int h_high_i_index = (max_high_active_outputs*(100-min_high_act_level)*11)/1000;
    int l_low_i_index = (max_low_active_outputs*(100-min_low_act_level)*9)/1000;
    int l_high_i_index = (max_low_active_outputs*(100-min_low_act_level)*12)/1000;
    int min_low_margin = max_low_active_outputs/64;
    int min_high_margin = max_high_active_outputs/64;
    if (l_low_i_index < min_low_margin)
      l_low_i_index = min_low_margin;
    if (l_high_i_index > max_low_active_outputs-min_low_margin)
      l_high_i_index = max_low_active_outputs-min_low_margin;
    if (h_low_i_index < min_high_margin)
      h_low_i_index = min_high_margin;
    if (h_high_i_index > max_high_active_outputs-min_high_margin)
      h_high_i_index = max_high_active_outputs-min_high_margin;
    if (do_high_sa)
      cout << "High activity level range is: " << 100-((100*h_high_i_index)/max_high_active_outputs) << "~" << 100-((100*h_low_i_index)/max_high_active_outputs) << " %" << endl;
    if (do_low_sa)
      cout << "Low activity level range is: " << 100-((100*l_high_i_index)/max_low_active_outputs) << "~" << 100-((100*l_low_i_index)/max_low_active_outputs) << " %" << endl;

    int cur_i = 0;
    if (do_high_sa) {
      while (cur_i < 10) {
        LiteralType lit_i = top_cnf.literal_str_map["yh0" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(-lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
      while (cur_i < h_low_i_index) {
        LiteralType lit_i = top_cnf.literal_str_map["yh" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(-lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
      cur_i = h_high_i_index;
      while (cur_i < max_high_active_outputs) {
        LiteralType lit_i = top_cnf.literal_str_map["yh" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
    }
    if (do_low_sa) {
      cur_i = 0;
      while (cur_i < 10) {
        LiteralType lit_i = top_cnf.literal_str_map["yl0" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(-lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
      while (cur_i < l_low_i_index) {
        LiteralType lit_i = top_cnf.literal_str_map["yl" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(-lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
      cur_i = l_high_i_index;
      while (cur_i < max_low_active_outputs) {
        LiteralType lit_i = top_cnf.literal_str_map["yl" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
    }


    top_cnf.WriteDimacs("temp.cnf");
    sat_act_gen.AddCnf(top_cnf.clauses, 0);

    bool sat_res = sat_act_gen.GetAnswer();
    if (!sat_res) {
      MSG("No Answer!");
      return -1;
    }
    MSG("Found the patterns! dumping to pattern files...");

    ValVector pi_vect1(new_spec_cir.inputs.size(), false);
    ValVector pi_vect2(new_spec_cir.inputs.size(), false);
    ValVector po_vect(new_spec_cir.outputs.size(), false);
    ValVector allsig_vect1(new_spec_cir.all_nodes.size(), false);
    ValVector allsig_vect2(new_spec_cir.all_nodes.size(), false);

    ofstream pat1file("pattern1");
    ofstream pat2file("pattern2");
    if (pat1file.good() && pat2file.good()) {
      for (int pi_k = 0; pi_k < new_spec_cir.inputs.size(); pi_k++) {
        LiteralType lit1 = top_cnf.GetLiteral(new_spec_cir.inputs[pi_k]->name);
        bool sat_res_var1 = sat_act_gen.ModelNum(lit1) == l__True;
        pat1file << sat_res_var1 ? '1' : '0';
        if (lit1 <= 0)
          pat1file << '?';
        pi_vect1[pi_k] = sat_res_var1;
        LiteralType lit2 = top_cnf.GetLiteral("cp_"+new_spec_cir.inputs[pi_k]->name);
        bool sat_res_var2 = sat_act_gen.ModelNum(lit2) == l__True;
        pat2file << sat_res_var2 ? '1' : '0';
        pi_vect2[pi_k] = sat_res_var2;
        if (lit2 <= 0)
          pat2file << '?';
      }
      pat1file << endl;
      pat2file << endl;
      pat1file.close();
      pat2file.close();
    }

    new_spec_cir.Simulate(pi_vect1, po_vect, allsig_vect1);
    new_spec_cir.Simulate(pi_vect2, po_vect, allsig_vect2);


    int num_diff_nodes = 0;
    bit64 diff_count = 0;
    for (int i = 0; i < allsig_vect1.size(); i++) {
      Node* cur_node = new_spec_cir.all_nodes[i];
      if (!cur_node->is_input && cur_node->type != NODE_ZERO && cur_node->type != NODE_ONE) {
        num_diff_nodes++;
        if (allsig_vect1[i] != allsig_vect2[i])
          diff_count++;
      }
    }
    bit64 diff_count_h = 0;
    for (int i = 0; i < sa_high_nodes.size(); i++) {
      index = sa_high_nodes[i]->index;
      if (allsig_vect1[index] != allsig_vect2[index])
        diff_count_h++;
    }
    bit64 diff_count_l = 0;
    for (int i = 0; i < sa_low_nodes.size(); i++) {
      index = sa_low_nodes[i]->index;
      if (allsig_vect1[index] != allsig_vect2[index])
        diff_count_l++;
    }
    cout << "Activity rate is: " << setprecision(2) << (double)diff_count/num_diff_nodes << " (" << num_diff_nodes <<") :::";
    if (do_high_sa)
      cout << " High=" << (double)diff_count_h/max_high_active_outputs;
    if (do_low_sa)
      cout << " Low=" << (double)diff_count_l/max_low_active_outputs;
    cout << endl;

    auto end_time = chrono::system_clock::now();
    chrono::duration<double> diff_time = end_time - start_time;
    cout << endl << "Pattern creation by SAT time is: " << setprecision(3) << diff_time.count() << " s" << endl;

  }
  else {
    ERR("Cannot write top level module!");
  }

  return 0;
}

int SatSatIncActivity1kHL(Circuit &new_spec_cir, int min_act_level) {
  return SatSatIncActivity(new_spec_cir, min_act_level, min_act_level/2, 1024, 1024);
}

int SatSatIncActivity1k(Circuit &new_spec_cir, int min_act_level) {
  return SatSatIncActivity(new_spec_cir, min_act_level, 0, 1024, 0);
}

int SatSatIncActivity4k(Circuit &new_spec_cir, int min_act_level) {
  return SatSatIncActivity(new_spec_cir, min_act_level, 0, 4096, 0);
}

int SatSimIncActivity(Circuit &new_spec_cir, int high_act_level, int low_act_level, int num_cycles) {

  ValVector pi_vect(new_spec_cir.inputs.size(), false);
  ValVector po_vect(new_spec_cir.outputs.size(), false);
  ValVector allsig_vect(new_spec_cir.all_nodes.size(), false);

  ofstream patfile("patterns.csv", ofstream::out);
  if (patfile.good()) {
    patfile << num_cycles << ",0,";
    std::srand(std::time(0));
    for (int k = 0; k < pi_vect.size(); k++) {
      pi_vect[k] = (std::rand() % 2 == 1);
      patfile << pi_vect[k];
    }
    patfile << endl;
    patfile.close();
  }
  else {
    ERR("cannot open pattern file for writing the results!");
    return 0;
  }

  new_spec_cir.Simulate(pi_vect, po_vect, allsig_vect);

  return SatSimIncActivity(new_spec_cir, allsig_vect, high_act_level > 0 ? high_act_level : 99, low_act_level, high_act_level > 0 ? 1024 : 0, low_act_level > 0 ? 1024 : 0, num_cycles-1);
}

int SatIncActivity(string spec_filename, int high_act_level, int low_act_level, int num_cycles) {

  if (high_act_level < 1 || high_act_level > 99) {
    ERR("Activity level should be between 1 and 99 (percent)");
    return -1;
  }
  if (low_act_level < 1 || low_act_level > 99) {
    ERR("Activity level should be between 1 and 99 (percent)");
    return -1;
  }

  Circuit new_spec_cir;

  new_spec_cir.ReadBlif(spec_filename);


  //return SatSatIncActivity1kHL(new_spec_cir, min_act_level);
  //return SatSatIncActivity1k(new_spec_cir, min_act_level);
  //return SatSatIncActivity4k(new_spec_cir, min_act_level);
  return SatSimIncActivity(new_spec_cir, high_act_level, low_act_level, num_cycles);

  //return SatSimIncActivity(new_spec_cir, min_act_level);

  //if (new_spec_cir.all_nodes.size()-new_spec_cir.inputs.size() > 4096)
  //  return SatSatIncActivity4k(new_spec_cir, min_act_level);
  //else
  //  return SatSatIncActivity1k(new_spec_cir, min_act_level);
}



int SatSimIncActivity(Circuit &new_spec_cir, vector<NodeVector> &all_regions, NodeVector &sa_high_nodes, NodeVector &sa_low_nodes, ValVector& allsig_vals, int ref_high_act_level, int ref_low_act_level, int min_high_act_level, int min_low_act_level, int max_high_active_outputs, int max_low_active_outputs, PatternVector2D &patterns, int counter = 0) {
  const int num_sat_solutions = 5;
  ValVector allsig_vect[num_sat_solutions];
  PatternVector cur_patterns;
  PatternVector2D nxt_patterns;
  patterns.clear();

  bool do_high_sa = (max_high_active_outputs > 0);
  bool do_low_sa = (max_low_active_outputs > 0);

  if (!do_high_sa && !do_low_sa) {
    return -1;
  }
  if (do_high_sa && min_high_act_level <= min_low_act_level) {
    ERR("high activity should be greater than low activity!");
    return -1;
  }

  if (!do_low_sa)
    min_low_act_level = ref_low_act_level;
  if (!do_high_sa)
    min_high_act_level = ref_high_act_level;

  if (do_high_sa && min_high_act_level < ref_high_act_level*75/100)
    return -1;
  if (do_low_sa && min_low_act_level > ref_low_act_level*130/100)
    return -1;

  cout << "preparing files...." << endl;

  int num_high_nodes = sa_high_nodes.size();
  if (num_high_nodes > 4096) {
    sa_high_nodes.resize(4096);
    MSG("Cutting the high SA nodes to: 4096");
  }
  else if (num_high_nodes == 4096) {
  }
  else if (num_high_nodes > 1024) {
    sa_high_nodes.resize(1024);
    MSG("Cutting the high SA nodes to: 1024");
  }
  /*else if (num_high_nodes%128) {
    // scale the high SA
    min_high_act_level = (min_high_act_level*num_high_nodes)/(num_high_nodes+128-num_high_nodes%128);
    MSG("scaled high SA to: "+std::to_string(min_high_act_level));
  }*/
  int num_low_nodes = sa_low_nodes.size();
  if (num_low_nodes > 4096) {
    sa_low_nodes.resize(4096);
    MSG("Cutting the low SA nodes to: 4096");
  }
  else if (num_low_nodes == 4096) {
  }
  else if (num_low_nodes > 1024) {
    sa_low_nodes.resize(1024);
    MSG("Cutting the low SA nodes to: 1024");
  }
  /*else if (num_low_nodes%128) {
    // scale the low SA
    min_low_act_level = (min_low_act_level*num_low_nodes)/(num_low_nodes+128-num_low_nodes%128);
    MSG("scaled low SA to: "+std::to_string(min_low_act_level));
  }*/

  // create a new circuit
  Circuit *diff_cir = GetDiffCircuit(new_spec_cir, sa_high_nodes, sa_low_nodes, allsig_vals);
  Circuit *top_cir = CreateTopSatCircuit(diff_cir, sa_high_nodes.size(), sa_low_nodes.size(), max_high_active_outputs, max_low_active_outputs, min_high_act_level, min_low_act_level);

  if (top_cir) {
    cout << "preparing for SAT..." <<endl;

    auto start_time = chrono::system_clock::now();

    GInterface sat_act_gen;

    Cnf top_cnf;
    top_cnf.Convert2Cnf(top_cir, NULL, NULL, false); // do not optimize as it is very time consuming, and useless!

    int h_low_i_index = max_high_active_outputs*(1000-15*min_high_act_level)/1000;
    int h_high_i_index = max_high_active_outputs*(1000-5*min_high_act_level)/1000;
    int l_low_i_index = max_low_active_outputs*(1000-15*min_low_act_level)/1000;
    int l_high_i_index = max_low_active_outputs*(1000-5*min_low_act_level)/1000;
    int min_low_margin = max_low_active_outputs/64;
    int min_high_margin = max_high_active_outputs/64;
    if (l_low_i_index < min_low_margin)
      l_low_i_index = min_low_margin;
    if (l_high_i_index > max_low_active_outputs-min_low_margin)
      l_high_i_index = max_low_active_outputs-min_low_margin;
    if (h_low_i_index < min_high_margin)
      h_low_i_index = min_high_margin;
    if (h_high_i_index > max_high_active_outputs-min_high_margin)
      h_high_i_index = max_high_active_outputs-min_high_margin;

    int cur_i = 0;
    if (do_high_sa) {
      while (cur_i < 10) {
        LiteralType lit_i = top_cnf.literal_str_map["yh0" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(-lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
      while (cur_i < h_low_i_index) {
        LiteralType lit_i = top_cnf.literal_str_map["yh" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(-lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
      cur_i = h_high_i_index;
      while (cur_i < max_high_active_outputs) {
        LiteralType lit_i = top_cnf.literal_str_map["yh" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
    }
    if (do_low_sa) {
      cur_i = 0;
      while (cur_i < 10) {
        LiteralType lit_i = top_cnf.literal_str_map["yl0" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(-lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
      while (cur_i < l_low_i_index) {
        LiteralType lit_i = top_cnf.literal_str_map["yl" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(-lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
      cur_i = l_high_i_index;
      while (cur_i < max_low_active_outputs) {
        LiteralType lit_i = top_cnf.literal_str_map["yl" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
    }

    // add ff constraints
    for (int ffn = 0; ffn < new_spec_cir.ffs.size(); ffn++) {
      Node* ffnode = new_spec_cir.ffs[ffn];
      LiteralType lit_i = top_cnf.GetLiteral(ffnode->name);
      LiteralVector clause_i;
      clause_i.push_back(allsig_vals[ffnode->index]?lit_i:-lit_i);
      top_cnf.clauses.push_back(clause_i);
    }


    sat_act_gen.AddCnf(top_cnf.clauses, 0);

    for (int sol_num = 0; sol_num < num_sat_solutions; sol_num++) {
      if (do_high_sa)
        cout << "High activity level range is: " << 100-((100*h_high_i_index)/max_high_active_outputs) << "~" << 100-((100*h_low_i_index)/max_high_active_outputs) << " %" << endl;
      if (do_low_sa)
        cout << "Low activity level range is: " << 100-((100*l_high_i_index)/max_low_active_outputs) << "~" << 100-((100*l_low_i_index)/max_low_active_outputs) << " %" << endl;

      //top_cnf.WriteDimacs("run"+std::to_string(counter)+"_"+std::to_string(sol_num)+".cnf");

      bool sat_res = sat_act_gen.GetAnswer();
      if (!sat_res) {
        top_cnf.clauses.clear();
        top_cnf.clauses.shrink_to_fit();
        top_cnf.literal_str_map.clear();
        top_cnf.literal_node_map.clear();
        if (sol_num == 0) {
          MSG("No Answer!");
          delete top_cir;
          delete diff_cir;
          return -1;
        }
        MSG("No more answers!");
        break;
      }
      MSG("Found the patterns! dumping to pattern files...");

      LiteralVector clause_block_pi;
      LiteralVector clause_block_ff;
      clause_block_pi.reserve(new_spec_cir.inputs.size());
      clause_block_ff.reserve(new_spec_cir.ffs.size());

      ValVector pi_vect(new_spec_cir.inputs.size(), false);
      ValVector po_vect(new_spec_cir.outputs.size(), false);
      ValVector ff_vect(new_spec_cir.ffs.size(), false);

      ofstream patfile("patterns.csv", ofstream::out|ofstream::app);
      if (patfile.good()) {
        patfile << counter << "," << sol_num << ",";
        for (int pi_k = 0; pi_k < new_spec_cir.inputs.size(); pi_k++) {
          LiteralType lit_block = top_cnf.GetLiteral(new_spec_cir.inputs[pi_k]->name);
          bool sat_res_var1 = sat_act_gen.ModelNum(lit_block) == l__True;
          patfile << sat_res_var1 ? '1' : '0';
          if (lit_block <= 0)
            patfile << '?';
          pi_vect[pi_k] = sat_res_var1;
          clause_block_pi.push_back(sat_res_var1 ? -lit_block : lit_block);
        }
        patfile << ",";
        for (int ff_k = 0; ff_k < new_spec_cir.ffs.size(); ff_k++) {
          LiteralType lit_block = top_cnf.GetLiteral(new_spec_cir.ffs[ff_k]->name);
          bool sat_res_var1 = sat_act_gen.ModelNum(lit_block) == l__True;
          patfile << sat_res_var1 ? '1' : '0';
          ff_vect[ff_k] = sat_res_var1;
          if (lit_block <= 0)
            patfile << '?';
          else if (sat_res_var1 != allsig_vals[new_spec_cir.ffs[ff_k]->index])
            patfile << 'X';
        }
        //patfile << endl;
        //patfile.close();
      }
      else {
        ERR("Error opening the pattern file to write the results!");
        return -1;
      }

      cur_patterns.resize(sol_num+1);
      cur_patterns[sol_num].pis = pi_vect;
      cur_patterns[sol_num].ff_cur = ff_vect;
      allsig_vect[sol_num].resize(new_spec_cir.all_nodes.size(), false);
      new_spec_cir.Simulate(pi_vect, po_vect, ff_vect, allsig_vect[sol_num]);
      cur_patterns[sol_num].ff_nxt = ff_vect;

      patfile << ",";
      for (int ff_k = 0; ff_k < new_spec_cir.ffs.size(); ff_k++) {
        patfile << ff_vect[ff_k] ? '1' : '0';
        LiteralType lit_block = top_cnf.GetLiteral(new_spec_cir.ffs[ff_k]->inputs[0]->name);
        bool sat_res_var1 = sat_act_gen.ModelNum(lit_block) == l__True;
        if (lit_block <= 0)
          patfile << '?';
        else if (sat_res_var1 != ff_vect[ff_k])
          patfile << 'X';
        clause_block_ff.push_back(sat_res_var1 ? -lit_block : lit_block);
      }
      //patfile << endl;
      //patfile.close();

      int num_diff_nodes = 0;
      bit64 diff_count = 0;
      for (int i = 0; i < allsig_vect[sol_num].size(); i++) {
        Node *cur_node = new_spec_cir.all_nodes[i];
        if (!cur_node->is_input && cur_node->type != NODE_ZERO && cur_node->type != NODE_ONE) {
          num_diff_nodes++;
          if (allsig_vect[sol_num][i] != allsig_vals[i])
            diff_count++;
        }
      }
      vector<bit64> diff_count_all(all_regions.size(), 0);
      for (int k = 0; k < all_regions.size(); k++) {
        for (int i = 0; i < all_regions[k].size(); i++) {
          int index = all_regions[k][i]->index;
          if (allsig_vect[sol_num][index] != allsig_vals[index])
            diff_count_all[k]++;
        }
      }
      bit64 diff_count_h = 0;
      for (int i = 0; i < sa_high_nodes.size(); i++) {
        int index = sa_high_nodes[i]->index;
        if (allsig_vect[sol_num][index] != allsig_vals[index])
          diff_count_h++;
      }
      bit64 diff_count_l = 0;
      for (int i = 0; i < sa_low_nodes.size(); i++) {
        int index = sa_low_nodes[i]->index;
        if (allsig_vect[sol_num][index] != allsig_vals[index])
          diff_count_l++;
      }
      cout << "Activity rate is: " << setprecision(2) << (double) 100*diff_count / num_diff_nodes << " (" << num_diff_nodes << ") :::";
      if (do_high_sa) {
        cur_patterns[sol_num].high_sa = 100.0 * diff_count_h / sa_high_nodes.size();
        cout << " High=" << (double) 100 * diff_count_h / sa_high_nodes.size();
        patfile << "," << (double) 100 * diff_count_h / sa_high_nodes.size();
      }
      if (do_low_sa) {
        cur_patterns[sol_num].low_sa = 100.0 * diff_count_l / sa_low_nodes.size();
        cout << " Low=" << (double) 100 * diff_count_l / sa_low_nodes.size();
        patfile << "," << (double) 100 * diff_count_l / sa_low_nodes.size();
      }
      cout << endl;
      cur_patterns[sol_num].all_sas.resize(all_regions.size());
      for (int i = 0; i < all_regions.size(); i++) {
        cur_patterns[sol_num].all_sas[i] = 100.0 * diff_count_all[i] / all_regions[i].size();
        cout << "\t" << i << "=" << (double) 100 * diff_count_all[i] / all_regions[i].size();
        patfile << "," << (double) 100 * diff_count_all[i] / all_regions[i].size();
      }
      cout << endl;

      patfile << endl;
      patfile.close();

      if (do_high_sa) {
        ClauseVector blocks_clauses;
        blocks_clauses.push_back(clause_block_pi);
        blocks_clauses.push_back(clause_block_ff);
        int cur_i = h_high_i_index;
        bool sat_res_var = false;
        do {
          LiteralType lit_i = top_cnf.literal_str_map["yh" + std::to_string(cur_i)];
          sat_res_var = sat_act_gen.ModelNum(lit_i) == l__True;
          LiteralVector clause_i;
          clause_i.push_back(lit_i);
          blocks_clauses.push_back(clause_i);
          cur_i--;
        } while (sat_res_var && cur_i > h_low_i_index);
        //sat_act_gen.AddCnf(blocks_clauses, top_cnf.clauses.size());
        sat_act_gen.AddCnf(blocks_clauses, 0);
        top_cnf.clauses.insert(top_cnf.clauses.end(), blocks_clauses.begin(), blocks_clauses.end());
        h_high_i_index = cur_i;
      }
      if (do_low_sa) {
        ClauseVector blocks_clauses;
        if (!do_high_sa) {
          blocks_clauses.push_back(clause_block_pi);
          blocks_clauses.push_back(clause_block_ff);
        }
        int cur_i = l_low_i_index;
        bool sat_res_var = false;
        do {
          LiteralType lit_i = top_cnf.literal_str_map["yl" + std::to_string(cur_i)];
          sat_res_var = sat_act_gen.ModelNum(lit_i) == l__True;
          LiteralVector clause_i;
          clause_i.push_back(-lit_i);
          blocks_clauses.push_back(clause_i);
          cur_i++;
        } while (!sat_res_var && cur_i > l_low_i_index);
        //sat_act_gen.AddCnf(blocks_clauses, top_cnf.clauses.size());
        sat_act_gen.AddCnf(blocks_clauses, 0);
        top_cnf.clauses.insert(top_cnf.clauses.end(), blocks_clauses.begin(), blocks_clauses.end());
        l_low_i_index = cur_i;
      }

    }

    auto end_time = chrono::system_clock::now();
    chrono::duration<double> diff_time = end_time - start_time;
    cout << endl << "Pattern creation by SAT time is: " << setprecision(3) << diff_time.count() << " s" << endl;

    top_cnf.clauses.clear();
    top_cnf.clauses.shrink_to_fit();
    top_cnf.literal_str_map.clear();
    top_cnf.literal_node_map.clear();
  }
  else {
    ERR("Cannot write top level module!");
  }

  delete top_cir;
  delete diff_cir;

  if (counter > 0) {
    int min_ret = 0;
    int index = num_sat_solutions-1;
    while (allsig_vect[index].size() == 0)
      index--;
    while (index >= 0) {
      int ret = SatSimIncActivity(new_spec_cir, all_regions, sa_high_nodes, sa_low_nodes, allsig_vect[index], ref_high_act_level, ref_low_act_level, min_high_act_level, min_low_act_level, max_high_active_outputs, max_low_active_outputs, patterns, counter - 1);
      if (ret >= 0) {
        for (int i = 0; i < patterns.size(); i++) {
          patterns[i].push_back(cur_patterns[index]);
        }
        return ret + 1;
      }
      else if (ret < min_ret)
        min_ret = ret;
      //allsig_vect[index].clear();
      //allsig_vect[index].shrink_to_fit();
      index--;
    }
    // try again with higher SA!
    MSG("trying again with higher/lower SA!");
    index = num_sat_solutions-1;
    while (allsig_vect[index].size() == 0)
      index--;
    while (index >= 0) {
      int ret = SatSimIncActivity(new_spec_cir, all_regions, sa_high_nodes, sa_low_nodes, allsig_vect[index], ref_high_act_level, ref_low_act_level, min_high_act_level*115/100, min_low_act_level*85/100, max_high_active_outputs, max_low_active_outputs, patterns, counter - 1);
      if (ret >= 0) {
        for (int i = 0; i < patterns.size(); i++) {
          patterns[i].push_back(cur_patterns[index]);
        }
        return ret + 1;
      }
      else if (ret < min_ret)
        min_ret = ret;
      //allsig_vect[index].clear();
      //allsig_vect[index].shrink_to_fit();
      index--;
    }
    // try again with lower SA!
    MSG("trying again with lower/higher SA!");
    index = num_sat_solutions-1;
    while (allsig_vect[index].size() == 0)
      index--;
    while (index >= 0) {
      int ret = SatSimIncActivity(new_spec_cir, all_regions, sa_high_nodes, sa_low_nodes, allsig_vect[index], ref_high_act_level, ref_low_act_level, min_high_act_level*85/100, min_low_act_level*115/100, max_high_active_outputs, max_low_active_outputs, patterns, counter - 1);
      if (ret >= 0) {
        for (int i = 0; i < patterns.size(); i++) {
          patterns[i].push_back(cur_patterns[index]);
        }
        return ret + 1;
      }
      else if (ret < min_ret)
        min_ret = ret;
      //allsig_vect[index].clear();
      //allsig_vect[index].shrink_to_fit();
      index--;
    }
    MSG("trying again with higher+/lower- SA!");
    index = num_sat_solutions-1;
    while (allsig_vect[index].size() == 0)
      index--;
    while (index >= 0) {
      int ret = SatSimIncActivity(new_spec_cir, all_regions, sa_high_nodes, sa_low_nodes, allsig_vect[index], ref_high_act_level, ref_low_act_level, min_high_act_level*13/10, min_low_act_level*75/100, max_high_active_outputs, max_low_active_outputs, patterns, counter - 1);
      if (ret >= 0) {
        for (int i = 0; i < patterns.size(); i++) {
          patterns[i].push_back(cur_patterns[index]);
        }
        return ret + 1;
      }
      else if (ret < min_ret)
        min_ret = ret;
      //allsig_vect[index].clear();
      //allsig_vect[index].shrink_to_fit();
      index--;
    }
    // try again with lower SA!
    MSG("trying again with lower-/higher+ SA!");
    index = num_sat_solutions-1;
    while (allsig_vect[index].size() == 0)
      index--;
    while (index >= 0) {
      int ret = SatSimIncActivity(new_spec_cir, all_regions, sa_high_nodes, sa_low_nodes, allsig_vect[index], ref_high_act_level, ref_low_act_level, min_high_act_level*75/100, min_low_act_level*13/10, max_high_active_outputs, max_low_active_outputs, patterns, counter - 1);
      if (ret >= 0) {
        for (int i = 0; i < patterns.size(); i++) {
          patterns[i].push_back(cur_patterns[index]);
        }
        return ret + 1;
      }
      else if (ret < min_ret)
        min_ret = ret;
      //allsig_vect[index].clear();
      //allsig_vect[index].shrink_to_fit();
      index--;
    }
    return min_ret-1;
  }
  if (cur_patterns.size() > 0) {
    patterns.resize(cur_patterns.size());
    for (int i = 0; i < cur_patterns.size(); i++) {
      patterns[i].push_back(cur_patterns[i]);
    }
  }
  return 0;
}


int SatSimIncActivity(Circuit &new_spec_cir, vector<NodeVector> &all_regions, NodeVector &sa_high_nodes, NodeVector &sa_low_nodes, int ref_high_act_level, int ref_low_act_level, int min_high_act_level, int min_low_act_level, int max_high_active_outputs, int max_low_active_outputs, PatternVector2D &patterns, int counter = 0) {
  const int num_sat_solutions = 10;
  ValVector allsig_vect[num_sat_solutions];
  ValVector allsig_vect_init[num_sat_solutions];
  PatternVector init_patterns;
  PatternVector cur_patterns;
  PatternVector2D nxt_patterns;
  patterns.clear();

  bool do_high_sa = (max_high_active_outputs > 0);
  bool do_low_sa = (max_low_active_outputs > 0);

  if (!do_high_sa && !do_low_sa) {
    return -1;
  }
  if (do_high_sa && min_high_act_level <= min_low_act_level) {
    ERR("high activity should be greater than low activity!");
    return -1;
  }

  cout << "preparing files...." << endl;

  int num_high_nodes = sa_high_nodes.size();
  if (num_high_nodes > 4096) {
    sa_high_nodes.resize(4096);
    MSG("Cutting the high SA nodes to: 4096");
  }
  else if (num_high_nodes == 4096) {
  }
  else if (num_high_nodes > 1024) {
    sa_high_nodes.resize(1024);
    MSG("Cutting the high SA nodes to: 1024");
  }
  /*else if (num_high_nodes%128) {
    // scale the high SA
    min_high_act_level = (min_high_act_level*num_high_nodes)/(num_high_nodes+128-num_high_nodes%128);
    MSG("scaled high SA to: "+std::to_string(min_high_act_level));
  }*/
  int num_low_nodes = sa_low_nodes.size();
  if (num_low_nodes > 4096) {
    sa_low_nodes.resize(4096);
    MSG("Cutting the low SA nodes to: 4096");
  }
  else if (num_low_nodes == 4096) {
  }
  else if (num_low_nodes > 1024) {
    sa_low_nodes.resize(1024);
    MSG("Cutting the low SA nodes to: 1024");
  }
  /*else if (num_low_nodes%128) {
    // scale the low SA
    min_low_act_level = (min_low_act_level*num_low_nodes)/(num_low_nodes+128-num_low_nodes%128);
    MSG("scaled low SA to: "+std::to_string(min_low_act_level));
  }*/

  // create a new circuit
  Circuit *diff_cir = GetDiffCircuit(new_spec_cir, sa_high_nodes, sa_low_nodes, false);
  Circuit *top_cir = CreateTopSatCircuit(diff_cir, sa_high_nodes.size(), sa_low_nodes.size(), max_high_active_outputs, max_low_active_outputs, min_high_act_level, min_low_act_level);

  if (top_cir) {
    cout << "preparing for SAT..." <<endl;

    auto start_time = chrono::system_clock::now();

    GInterface sat_act_gen;

    Cnf top_cnf;
    top_cnf.Convert2Cnf(top_cir, NULL, NULL, false); // do not optimize as it is very time consuming, and useless!

    int h_low_i_index = max_high_active_outputs*(1000-15*min_high_act_level)/1000;
    int h_high_i_index = max_high_active_outputs*(1000-5*min_high_act_level)/1000;
    int l_low_i_index = max_low_active_outputs*(1000-15*min_low_act_level)/1000;
    int l_high_i_index = max_low_active_outputs*(1000-5*min_low_act_level)/1000;
    int min_low_margin = max_low_active_outputs/64;
    int min_high_margin = max_high_active_outputs/64;
    if (l_low_i_index < min_low_margin)
      l_low_i_index = min_low_margin;
    if (l_high_i_index > max_low_active_outputs-min_low_margin)
      l_high_i_index = max_low_active_outputs-min_low_margin;
    if (h_low_i_index < min_high_margin)
      h_low_i_index = min_high_margin;
    if (h_high_i_index > max_high_active_outputs-min_high_margin)
      h_high_i_index = max_high_active_outputs-min_high_margin;

    int cur_i = 0;
    if (do_high_sa) {
      while (cur_i < 10) {
        LiteralType lit_i = top_cnf.literal_str_map["yh0" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(-lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
      while (cur_i < h_low_i_index) {
        LiteralType lit_i = top_cnf.literal_str_map["yh" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(-lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
      cur_i = h_high_i_index;
      while (cur_i < max_high_active_outputs) {
        LiteralType lit_i = top_cnf.literal_str_map["yh" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
    }
    if (do_low_sa) {
      cur_i = 0;
      while (cur_i < 10) {
        LiteralType lit_i = top_cnf.literal_str_map["yl0" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(-lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
      while (cur_i < l_low_i_index) {
        LiteralType lit_i = top_cnf.literal_str_map["yl" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(-lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
      cur_i = l_high_i_index;
      while (cur_i < max_low_active_outputs) {
        LiteralType lit_i = top_cnf.literal_str_map["yl" + std::to_string(cur_i)];
        LiteralVector clause_i;
        clause_i.push_back(lit_i);
        top_cnf.clauses.push_back(clause_i);
        cur_i++;
      }
    }

    sat_act_gen.AddCnf(top_cnf.clauses, 0);

    for (int sol_num = 0; sol_num < num_sat_solutions; sol_num++) {
      if (do_high_sa)
        cout << "High activity level range is: " << 100-((100*h_high_i_index)/max_high_active_outputs) << "~" << 100-((100*h_low_i_index)/max_high_active_outputs) << " %" << endl;
      if (do_low_sa)
        cout << "Low activity level range is: " << 100-((100*l_high_i_index)/max_low_active_outputs) << "~" << 100-((100*l_low_i_index)/max_low_active_outputs) << " %" << endl;

      //top_cnf.WriteDimacs("run"+std::to_string(counter)+"_"+std::to_string(sol_num)+".cnf");

      bool sat_res = sat_act_gen.GetAnswer();
      if (!sat_res) {
        if (sol_num == 0) {
          MSG("No Answer!");
          delete top_cir;
          delete diff_cir;
          return -1;
        }
        MSG("No more answers!");
        for (int idx = 0; idx < num_sat_solutions; idx++) {
          allsig_vect_init[idx].clear();
          allsig_vect_init[idx].shrink_to_fit();
        }
        top_cnf.clauses.clear();
        top_cnf.clauses.shrink_to_fit();
        top_cnf.literal_str_map.clear();
        top_cnf.literal_node_map.clear();
        break;
      }
      MSG("Found the patterns! dumping to pattern files...");

      LiteralVector clause_block_pi;
      LiteralVector clause_block_ff;
      LiteralVector clause_block_pi_init;
      LiteralVector clause_block_ff_init;
      clause_block_pi.reserve(new_spec_cir.inputs.size());
      clause_block_ff.reserve(new_spec_cir.ffs.size());
      clause_block_pi_init.reserve(new_spec_cir.inputs.size());
      clause_block_ff_init.reserve(new_spec_cir.ffs.size());

      ValVector pi_vect(new_spec_cir.inputs.size(), false);
      ValVector po_vect(new_spec_cir.outputs.size(), false);
      ValVector ff_vect(new_spec_cir.ffs.size(), false);

      ofstream patfile("patterns.csv", ofstream::out|ofstream::app);

      if (patfile.good()) {
        patfile << counter+1 << "," << sol_num << ",";
        for (int pi_k = 0; pi_k < new_spec_cir.inputs.size(); pi_k++) {
          LiteralType lit_block = top_cnf.GetLiteral("cp_"+new_spec_cir.inputs[pi_k]->name);
          bool sat_res_var1 = sat_act_gen.ModelNum(lit_block) == l__True;
          patfile << sat_res_var1 ? '1' : '0';
          if (lit_block <= 0)
            patfile << '?';
          pi_vect[pi_k] = sat_res_var1;
          clause_block_pi_init.push_back(sat_res_var1 ? -lit_block : lit_block);
        }
        patfile << ",";
        for (int ff_k = 0; ff_k < new_spec_cir.ffs.size(); ff_k++) {
          LiteralType lit_block = top_cnf.GetLiteral("cp_"+new_spec_cir.ffs[ff_k]->name);
          bool sat_res_var1 = sat_act_gen.ModelNum(lit_block) == l__True;
          patfile << sat_res_var1 ? '1' : '0';
          ff_vect[ff_k] = sat_res_var1;
          if (lit_block <= 0)
            patfile << '?';
        }
        //patfile << endl;
        //patfile.close();
      }
      else {
        ERR("Error opening the pattern file to write the results!");
        return -1;
      }

      init_patterns.resize(sol_num+1);
      init_patterns[sol_num].pis = pi_vect;
      init_patterns[sol_num].ff_cur = ff_vect;
      allsig_vect_init[sol_num].resize(new_spec_cir.all_nodes.size(), false);
      new_spec_cir.Simulate(pi_vect, po_vect, ff_vect, allsig_vect_init[sol_num]);
      init_patterns[sol_num].ff_nxt = ff_vect;

      patfile << ",";
      for (int ff_k = 0; ff_k < new_spec_cir.ffs.size(); ff_k++) {
        patfile << ff_vect[ff_k] ? '1' : '0';
      }

      patfile << ","; // for high SA region
      for (int i = 0; i < all_regions.size(); i++)
        patfile << ",";
      patfile << endl;

      if (patfile.good()) {
        patfile << counter << "," << sol_num << ",";
        for (int pi_k = 0; pi_k < new_spec_cir.inputs.size(); pi_k++) {
          LiteralType lit_block = top_cnf.GetLiteral(new_spec_cir.inputs[pi_k]->name);
          bool sat_res_var1 = sat_act_gen.ModelNum(lit_block) == l__True;
          patfile << sat_res_var1 ? '1' : '0';
          if (lit_block <= 0)
            patfile << '?';
          pi_vect[pi_k] = sat_res_var1;
          clause_block_pi.push_back(sat_res_var1 ? -lit_block : lit_block);
        }
        patfile << ",";
        for (int ff_k = 0; ff_k < new_spec_cir.ffs.size(); ff_k++) {
          patfile << ff_vect[ff_k] ? '1' : '0';
          /*LiteralType lit_block = top_cnf.GetLiteral(new_spec_cir.ffs[ff_k]->name);
          bool sat_res_var1 = sat_act_gen.ModelNum(lit_block) == l__True;
          patfile << sat_res_var1 ? '1' : '0';
          ff_vect[ff_k] = sat_res_var1;
          if (lit_block <= 0)
            patfile << '?';
          else if (sat_res_var1 != allsig_vect_init[sol_num][new_spec_cir.ffs[ff_k]->index])
            patfile << 'X';
          clause_block_ff_init.push_back(sat_res_var1 ? -lit_block : lit_block);*/
        }
        //patfile << endl;
        //patfile.close();
      }
      else {
        ERR("Error opening the pattern file to write the results!");
        return -1;
      }

      cur_patterns.resize(sol_num+1);
      cur_patterns[sol_num].pis = pi_vect;
      cur_patterns[sol_num].ff_cur = ff_vect;
      allsig_vect[sol_num].resize(new_spec_cir.all_nodes.size(), false);
      new_spec_cir.Simulate(pi_vect, po_vect, ff_vect, allsig_vect[sol_num]);
      cur_patterns[sol_num].ff_nxt = ff_vect;

      patfile << ",";
      for (int ff_k = 0; ff_k < new_spec_cir.ffs.size(); ff_k++) {
        patfile << ff_vect[ff_k] ? '1' : '0';
        LiteralType lit_block = top_cnf.GetLiteral(new_spec_cir.ffs[ff_k]->inputs[0]->name);
        bool sat_res_var1 = sat_act_gen.ModelNum(lit_block) == l__True;
        if (lit_block <= 0)
          patfile << '?';
        else if (sat_res_var1 != ff_vect[ff_k])
          patfile << 'X';
        clause_block_ff.push_back(sat_res_var1 ? -lit_block : lit_block);
      }
      //patfile << endl;
      //patfile.close();

      int num_diff_nodes = 0;
      bit64 diff_count = 0;
      for (int i = 0; i < allsig_vect[sol_num].size(); i++) {
        Node *cur_node = new_spec_cir.all_nodes[i];
        if (!cur_node->is_input && cur_node->type != NODE_ZERO && cur_node->type != NODE_ONE && cur_node->type != NODE_DFF) {
          num_diff_nodes++;
          if (allsig_vect[sol_num][i] != allsig_vect_init[sol_num][i])
            diff_count++;
        }
      }
      vector<bit64> diff_count_all(all_regions.size(), 0);
      for (int k = 0; k < all_regions.size(); k++) {
        for (int i = 0; i < all_regions[k].size(); i++) {
          int index = all_regions[k][i]->index;
          if (allsig_vect[sol_num][index] != allsig_vect_init[sol_num][index])
            diff_count_all[k]++;
        }
      }
      bit64 diff_count_h = 0;
      for (int i = 0; i < sa_high_nodes.size(); i++) {
        int index = sa_high_nodes[i]->index;
        if (allsig_vect[sol_num][index] != allsig_vect_init[sol_num][index])
          diff_count_h++;
      }
      bit64 diff_count_l = 0;
      for (int i = 0; i < sa_low_nodes.size(); i++) {
        int index = sa_low_nodes[i]->index;
        if (allsig_vect[sol_num][index] != allsig_vect_init[sol_num][index])
          diff_count_l++;
      }
      cout << "Activity rate is: " << setprecision(2) << (double) 100*diff_count / num_diff_nodes << " (" << num_diff_nodes << ") :::";
      if (do_high_sa) {
        cur_patterns[sol_num].high_sa = 100.0 * diff_count_h / sa_high_nodes.size();
        cout << " High=" << (double) 100 * diff_count_h / sa_high_nodes.size();
        patfile << "," << (double) 100*diff_count_h / sa_high_nodes.size();;
      }
      if (do_low_sa) {
        cur_patterns[sol_num].low_sa = 100.0 * diff_count_l / sa_low_nodes.size();
        cout << " Low=" << (double) 100 * diff_count_l / sa_low_nodes.size();
        patfile << "," << (double) 100*diff_count_l / sa_low_nodes.size();;
      }
      cout << endl;
      cur_patterns[sol_num].all_sas.resize(all_regions.size());
      for (int i = 0; i < all_regions.size(); i++) {
        cur_patterns[sol_num].all_sas[i] = 100.0 * diff_count_all[i] / all_regions[i].size();
        cout << "\t" << i << "=" << (double) 100 * diff_count_all[i] / all_regions[i].size();
        patfile << "," << (double) 100 * diff_count_all[i] / all_regions[i].size();
      }
      cout << endl;

      patfile << endl;
      patfile.close();

      if (do_high_sa) {
        ClauseVector blocks_clauses;
        blocks_clauses.push_back(clause_block_pi);
        blocks_clauses.push_back(clause_block_pi_init);
        blocks_clauses.push_back(clause_block_ff);
        //blocks_clauses.push_back(clause_block_ff_init);
        int cur_i = h_high_i_index;
        bool sat_res_var = false;
        do {
          LiteralType lit_i = top_cnf.literal_str_map["yh" + std::to_string(cur_i)];
          sat_res_var = sat_act_gen.ModelNum(lit_i) == l__True;
          LiteralVector clause_i;
          clause_i.push_back(lit_i);
          blocks_clauses.push_back(clause_i);
          cur_i--;
        } while (sat_res_var && cur_i > h_low_i_index);
        //sat_act_gen.AddCnf(blocks_clauses, top_cnf.clauses.size());
        sat_act_gen.AddCnf(blocks_clauses, 0);
        top_cnf.clauses.insert(top_cnf.clauses.end(), blocks_clauses.begin(), blocks_clauses.end());
        h_high_i_index = cur_i;
      }
      if (do_low_sa) {
        ClauseVector blocks_clauses;
        if (!do_high_sa) {
          blocks_clauses.push_back(clause_block_pi);
          blocks_clauses.push_back(clause_block_pi_init);
          blocks_clauses.push_back(clause_block_ff);
          //blocks_clauses.push_back(clause_block_ff_init);
        }
        int cur_i = l_low_i_index;
        bool sat_res_var = false;
        do {
          LiteralType lit_i = top_cnf.literal_str_map["yl" + std::to_string(cur_i)];
          sat_res_var = sat_act_gen.ModelNum(lit_i) == l__True;
          LiteralVector clause_i;
          clause_i.push_back(-lit_i);
          blocks_clauses.push_back(clause_i);
          cur_i++;
        } while (!sat_res_var && cur_i > l_low_i_index);
        //sat_act_gen.AddCnf(blocks_clauses, top_cnf.clauses.size());
        sat_act_gen.AddCnf(blocks_clauses, 0);
        top_cnf.clauses.insert(top_cnf.clauses.end(), blocks_clauses.begin(), blocks_clauses.end());
        l_low_i_index = cur_i;
      }

    }

    auto end_time = chrono::system_clock::now();
    chrono::duration<double> diff_time = end_time - start_time;
    cout << endl << "Pattern creation by SAT time is: " << setprecision(3) << diff_time.count() << " s" << endl;

    for (int idx = 0; idx < num_sat_solutions; idx++) {
      allsig_vect_init[idx].clear();
      allsig_vect_init[idx].shrink_to_fit();
    }
    top_cnf.clauses.clear();
    top_cnf.clauses.shrink_to_fit();
    top_cnf.literal_str_map.clear();
    top_cnf.literal_node_map.clear();
  }
  else {
    ERR("Cannot write top level module!");
  }

  delete top_cir;
  delete diff_cir;

  if (counter > 0) {
    int min_ret = 0;
    int index = num_sat_solutions-1;
    while (allsig_vect[index].size() == 0)
      index--;
    while (index >= 0) {
      int ret = SatSimIncActivity(new_spec_cir, all_regions, sa_high_nodes, sa_low_nodes, allsig_vect[index], ref_high_act_level, ref_low_act_level, min_high_act_level, min_low_act_level, max_high_active_outputs, max_low_active_outputs, nxt_patterns, counter - 1);
      if (ret >= 0) {
        for (int i = 0; i < nxt_patterns.size(); i++) {
          nxt_patterns[i].push_back(cur_patterns[index]);
          nxt_patterns[i].push_back(init_patterns[index]);
        }
        patterns.insert(patterns.end(), nxt_patterns.begin(), nxt_patterns.end());
        //return ret + 1;
      }
      else if (ret < min_ret)
        min_ret = ret;
      //allsig_vect[index].clear();
      //allsig_vect[index].shrink_to_fit();
      index--;
    }
    if (patterns.size() > num_sat_solutions/2)
      return 0;
    // try again with higher SA!
    MSG("trying again with higher/lower SA!");
    index = num_sat_solutions-1;
    while (allsig_vect[index].size() == 0)
      index--;
    while (index >= 0) {
      int ret = SatSimIncActivity(new_spec_cir, all_regions, sa_high_nodes, sa_low_nodes, allsig_vect[index], ref_high_act_level, ref_low_act_level, min_high_act_level*115/100, min_low_act_level*85/100, max_high_active_outputs, max_low_active_outputs, nxt_patterns, counter - 1);
      if (ret >= 0) {
        for (int i = 0; i < nxt_patterns.size(); i++) {
          nxt_patterns[i].push_back(cur_patterns[index]);
          nxt_patterns[i].push_back(init_patterns[index]);
        }
        patterns.insert(patterns.end(), nxt_patterns.begin(), nxt_patterns.end());
        //return ret + 1;
      }
      else if (ret < min_ret)
        min_ret = ret;
      //allsig_vect[index].clear();
      //allsig_vect[index].shrink_to_fit();
      index--;
    }
    if (patterns.size() > num_sat_solutions/2)
      return 0;
    // try again with lower SA!
    MSG("trying again with lower/higher SA!");
    index = num_sat_solutions-1;
    while (allsig_vect[index].size() == 0)
      index--;
    while (index >= 0) {
      int ret = SatSimIncActivity(new_spec_cir, all_regions, sa_high_nodes, sa_low_nodes, allsig_vect[index], ref_high_act_level, ref_low_act_level, min_high_act_level*85/100, min_low_act_level*115/100, max_high_active_outputs, max_low_active_outputs, nxt_patterns, counter - 1);
      if (ret >= 0) {
        for (int i = 0; i < nxt_patterns.size(); i++) {
          nxt_patterns[i].push_back(cur_patterns[index]);
          nxt_patterns[i].push_back(init_patterns[index]);
        }
        patterns.insert(patterns.end(), nxt_patterns.begin(), nxt_patterns.end());
        //return ret + 1;
      }
      else if (ret < min_ret)
        min_ret = ret;
      //allsig_vect[index].clear();
      //allsig_vect[index].shrink_to_fit();
      index--;
    }
    if (patterns.size() > num_sat_solutions/2)
      return 0;
    MSG("trying again with higher+/lower- SA!");
    index = num_sat_solutions-1;
    while (allsig_vect[index].size() == 0)
      index--;
    while (index >= 0) {
      int ret = SatSimIncActivity(new_spec_cir, all_regions, sa_high_nodes, sa_low_nodes, allsig_vect[index], ref_high_act_level, ref_low_act_level, min_high_act_level*13/10, min_low_act_level*75/100, max_high_active_outputs, max_low_active_outputs, nxt_patterns, counter - 1);
      if (ret >= 0) {
        for (int i = 0; i < nxt_patterns.size(); i++) {
          nxt_patterns[i].push_back(cur_patterns[index]);
          nxt_patterns[i].push_back(init_patterns[index]);
        }
        patterns.insert(patterns.end(), nxt_patterns.begin(), nxt_patterns.end());
        //return ret + 1;
      }
      else if (ret < min_ret)
        min_ret = ret;
      //allsig_vect[index].clear();
      //allsig_vect[index].shrink_to_fit();
      index--;
    }
    if (patterns.size() > num_sat_solutions/2)
      return 0;
    // try again with lower SA!
    MSG("trying again with lower-/higher+ SA!");
    index = num_sat_solutions-1;
    while (allsig_vect[index].size() == 0)
      index--;
    while (index >= 0) {
      int ret = SatSimIncActivity(new_spec_cir, all_regions, sa_high_nodes, sa_low_nodes, allsig_vect[index], ref_high_act_level, ref_low_act_level, min_high_act_level*75/100, min_low_act_level*13/10, max_high_active_outputs, max_low_active_outputs, nxt_patterns, counter - 1);
      if (ret >= 0) {
        for (int i = 0; i < nxt_patterns.size(); i++) {
          nxt_patterns[i].push_back(cur_patterns[index]);
          nxt_patterns[i].push_back(init_patterns[index]);
        }
        patterns.insert(patterns.end(), nxt_patterns.begin(), nxt_patterns.end());
        //return ret + 1;
      }
      else if (ret < min_ret)
        min_ret = ret;
      //allsig_vect[index].clear();
      //allsig_vect[index].shrink_to_fit();
      index--;
    }
    if (patterns.size() > 0)
      return 0;
    return min_ret-1;
  }
  if (cur_patterns.size() > 0) {
    patterns.resize(cur_patterns.size());
    for (int i = 0; i < cur_patterns.size(); i++) {
      patterns[i].push_back(cur_patterns[i]);
      patterns[i].push_back(init_patterns[i]);
    }
  }
  return 0;
}

int WritePatterns(PatternVector2D patterns, string filename_prefix) {
  for (int i = 0; i < patterns.size(); i++) {
    ofstream pat_file(filename_prefix+"-"+std::to_string(i)+".csv");
    for (int j = patterns[i].size()-1; j >= 0; j--) {
      for (int k = 0; k < patterns[i][j].pis.size(); k++)
        pat_file << patterns[i][j].pis[k]?'1':'0';
      pat_file << ',';
      for (int k = 0; k < patterns[i][j].ff_cur.size(); k++)
        pat_file << patterns[i][j].ff_cur[k]?'1':'0';
      pat_file << ',';
      for (int k = 0; k < patterns[i][j].ff_nxt.size(); k++)
        pat_file << patterns[i][j].ff_nxt[k]?'1':'0';
      pat_file << ',';
      pat_file << patterns[i][j].low_sa;
      pat_file << ',';
      pat_file << patterns[i][j].high_sa;
      if (patterns[i][j].all_sas.size() > 0)
        for (int k = 0; k < patterns[i][j].all_sas.size(); k++) {
          pat_file << ',';
          pat_file << patterns[i][j].all_sas[k];
        }
      else
        for (int k = 0; k < patterns[i][j-1].all_sas.size(); k++) {
          pat_file << ',';
        }
      pat_file << endl;
    }
    pat_file.close();
  }

  return 0;
}

int ReadRegionSignals(string region_filename, Circuit& spec_cir, NodeVector& nodes) {
  ifstream reg_file(region_filename);
  if (!reg_file.good())
    return 1;

  char line[128];
  reg_file.getline(line, 127); // read the header
  reg_file.getline(line, 127); // read the next line
  while (reg_file.good()) {
    char* signame = line;
    while(*signame != ' ' && *signame != '\t') signame++;
    while(*signame == ' ' || *signame == '\t') signame++;
    Node* node = spec_cir.GetNode(signame);
    if (node) {
      if (node->type == NODE_DFF)
        MSG("ignoring FF node: "+node->name+" selected in region: "+region_filename);
      else
        nodes.push_back(node);
    }
    else
      ERR("Node not found: "+string(signame));
    reg_file.getline(line, 127); // read the next line
  }


  return 0;
}


int IncActivity(string spec_filename, string allregions_filename, string highregion_filename, int high_act_level, int num_cycles, bool rndsim) {
  if (high_act_level < 1 || high_act_level > 99) {
    ERR("Activity level should be between 1 and 99 (percent)");
    return -1;
  }

  Circuit spec_cir;

  spec_cir.ReadVerilogGL(spec_filename, false, false, false);

  for (int i = 0; i < spec_cir.outputs.size(); i++) {
    Node* ffnode = spec_cir.outputs[i];
    if (ffnode->type == NODE_DFF) {
      MSG("FF is output: " + ffnode->name);
      Node* newout = new Node(NODE_BUF);
      newout->is_output = true;
      newout->name = ffnode->name+"_newout";
      newout->inputs.push_back(ffnode);
      ffnode->outputs.push_back(newout);
      ffnode->is_output = false;
      spec_cir.outputs[i] = newout;
      spec_cir.all_nodes.push_back(newout);
      spec_cir.all_nodes_map[newout->name] = newout;
    }
  }

  for (int i = 0; i < spec_cir.ffs.size(); i++) {
    Node* ffnode = spec_cir.ffs[i];
    Node* ffinnode = ffnode->inputs[0];
    if (ffinnode->type == NODE_DFF) {
      MSG("FF input is FF: " + ffnode->name);
      Node* newin = new Node(NODE_BUF);
      newin->name = ffnode->name+"_newin";
      newin->inputs.push_back(ffinnode);
      newin->outputs.push_back(ffnode);
      ffnode->inputs[0] = newin;
      for (int j = 0; j < ffinnode->outputs.size(); j++) {
        if (ffinnode->outputs[j] == ffnode) {
          ffinnode->outputs[j] = newin;
          break;
        }
      }
      spec_cir.all_nodes.push_back(newin);
      spec_cir.all_nodes_map[newin->name] = newin;
    }
  }

  spec_cir.LevelizeSortTopological(false);
  spec_cir.SetIndexes();

  vector<NodeVector> all_regions;
  ifstream all_files(allregions_filename);
  char reg_filename[128];
  all_files.getline(reg_filename, 127);
  while (all_files.good()) {
    int num = 0;
    char* cc = reg_filename;
    while (*cc < '0' || *cc > '9')
      cc++;
    while (*cc >= '0' && *cc <= '9') {
      num *= 10;
      num += *cc-'0';
      cc++;
    }
    if (all_regions.size() < num+1)
      all_regions.resize(num+1);
    ReadRegionSignals(reg_filename, spec_cir, all_regions[num]);
    all_files.getline(reg_filename, 127);
  }

  NodeVector high_region;
  NodeVector low_region;
  ReadRegionSignals(highregion_filename, spec_cir, high_region);

  for (int i = 0; i < all_regions.size(); i++)
    cout << "region " << i << " has " << all_regions[i].size() << " signals." << endl;
  cout << "high region has " << high_region.size() << " signals." << endl;


  int num_high_nodes = high_region.size();
  int num_high_nodes_plus = num_high_nodes;
  if (num_high_nodes % 128) {
    num_high_nodes_plus += 128 - num_high_nodes % 128;
    // scale the high SA
    high_act_level = (high_act_level * num_high_nodes) / (num_high_nodes_plus);
    MSG("scaled high SA to: " + std::to_string(high_act_level));
  }

  PatternVector2D patterns;

  if (rndsim) {
    // Do random simulation for the initial state of the system
    ValVector pi_vect(spec_cir.inputs.size(), false);
    ValVector po_vect(spec_cir.outputs.size(), false);
    ValVector ff_vect(spec_cir.ffs.size(), false);
    ValVector allsig_vect(spec_cir.all_nodes.size(), false);

    for (int k = 0; k < pi_vect.size(); k++) {
      pi_vect[k] = (std::rand() % 2 == 1);
    }
    for (int k = 0; k < ff_vect.size(); k++) {
      ff_vect[k] = (std::rand() % 2 == 1);
    }

    ofstream patfile("patterns.csv", ofstream::out);
    if (patfile.good()) {
      patfile << num_cycles << ",0,";
      std::srand(std::time(0));
      for (int k = 0; k < pi_vect.size(); k++) {
        patfile << pi_vect[k];
      }
      patfile << ",";
      for (int k = 0; k < ff_vect.size(); k++) {
        patfile << ff_vect[k];
      }
      //patfile << endl;
      //patfile.close();
    }
    else {
      ERR("cannot open pattern file for writing the results!");
      return 0;
    }

    PatternType cur_pattern;
    cur_pattern.pis = pi_vect;
    cur_pattern.ff_cur = ff_vect;
    spec_cir.Simulate(pi_vect, po_vect, ff_vect, allsig_vect);
    cur_pattern.ff_nxt = ff_vect;

    patfile << ",";
    for (int k = 0; k < ff_vect.size(); k++) {
      patfile << ff_vect[k];
    }

    patfile << ","; // for high SA region
    for (int i = 0; i < all_regions.size(); i++)
      patfile << ",";

    patfile << endl;
    patfile.close();

    // TODO: low SA is considered 0, for now!
    int ret_val = SatSimIncActivity(spec_cir, all_regions, high_region, low_region, allsig_vect, high_act_level, 0, high_act_level, 0, num_high_nodes_plus, 0, patterns, num_cycles - 1);
    for (int i = 0; i < patterns.size(); i++)
      patterns[i].push_back(cur_pattern);
    cout << "Total number of solutions: " << patterns.size() << endl ;
    WritePatterns(patterns, "patterns"+std::to_string(num_cycles));
    return ret_val;
  }
  else {
    ofstream patfile("patterns.csv", ofstream::out);
    patfile << endl;
    patfile.close();
    int ret_val = SatSimIncActivity(spec_cir, all_regions, high_region, low_region, high_act_level, 0, high_act_level, 0, num_high_nodes_plus, 0, patterns, num_cycles - 1);
    cout << "Total number of solutions: " << patterns.size() << endl ;
    WritePatterns(patterns, "patterns"+std::to_string(num_cycles));
    return ret_val;
  }

  return 0;
}

int DecActivity(string spec_filename, string allregions_filename, string lowregion_filename, int low_act_level, int num_cycles, bool rndsim) {
  if (low_act_level < 1 || low_act_level > 99) {
    ERR("Activity level should be between 1 and 99 (percent)");
    return -1;
  }

  Circuit spec_cir;

  spec_cir.ReadVerilogGL(spec_filename, false, false, false);

  for (int i = 0; i < spec_cir.outputs.size(); i++) {
    Node* ffnode = spec_cir.outputs[i];
    if (ffnode->type == NODE_DFF) {
      MSG("FF is output: " + ffnode->name);
      Node* newout = new Node(NODE_BUF);
      newout->is_output = true;
      newout->name = ffnode->name+"_newout";
      newout->inputs.push_back(ffnode);
      ffnode->outputs.push_back(newout);
      ffnode->is_output = false;
      spec_cir.outputs[i] = newout;
      spec_cir.all_nodes.push_back(newout);
      spec_cir.all_nodes_map[newout->name] = newout;
    }
  }

  for (int i = 0; i < spec_cir.ffs.size(); i++) {
    Node* ffnode = spec_cir.ffs[i];
    Node* ffinnode = ffnode->inputs[0];
    if (ffinnode->type == NODE_DFF) {
      MSG("FF input is FF: " + ffnode->name);
      Node* newin = new Node(NODE_BUF);
      newin->name = ffnode->name+"_newin";
      newin->inputs.push_back(ffinnode);
      newin->outputs.push_back(ffnode);
      ffnode->inputs[0] = newin;
      for (int j = 0; j < ffinnode->outputs.size(); j++) {
        if (ffinnode->outputs[j] == ffnode) {
          ffinnode->outputs[j] = newin;
          break;
        }
      }
      spec_cir.all_nodes.push_back(newin);
      spec_cir.all_nodes_map[newin->name] = newin;
    }
  }

  spec_cir.LevelizeSortTopological(false);
  spec_cir.SetIndexes();

  vector<NodeVector> all_regions;
  ifstream all_files(allregions_filename);
  char reg_filename[128];
  all_files.getline(reg_filename, 127);
  while (all_files.good()) {
    int num = 0;
    char* cc = reg_filename;
    while (*cc < '0' || *cc > '9')
      cc++;
    while (*cc >= '0' && *cc <= '9') {
      num *= 10;
      num += *cc-'0';
      cc++;
    }
    if (all_regions.size() < num+1)
      all_regions.resize(num+1);
    ReadRegionSignals(reg_filename, spec_cir, all_regions[num]);
    all_files.getline(reg_filename, 127);
  }

  NodeVector high_region;
  NodeVector low_region;
  ReadRegionSignals(lowregion_filename, spec_cir, low_region);

  for (int i = 0; i < all_regions.size(); i++)
    cout << "region " << i << " has " << all_regions[i].size() << " signals." << endl;
  cout << "low region has " << low_region.size() << " signals." << endl;


  int num_low_nodes = low_region.size();
  int num_low_nodes_plus = num_low_nodes;
  if (num_low_nodes % 128) {
    num_low_nodes_plus += 128 - num_low_nodes % 128;
    // scale the high SA
    low_act_level = (low_act_level * num_low_nodes) / (num_low_nodes_plus);
    MSG("scaled low SA to: " + std::to_string(low_act_level));
  }

  PatternVector2D patterns;

  if (rndsim) {
    // Do random simulation for the initial state of the system
    ValVector pi_vect(spec_cir.inputs.size(), false);
    ValVector po_vect(spec_cir.outputs.size(), false);
    ValVector ff_vect(spec_cir.ffs.size(), false);
    ValVector allsig_vect(spec_cir.all_nodes.size(), false);

    ofstream patfile("patterns.csv", ofstream::out);
    if (patfile.good()) {
      patfile << num_cycles << ",0,";
      std::srand(std::time(0));
      for (int k = 0; k < pi_vect.size(); k++) {
        pi_vect[k] = (std::rand() % 2 == 1);
        patfile << pi_vect[k];
      }
      patfile << ",";
      for (int k = 0; k < ff_vect.size(); k++) {
        ff_vect[k] = (std::rand() % 2 == 1);
        patfile << ff_vect[k];
      }
      //patfile << endl;
      //patfile.close();
    }
    else {
      ERR("cannot open pattern file for writing the results!");
      return 0;
    }

    PatternType cur_pattern;
    cur_pattern.pis = pi_vect;
    cur_pattern.ff_cur = ff_vect;
    spec_cir.Simulate(pi_vect, po_vect, ff_vect, allsig_vect);
    cur_pattern.ff_nxt = ff_vect;

    patfile << ",";
    for (int k = 0; k < ff_vect.size(); k++) {
      patfile << ff_vect[k];
    }

    patfile << ","; // for high SA region
    for (int i = 0; i < all_regions.size(); i++)
      patfile << ",";

    patfile << endl;
    patfile.close();

    // TODO: low SA is considered 0, for now!
    int ret_val = SatSimIncActivity(spec_cir, all_regions, high_region, low_region, allsig_vect, 99, low_act_level, 99, low_act_level, 0, num_low_nodes_plus, patterns, num_cycles - 1);
    for (int i = 0; i < patterns.size(); i++)
      patterns[i].push_back(cur_pattern);
    cout << "Total number of solutions: " << patterns.size() << endl ;
    WritePatterns(patterns, "patterns"+std::to_string(num_cycles));
    return ret_val;
  }
  else {
    ofstream patfile("patterns.csv", ofstream::out);
    patfile << endl;
    patfile.close();
    int ret_val = SatSimIncActivity(spec_cir, all_regions, high_region, low_region, 99, low_act_level, 99, low_act_level, 0, num_low_nodes_plus, patterns, num_cycles - 1);
    cout << "Total number of solutions: " << patterns.size() << endl ;
    WritePatterns(patterns, "patterns"+std::to_string(num_cycles));
    return ret_val;
  }

  return 0;
}

int SandActivity(string spec_filename, string allregions_filename, string sandregion_filename, int low_act_level, int high_act_level, int num_cycles_low, int num_cycles_high) {
  if (high_act_level < 1 || high_act_level > 99) {
    ERR("Activity level should be between 1 and 99 (percent)");
    return -1;
  }
  if (low_act_level < 1 || low_act_level > 99) {
    ERR("Activity level should be between 1 and 99 (percent)");
    return -1;
  }
  if (high_act_level < low_act_level) {
    MSG("HighSA should be higher than LowSA!");
    return -1;
  }

  Circuit spec_cir;

  spec_cir.ReadVerilogGL(spec_filename, false, false, false);

  for (int i = 0; i < spec_cir.outputs.size(); i++) {
    Node* ffnode = spec_cir.outputs[i];
    if (ffnode->type == NODE_DFF) {
      MSG("FF is output: " + ffnode->name);
      Node* newout = new Node(NODE_BUF);
      newout->is_output = true;
      newout->name = ffnode->name+"_newout";
      newout->inputs.push_back(ffnode);
      ffnode->outputs.push_back(newout);
      ffnode->is_output = false;
      spec_cir.outputs[i] = newout;
      spec_cir.all_nodes.push_back(newout);
      spec_cir.all_nodes_map[newout->name] = newout;
    }
  }

  for (int i = 0; i < spec_cir.ffs.size(); i++) {
    Node* ffnode = spec_cir.ffs[i];
    Node* ffinnode = ffnode->inputs[0];
    if (ffinnode->type == NODE_DFF) {
      MSG("FF input is FF: " + ffnode->name);
      Node* newin = new Node(NODE_BUF);
      newin->name = ffnode->name+"_newin";
      newin->inputs.push_back(ffinnode);
      newin->outputs.push_back(ffnode);
      ffnode->inputs[0] = newin;
      for (int j = 0; j < ffinnode->outputs.size(); j++) {
        if (ffinnode->outputs[j] == ffnode) {
          ffinnode->outputs[j] = newin;
          break;
        }
      }
      spec_cir.all_nodes.push_back(newin);
      spec_cir.all_nodes_map[newin->name] = newin;
    }
  }

  spec_cir.LevelizeSortTopological(false);
  spec_cir.SetIndexes();

  vector<NodeVector> all_regions;
  ifstream all_files(allregions_filename);
  char reg_filename[128];
  all_files.getline(reg_filename, 127);
  while (all_files.good()) {
    int num = 0;
    char* cc = reg_filename;
    while (*cc < '0' || *cc > '9')
      cc++;
    while (*cc >= '0' && *cc <= '9') {
      num *= 10;
      num += *cc-'0';
      cc++;
    }
    if (all_regions.size() < num+1)
      all_regions.resize(num+1);
    ReadRegionSignals(reg_filename, spec_cir, all_regions[num]);
    all_files.getline(reg_filename, 127);
  }

  NodeVector sand_region;
  NodeVector temp_region;
  ReadRegionSignals(sandregion_filename, spec_cir, sand_region);

  for (int i = 0; i < all_regions.size(); i++)
    cout << "region " << i << " has " << all_regions[i].size() << " signals." << endl;
  cout << "sand region has " << sand_region.size() << " signals." << endl;


  int num_sand_nodes = sand_region.size();
  int num_sand_nodes_plus = num_sand_nodes;
  if (num_sand_nodes % 128) {
    num_sand_nodes_plus += 128 - num_sand_nodes % 128;
    // scale the high SA
    high_act_level = (high_act_level * num_sand_nodes) / (num_sand_nodes_plus);
    MSG("scaled high SA to: " + std::to_string(high_act_level));
    low_act_level = (low_act_level * num_sand_nodes) / (num_sand_nodes_plus);
    MSG("scaled low SA to: " + std::to_string(low_act_level));
  }

  PatternVector2D patternsLpre;
  PatternVector2D patternsHigh;
  PatternVector2D patternsLpost;

  ofstream patfile("patterns.csv", ofstream::out);
  patfile << endl;
  patfile.close();

  do {
    MSG("********** Low Region Pre **********");

    /*
    // Do random simulation for the initial state of the system
    ValVector pi_vect(spec_cir.inputs.size(), false);
    ValVector po_vect(spec_cir.outputs.size(), false);
    ValVector ff_vect(spec_cir.ffs.size(), false);
    ValVector allsig_vect(spec_cir.all_nodes.size(), false);

    for (int k = 0; k < pi_vect.size(); k++) {
      pi_vect[k] = (std::rand() % 2 == 1);
    }
    for (int k = 0; k < ff_vect.size(); k++) {
      ff_vect[k] = (std::rand() % 2 == 1);
    }

    ofstream patfile("patterns.csv", ofstream::out|ofstream::app);
    if (patfile.good()) {
      patfile << num_cycles_low*12/10 << ",0,";
      std::srand(std::time(0));
      for (int k = 0; k < pi_vect.size(); k++) {
        patfile << pi_vect[k];
      }
      patfile << ",";
      for (int k = 0; k < ff_vect.size(); k++) {
        patfile << ff_vect[k];
      }
      //patfile << endl;
      //patfile.close();
    }
    else {
      ERR("cannot open pattern file for writing the results!");
      return 0;
    }

    PatternType cur_pattern;
    cur_pattern.pis = pi_vect;
    cur_pattern.ff_cur = ff_vect;
    spec_cir.Simulate(pi_vect, po_vect, ff_vect, allsig_vect);
    cur_pattern.ff_nxt = ff_vect;

    patfile << ",";
    for (int k = 0; k < ff_vect.size(); k++) {
      patfile << ff_vect[k];
    }

    patfile << ","; // for high SA region
    for (int i = 0; i < all_regions.size(); i++)
      patfile << ",";

    patfile << endl;
    patfile.close();

    patternsLpre.clear();
    // TODO: low SA is considered 0, for now!
    int ret_val = SatSimIncActivity(spec_cir, all_regions, sand_region, temp_region, allsig_vect, low_act_level, 0, low_act_level, 0, num_sand_nodes_plus, 0, patternsLpre, num_cycles_low*12/10 - 1);
    for (int i = 0; i < patternsLpre.size(); i++)
      patternsLpost[i].push_back(cur_pattern);
    cout << "Total number of solutions: " << patternsLpre.size() << endl ;
    WritePatterns(patternsLpre, "patterns_lpre_"+std::to_string(num_cycles_low*12/10));
    */

    patternsLpre.clear();
    // TODO: low SA is considered 0, for now!
    int ret_val = SatSimIncActivity(spec_cir, all_regions, sand_region, temp_region, low_act_level, 0, low_act_level, 0, num_sand_nodes_plus, 0, patternsLpre, num_cycles_low*12/10 - 1);
    cout << "Total number of solutions: " << patternsLpre.size() << endl ;
    WritePatterns(patternsLpre, "patterns_lpre_"+std::to_string(num_cycles_low*12/10));

  } while (patternsLpre.size() <= 0);

  for (int pl1 = 0; pl1 < patternsLpre.size() && patternsHigh.size() <= 0; pl1 ++) {
    for (int pl2 = 0; pl2 < patternsLpre[pl1].size() && patternsHigh.size() <= 0; pl2++) {
      MSG("********** High Region sand **********");
      // Do random simulation for the initial state of the system
      ValVector pi_vect(spec_cir.inputs.size(), false);
      ValVector po_vect(spec_cir.outputs.size(), false);
      ValVector ff_vect(spec_cir.ffs.size(), false);
      ValVector allsig_vect(spec_cir.all_nodes.size(), false);

      pi_vect = patternsLpre[pl1][pl2].pis;
      ff_vect = patternsLpre[pl1][pl2].ff_cur;

      ofstream patfile("patterns.csv", ofstream::out|ofstream::app);
      if (patfile.good()) {
        patfile << num_cycles_low*12/10 << ",0,";
        std::srand(std::time(0));
        for (int k = 0; k < pi_vect.size(); k++) {
          patfile << pi_vect[k];
        }
        patfile << ",";
        for (int k = 0; k < ff_vect.size(); k++) {
          patfile << ff_vect[k];
        }
        //patfile << endl;
        //patfile.close();
      }
      else {
        ERR("cannot open pattern file for writing the results!");
        return 0;
      }

      PatternType cur_pattern;
      cur_pattern.pis = pi_vect;
      cur_pattern.ff_cur = ff_vect;
      spec_cir.Simulate(pi_vect, po_vect, ff_vect, allsig_vect);
      cur_pattern.ff_nxt = ff_vect;

      patfile << ",";
      for (int k = 0; k < ff_vect.size(); k++) {
        patfile << ff_vect[k];
      }

      patfile << ","; // for high SA region
      for (int i = 0; i < all_regions.size(); i++)
        patfile << ",";

      patfile << endl;
      patfile.close();

      patternsHigh.clear();
      // TODO: low SA is considered 0, for now!
      int ret_val = SatSimIncActivity(spec_cir, all_regions, sand_region, temp_region, allsig_vect, high_act_level, 0, high_act_level, 0, num_sand_nodes_plus, 0, patternsHigh, num_cycles_high - 1);
      for (int i = 0; i < patternsHigh.size(); i++)
        patternsHigh[i].push_back(cur_pattern);
      cout << "Total number of solutions: " << patternsHigh.size() << endl ;
      WritePatterns(patternsHigh, "patterns_high_"+std::to_string(num_cycles_high));
    }
  }

  for (int pl1 = 0; pl1 < patternsHigh.size() && patternsLpost.size() <= 0; pl1 ++) {
    for (int pl2 = 0; pl2 < patternsHigh[pl1].size() && patternsLpost.size() <= 0; pl2++) {
      MSG("********** Low Region Post **********");
      // Do random simulation for the initial state of the system
      ValVector pi_vect(spec_cir.inputs.size(), false);
      ValVector po_vect(spec_cir.outputs.size(), false);
      ValVector ff_vect(spec_cir.ffs.size(), false);
      ValVector allsig_vect(spec_cir.all_nodes.size(), false);

      pi_vect = patternsHigh[pl1][pl2].pis;
      ff_vect = patternsHigh[pl1][pl2].ff_cur;

      ofstream patfile("patterns.csv", ofstream::out|ofstream::app);
      if (patfile.good()) {
        patfile << num_cycles_low*12/10 << ",0,";
        std::srand(std::time(0));
        for (int k = 0; k < pi_vect.size(); k++) {
          patfile << pi_vect[k];
        }
        patfile << ",";
        for (int k = 0; k < ff_vect.size(); k++) {
          patfile << ff_vect[k];
        }
        //patfile << endl;
        //patfile.close();
      }
      else {
        ERR("cannot open pattern file for writing the results!");
        return 0;
      }

      PatternType cur_pattern;
      cur_pattern.pis = pi_vect;
      cur_pattern.ff_cur = ff_vect;
      spec_cir.Simulate(pi_vect, po_vect, ff_vect, allsig_vect);
      cur_pattern.ff_nxt = ff_vect;

      patfile << ",";
      for (int k = 0; k < ff_vect.size(); k++) {
        patfile << ff_vect[k];
      }

      patfile << ","; // for high SA region
      for (int i = 0; i < all_regions.size(); i++)
        patfile << ",";

      patfile << endl;
      patfile.close();

      patternsLpost.clear();
      // TODO: low SA is considered 0, for now!
      int ret_val = SatSimIncActivity(spec_cir, all_regions, sand_region, temp_region, allsig_vect, low_act_level, 0, low_act_level, 0, num_sand_nodes_plus, 0, patternsLpost, num_cycles_low - 1);
      for (int i = 0; i < patternsLpost.size(); i++)
        patternsLpost[i].push_back(cur_pattern);
      cout << "Total number of solutions: " << patternsLpost.size() << endl ;
      WritePatterns(patternsLpost, "patterns_lpost_"+std::to_string(num_cycles_low));
    }
  }


  return 0;
}