#include "global.h"
#include <chrono>

#include <iomanip>

using namespace std;

// ABC headers
#include "base/abc/abc.h"
#include "base/main/main.h"
#include "base/main/mainInt.h"

#include "circuit.h"

using namespace nodecircuit;

extern int LutDebug(Circuit *spec_cir, Circuit *impl_cir, ValVector& lut_result, vector<ValVector> &init_patterns);

const char default_spec_name[] = "new-spec-final.blif";
const char default_impl_name[] = "new-impl-final.blif";

const char default_org_spec_name[] = "spec.blif";
const char default_mod_impl_name[] = "impl.blif";

const char lut_sig_prefix[] = "lutv_";


int PrepareDebugFiles(string filename, vector<string> &targets, string out_node_name, int num_min_change) {

  int filenamesize = filename.size();

  Circuit new_spec_cir, new_impl_cir;

  // prepare spec and implementation circuits
  Circuit impl_cir;

  if (filename.substr(filenamesize - 2, 2) == ".v") {
    //impl_cir.ReadVerilogG(filename, false, false, false);
    ERR("currently only BLIF is supported!");
    return 0;
  }
  else if (filename.substr(filenamesize - 5, 5) == ".blif") {
    impl_cir.ReadBlif(filename, false, false, false);
  }
  else {
    ERR("unknown input file: " + filename);
    ERR("only gate-level Verilog and Blif are supported!");
    return 0;
  }

  NodeVector selected_nodes;
  for (int i = 0; i < targets.size(); i++) {
    Node *selectednode = impl_cir.GetNode(targets[i]);
    if (!selectednode) {
      ERR("selected node not found: " + targets[i]);
      return 0;
    }
    if (selectednode->type != NODE_BLIF) {
      // TODO: check that the original node has been Blif node!
      //ERR("currently only BLIF nodes are supported!");
      //return 0;
      selectednode->type = NODE_BLIF; // it is already coming from BLIF, so change it again to BLIF!
    }
    selected_nodes.push_back(selectednode);
  }

  // find affected outputs
  NodeSet fanout_cone;
  NodeSet fanout_outputs;
  NodeSet process_list;
  for (int i = 0; i < selected_nodes.size(); i++)
    process_list.insert(selected_nodes[i]);
  while (process_list.size() > 0) {
    Node *cur_node = *process_list.begin();
    process_list.erase(process_list.begin());
    fanout_cone.insert(cur_node);
    if (cur_node->is_output)
      fanout_outputs.insert(cur_node);
    for (int i = 0; i < cur_node->outputs.size(); i++) {
      Node *out_node = cur_node->outputs[i];
      if (fanout_cone.find(out_node) == fanout_cone.end())
        process_list.insert(out_node);
    }
  }

  Node* out_node_change = NULL;
  if (num_min_change > 0) {
    out_node_change = impl_cir.GetNode(out_node_name);
    if (!out_node_change) {
      ERR("output node for change not found! ignoring the change on: "+out_node_name);
      num_min_change = -1;
    }
    else if (fanout_outputs.find(out_node_change) == fanout_outputs.end()) {
      ERR("output node is not affected by the LUT nodes! ignoring the change on: "+out_node_name);
      num_min_change = -1;
    }
  }
  if (num_min_change > 0) {
    // add the xor gate with random minterm change!
    BlifNode* minterm_node = new BlifNode;
    minterm_node->name = out_node_change->name + "_minterm";
    NodeSet node_pis, node_fanins, process_nodes;
    process_nodes.insert(out_node_change);
    while (process_nodes.size() > 0) {
      Node* cur_node = *process_nodes.begin();
      process_nodes.erase(cur_node);
      if (cur_node->is_input)
        node_pis.insert(cur_node);
      else {
        node_fanins.insert(cur_node);
        for (int k = 0; k < cur_node->inputs.size(); k++) {
          Node* in_node = cur_node->inputs[k];
          if (node_fanins.find(in_node) == node_fanins.end())
            process_nodes.insert(in_node);
        }
      }
    }
    string minterm_str;
    minterm_str.reserve(node_pis.size());
    std::srand(std::time(0));
    for (NodeSet::iterator ns_it = node_pis.begin(); ns_it != node_pis.end(); ++ns_it) {
      minterm_node->inputs.push_back(*ns_it);
      //(*ns_it)->outputs.push_back(minterm_node); // do not add, as removal is difficult!
      minterm_str.push_back((char)((std::rand()%2)+'0'));
    }
    // TODO: if number of minterms > 1
    minterm_str += " 1";
    minterm_node->str_values.push_back(minterm_str);

    Node* out_node_xored = new Node(NODE_XOR);
    out_node_xored->name = out_node_change->name;
    out_node_xored->inputs.push_back(minterm_node);
    out_node_xored->inputs.push_back(out_node_change);
    minterm_node->outputs.push_back(out_node_xored);
    out_node_change->outputs.push_back(out_node_xored);

    out_node_change->name += "_org";
    // TODO: out_node_change outputs is assumed to be none!

    impl_cir.all_nodes.push_back(minterm_node);
    impl_cir.all_nodes.push_back(out_node_xored);

    impl_cir.all_nodes_map[minterm_node->name] = minterm_node;
    impl_cir.all_nodes_map[out_node_xored->name] = out_node_xored;
    impl_cir.all_nodes_map[out_node_change->name] = out_node_change;

    out_node_change->is_output = false;
    out_node_xored->is_output = true;

    for (int g = 0; g < impl_cir.outputs.size(); g++)
      if (impl_cir.outputs[g] == out_node_change) {
        impl_cir.outputs[g] = out_node_xored;
        break;
      }
    fanout_outputs.insert(out_node_xored);

    impl_cir.WriteBlif(default_org_spec_name);
  }
  else {
    string cmd_str = "cp " + filename + " " + default_org_spec_name;
    system(cmd_str.c_str());
  }
  string cmd_str = "cp " + filename + " " + default_mod_impl_name;
  system(cmd_str.c_str());


  // remove not affected outputs
  NodeVector::iterator iter = impl_cir.outputs.begin();
  while (iter != impl_cir.outputs.end()) {
    if (fanout_outputs.find(*iter) == fanout_outputs.end()) {
      (*iter)->is_output = false;
      iter = impl_cir.outputs.erase(iter);
    }
    else
      ++iter;
  }
  impl_cir.WriteBlif("temp-spec.blif");

  if (num_min_change > 0) {
    // revert the addition of xor gate with random minterm change!
    Node* out_node_xored = out_node_change->outputs[out_node_change->outputs.size()-1];
    Node* minterm_node = out_node_xored->inputs[0];
    impl_cir.all_nodes.resize(impl_cir.all_nodes.size()-2);

    impl_cir.all_nodes_map.erase(minterm_node->name);
    impl_cir.all_nodes_map.erase(out_node_change->name);
    out_node_change->name = out_node_xored->name;
    impl_cir.all_nodes_map[out_node_change->name] = out_node_change;

    out_node_change->is_output = true;
    out_node_xored->is_output = false;
    out_node_change->outputs.resize(out_node_change->outputs.size()-1);
    for (int g = 0; g < impl_cir.outputs.size(); g++)
      if (impl_cir.outputs[g] == out_node_xored) {
        impl_cir.outputs[g] = out_node_change;
        break;
      }

    delete out_node_xored;
    delete minterm_node;
  }

  // replace nodes with LUT
  long num_added_outputs = 0;
  for (int k = 0; k < selected_nodes.size(); k++) {
    Node *selectednode = selected_nodes[k];
    if (!selectednode->is_output) {
      selectednode->is_output = true;
      impl_cir.outputs.push_back(selectednode);
      num_added_outputs++;
    }
    long num_org_inputs = selectednode->inputs.size();

    long num_lut_vars = 1 << num_org_inputs;
    for (long i = 0; i < num_lut_vars; i++) {
      string var_name = lut_sig_prefix + selectednode->name + "_" + to_string(i);
      Node *var_node = new Node;
      var_node->name = var_name;
      var_node->is_input = true;
      var_node->outputs.push_back(selectednode);
      selectednode->inputs.push_back(var_node);
      impl_cir.inputs.push_back(var_node);
      impl_cir.all_nodes.push_back(var_node);
      impl_cir.all_nodes_map[var_name] = var_node;
    }
    BlifNode *selectednode_blif = (BlifNode *) selectednode;
    selectednode_blif->str_values.clear();
    string dc_sub_str(num_lut_vars, '-');
    for (long i = 0; i < num_lut_vars; i++) {
      string bin_str(num_org_inputs, '0');
      long temp_num = i;
      for (long j = 0; j < num_org_inputs; j++) {
        if (temp_num % 2)
          bin_str[j] = '1';
        temp_num /= 2;
      }
      string dc_str = dc_sub_str;
      dc_str[i] = '1';
      string final_str = bin_str + dc_str + " 1";
      selectednode_blif->str_values.push_back(final_str);
    }

    impl_cir.WriteBlif("temp-impl.blif");
    ABC::Abc_Frame_t *pAbc;
    pAbc = ABC::Abc_FrameGetGlobalFrame();
    string command = "read temp-impl.blif; strash; dc2; dc2; dc2; write new-impl.blif";
    if (ABC::Cmd_CommandExecute(pAbc, command.c_str())) {
      ERR(string("ABC Command execution error: ") + command);
      return 1;
    }
    command = "read temp-spec.blif; strash; dc2; dc2; dc2; write new-spec.blif";
    if (ABC::Cmd_CommandExecute(pAbc, command.c_str())) {
      ERR(string("ABC Command execution error: ") + command);
      return 1;
    }

  }

  new_spec_cir.ReadBlif("new-spec.blif", true);
  new_spec_cir.name = "spec";
  new_spec_cir.WriteBlif(default_spec_name);

  new_impl_cir.ReadBlif("new-impl.blif", false, false, false);

  //remove the added node at the output
  long last_index = new_impl_cir.outputs.size() - 1;
  for (int k = 0; k < num_added_outputs; k++) {
    Node *ext_out_node = new_impl_cir.outputs[last_index - k];
    ext_out_node->is_output = false;
  }
  new_impl_cir.outputs.resize(last_index + 1 - num_added_outputs);

  new_impl_cir.LevelizeSortTopological(false);
  new_impl_cir.RemoveBufNot();
  new_impl_cir.Simplify();
  new_impl_cir.name = "impl";
  new_impl_cir.WriteBlif(default_impl_name);

  //LutDebug(&new_spec_cir, &new_impl_cir);

  system("rm temp-*");
  system("rm new-impl.*");
  system("rm new-spec.*");

  return 0;
}

int ModifySpecImpl(string specfilename, string implfilename, Circuit* impl_lut_cur, ValVector& lut_result) {
  if (lut_result.size() == 0) {
    ERR("cannot create impl file: no lut results!");
    return -1;
  }

  cout << "creating impl file after modification...";

  Circuit spec_cir;
  spec_cir.ReadBlif(specfilename, false, false, false);

  string cur_sig_name;
  BlifNode* cur_sig = NULL;
  ValVector cur_lut_vars;
  cur_lut_vars.reserve(lut_result.size());
  long lut_sig_index = 0;
  while (lut_sig_index < impl_lut_cur->inputs.size()) {
    if (impl_lut_cur->inputs[lut_sig_index]->name.substr(0,5) == lut_sig_prefix)
      break;
    lut_sig_index++;
  }
  for (long lut_index = 0; lut_index < lut_result.size(); lut_index++, lut_sig_index++) {
    string lut_name = impl_lut_cur->inputs[lut_sig_index]->name;
    long last_index = lut_name.size()-1;
    while (lut_name[last_index] != '_')
      last_index--;
    string sig_name = lut_name.substr(5, last_index-5);
    if (sig_name != cur_sig_name) {
      if (cur_sig) {
        bool sig_out = true;
        string sig_out_str = " 1";
        long true_cnt = 0;
        for (long i = 0; i < cur_lut_vars.size(); i++)
          if (cur_lut_vars[i])
            true_cnt++;
        if (true_cnt > cur_lut_vars.size()/2) {
          sig_out = false;
          sig_out_str = " 0";
        }
        cur_sig->result_is_one = sig_out;
        cur_sig->str_values.clear();
        cur_sig->coded_values.clear();
        if (true_cnt == 0) {
          cur_sig->str_values.push_back("0");
          cur_sig->inputs.clear();
          cur_sig->type == NODE_ZERO;
        }
        else if (true_cnt == cur_lut_vars.size()) {
          cur_sig->str_values.push_back("1");
          cur_sig->inputs.clear();
          cur_sig->type == NODE_ONE;
        }
        else
          for (long i = 0; i < cur_lut_vars.size(); i++) {
            if (cur_lut_vars[i] == sig_out) {
              string one_row;
              one_row.reserve(cur_sig->inputs.size()+2);
              long j = i;
              for (long k = 0; k < cur_sig->inputs.size(); k++) {
                one_row.push_back((j%2)?'1':'0');
                j /= 2;
              }
              one_row += sig_out_str;
              cur_sig->str_values.push_back(one_row);
            }
          }
      }
      cur_sig_name = sig_name;
      cur_sig = (BlifNode*)spec_cir.GetNode(sig_name);
      cur_lut_vars.clear();
    }
    cur_lut_vars.push_back(lut_result[lut_index]);
    if (lut_index == lut_result.size()-1) {
      bool sig_out = true;
      string sig_out_str = " 1";
      long true_cnt = 0;
      for (long i = 0; i < cur_lut_vars.size(); i++)
        if (cur_lut_vars[i])
          true_cnt++;
      if (true_cnt > cur_lut_vars.size()/2) {
        sig_out = false;
        sig_out_str = " 0";
      }
      cur_sig->result_is_one = sig_out;
      cur_sig->str_values.clear();
      cur_sig->coded_values.clear();
      if (true_cnt == 0) {
        cur_sig->str_values.push_back("0");
        cur_sig->inputs.clear();
        cur_sig->type == NODE_ZERO;
      }
      else if (true_cnt == cur_lut_vars.size()) {
        cur_sig->str_values.push_back("1");
        cur_sig->inputs.clear();
        cur_sig->type == NODE_ONE;
      }
      else
        for (long i = 0; i < cur_lut_vars.size(); i++) {
          if (cur_lut_vars[i] == sig_out) {
            string one_row;
            one_row.reserve(cur_sig->inputs.size()+2);
            long j = i;
            for (long k = 0; k < cur_sig->inputs.size(); k++) {
              one_row.push_back((j%2)?'1':'0');
              j /= 2;
            }
            one_row += sig_out_str;
            cur_sig->str_values.push_back(one_row);
          }
        }
    }
  }

  spec_cir.WriteBlif(implfilename.c_str());

  cout << "done" << endl;
  return 0;
}

int DoDebug(string spec_filename, string impl_filename, int num_init_patterns, string pattern_filename) {
  Circuit new_spec_cir, new_impl_cir;

  new_spec_cir.ReadBlif(spec_filename);
  new_impl_cir.ReadBlif(impl_filename);

  // prepare initial input vectors
  std::srand(std::time(0));

  vector<ValVector> input_patterns;
  ValVector cur_input;
  cur_input.resize(new_spec_cir.inputs.size());

  set<string> patterns_str;
  string cur_pattern_str(new_spec_cir.inputs.size(), '0');
  vector<int> pattern_mask; // unused inputs!
  for (int kk = 0; kk < new_spec_cir.inputs.size(); kk++) {
    if (new_spec_cir.inputs[kk]->outputs.size() == 0)
      pattern_mask.push_back(kk);
  }

  if (num_init_patterns < 0 && pattern_filename.size() > 0) { // read from pattern file
    MSG("reading inputs from file: "+pattern_filename);
    ifstream pat_file(pattern_filename);
    char c;
    pat_file.read(&c,1);
    while (pat_file.good()) {
      int i3;
      for (i3 = 0; i3 < cur_input.size() && pat_file.good(); i3++) {
        cur_pattern_str[i3] = c;
        if (c == '0')
          cur_input[i3] = false;
        else if (c == '1')
          cur_input[i3] = true;
        else
          i3--;
        pat_file.read(&c,1);
      }
      if (i3 == cur_input.size()) {
        for (int kk = 0; kk < pattern_mask.size(); kk++)
          cur_pattern_str[pattern_mask[kk]] = '0';
        if (patterns_str.find(cur_pattern_str) == patterns_str.end()) {
          input_patterns.push_back(cur_input);
          patterns_str.insert(cur_pattern_str);
        }
        else {
          MSG("skipping duplicate (effective) input pattern: "+cur_pattern_str);
        }
      }
      else
        ERR("values for all inputs are not specified in the pattern file! "+pattern_filename);
      pat_file.read(&c,1);
    }
    num_init_patterns = input_patterns.size();
    if (num_init_patterns <= 0) {
      ERR("error reading patterns from the file: "+pattern_filename);
      return 0;
    }
  }
  else {
    if (num_init_patterns <= 0) {
      MSG("number of initial patterns should be greater than 0 -> setting to 1 to continue...");
      num_init_patterns = 1;
    }
    input_patterns.resize(num_init_patterns);

    int i3 = 0; // how many times existing patterns generated!
    for (int i1 = 0; i1 < num_init_patterns && i3 < num_init_patterns; i1++) {
      for (int i2 = 0; i2 < cur_input.size(); i2++) {
        cur_input[i2] = std::rand() % 2;
        cur_pattern_str[i2] = cur_input[i2]?'1':'0';
      }
      for (int kk = 0; kk < pattern_mask.size(); kk++)
        cur_pattern_str[pattern_mask[kk]] = '0';
      if (patterns_str.find(cur_pattern_str) == patterns_str.end()) {
        input_patterns[i1] = cur_input;
        patterns_str.insert(cur_pattern_str);
      }
      else {
        i1--;
        i3++;
        MSG("skipping duplicate (effective) input pattern!");
      }

    }
  }

  ValVector lut_result;
  auto start_time = chrono::system_clock::now();
  LutDebug(&new_spec_cir, &new_impl_cir, lut_result, input_patterns);
  auto end_time = chrono::system_clock::now();
  chrono::duration<double> diff_time = end_time - start_time;
  cout << endl << "Debug time is: " << setprecision(3) << diff_time.count() << " s" << endl;

  //ModifySpecImpl(default_org_spec_name, default_mod_impl_name, &new_impl_cir, lut_result);
  ModifySpecImpl(default_mod_impl_name, default_mod_impl_name, &new_impl_cir, lut_result);

  return 0;
}

int CreateNodeMiter(string impl_filename, string node_name) {
  Circuit new_impl_cir;

  new_impl_cir.ReadBlif(impl_filename, false, false, false);

  Node* selected_node = new_impl_cir.GetNode(node_name);
  if (!selected_node) {
    ERR("selected node not found: " + node_name);
    return 0;
  }

  // find affected outputs
  NodeSet fanout_cone;
  NodeSet fanout_outputs;
  NodeSet process_list;
  process_list.insert(selected_node);
  while (process_list.size() > 0) {
    Node *cur_node = *process_list.begin();
    process_list.erase(process_list.begin());
    fanout_cone.insert(cur_node);
    if (cur_node->is_output)
      fanout_outputs.insert(cur_node);
    for (int i = 0; i < cur_node->outputs.size(); i++) {
      Node *out_node = cur_node->outputs[i];
      if (fanout_cone.find(out_node) == fanout_cone.end())
        process_list.insert(out_node);
    }
  }
  // remove all outputs
  for (int i = 0; i < new_impl_cir.outputs.size(); i++)
    new_impl_cir.outputs[i]->is_output = false;
  new_impl_cir.outputs.clear();
  // add selected node inputs as output
  for (int i = 0; i < selected_node->inputs.size(); i++) {
    selected_node->inputs[i]->is_output = true;
    new_impl_cir.outputs.push_back(selected_node->inputs[i]);
  }

  // duplicate and rename fanout cone nodes
  new_impl_cir.all_nodes.reserve(new_impl_cir.all_nodes.size()+fanout_cone.size());
  for (NodeSet::iterator node_it = fanout_cone.begin(); node_it != fanout_cone.end(); ++node_it) {
    Node* node0 = *node_it;
    Node* node1;
    if (node0->type == NODE_BLIF)
      node1 = new BlifNode;
    else
      node1 = new Node(node0->type);
    node1->is_output = node0->is_output;
    node1->is_input = node0->is_input;
    node1->name = node0->name+"_1";
    node1->inputs.reserve(node0->inputs.size());
    node1->outputs.reserve(node0->outputs.size());
    new_impl_cir.all_nodes.push_back(node1);
    new_impl_cir.all_nodes_map[node1->name] = node1;
  }
  for (NodeSet::iterator node_it = fanout_cone.begin(); node_it != fanout_cone.end(); ++node_it) {
    Node* node0 = *node_it;
    Node* node1 = new_impl_cir.all_nodes_map[node0->name+"_1"];
    for (int i = 0; i < node0->inputs.size(); i++)
      if (fanout_cone.find(node0->inputs[i]) == fanout_cone.end()) {
        Node* inode = node0->inputs[i];
        node1->inputs.push_back(new_impl_cir.all_nodes_map[inode->name]);
        inode->outputs.push_back(node1);
      }
      else
        node1->inputs.push_back(new_impl_cir.all_nodes_map[node0->inputs[i]->name+"_1"]);
    for (int i = 0; i < node0->outputs.size(); i++)
      node1->outputs.push_back(new_impl_cir.all_nodes_map[node0->outputs[i]->name+"_1"]);
  }
  for (NodeSet::iterator node_it = fanout_cone.begin(); node_it != fanout_cone.end(); ++node_it) {
    Node *node0 = *node_it;
    new_impl_cir.all_nodes_map.erase(node0->name);
    node0->name = node0->name+"_0";
    new_impl_cir.all_nodes_map[node0->name] = node0;
  }

  // create miter output
  Node* miter_out_node = new Node(NODE_OR);
  miter_out_node->is_input = false;
  miter_out_node->is_output = true;
  miter_out_node->name = "sp_miter_out";

  new_impl_cir.all_nodes.reserve(new_impl_cir.all_nodes.size()+fanout_outputs.size()+1);
  for (NodeSet::iterator node_it = fanout_outputs.begin(); node_it != fanout_outputs.end(); ++node_it) {
    Node* xor_node = new Node(NODE_XOR);
    xor_node->is_input = false;
    xor_node->is_output = false;
    xor_node->name = (*node_it)->name+"_1_xor";
    string n1name = (*node_it)->name;
    n1name[n1name.size()-1] = '1';
    Node* node0 = *node_it;
    Node* node1 = new_impl_cir.all_nodes_map[n1name];
    xor_node->inputs.push_back(node0);
    xor_node->inputs.push_back(node1);
    xor_node->outputs.push_back(miter_out_node);

    node0->outputs.push_back(xor_node);
    node1->outputs.push_back(xor_node);
    miter_out_node->inputs.push_back(xor_node);

    new_impl_cir.all_nodes.push_back(xor_node);
    new_impl_cir.all_nodes_map[xor_node->name] = xor_node;
  }

  new_impl_cir.outputs.push_back(miter_out_node);
  new_impl_cir.all_nodes.push_back(miter_out_node);
  new_impl_cir.all_nodes_map[miter_out_node->name] = miter_out_node;

  Node* selected_node_0 = new_impl_cir.all_nodes_map[node_name+"_0"];
  selected_node_0->inputs.clear();
  selected_node_0->type = NODE_ZERO;
  Node* selected_node_1 = new_impl_cir.all_nodes_map[node_name+"_1"];
  selected_node_1->inputs.clear();
  selected_node_1->type = NODE_ONE;


  // alternative approach without using abc for optimization
/*
  new_impl_cir.LevelizeSortTopological(false);
  new_impl_cir.Simplify();
  new_impl_cir.RemoveBufNot();
  new_impl_cir.LevelizeSortTopological(false);
  new_impl_cir.WriteBlif(selected_node->name+"_1_miter.blif");
*/


  new_impl_cir.WriteBlif("temp.blif");
  ABC::Abc_Frame_t *pAbc;
  pAbc = ABC::Abc_FrameGetGlobalFrame();
  string command = "read temp.blif; strash; dc2; dc2; dc2; write temp.blif";
  if (ABC::Cmd_CommandExecute(pAbc, command.c_str())) {
    ERR(string("ABC Command execution error: ") + command);
    return 1;
  }

  Circuit opt_impl_cir;
  opt_impl_cir.ReadBlif("temp.blif");
  opt_impl_cir.WriteBlif(selected_node->name+"_1_miter.blif");

  system("rm temp*");

  return 0;
}