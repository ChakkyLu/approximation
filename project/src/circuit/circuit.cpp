#include "circuit.h"
#include <queue>

using namespace std;

#include "verilog_gate_parser.h"
#include "verilog_glib_parser.h"
#include "blif_parser.h"

namespace nodecircuit {

  int Circuit::ApplyInOutSimplify(ValVector &input_vals, ValVector &output_vals) {
    if (input_vals.size() > inputs.size() /*|| output_vals.size() != outputs.size()*/)
      return -1;

    NodeVector::iterator in_iter = inputs.begin();
    for (long i = 0; i < input_vals.size(); i++) {
      Node *inode = *in_iter;
      inode->is_input = false;
      inode->inputs.clear();
      inode->type = input_vals[i] ? NODE_ONE : NODE_ZERO;
      ++in_iter;
    }
    inputs.erase(inputs.begin(), in_iter);

    return Simplify();
  }

  int Circuit::Simplify() {
    NodeVector::iterator it;
    NodeSet orphan_nodes;
    // TODO: add simplification for the new 2-input and/nand gates
    for (long index = 0; index < all_nodes.size(); index++) {
      Node *cur_node = all_nodes[index];
      if (cur_node->type == NODE_ZERO) {
        for (long k = 0; k < cur_node->outputs.size(); k++) {
          Node *out_node = cur_node->outputs[k];
          switch (out_node->type) {
            case NODE_BUF:
            case NODE_AND:
              for (long m = 0; m < out_node->inputs.size(); m++) {
                Node *n = out_node->inputs[m];
                if (n == cur_node)
                  continue;
                it = n->outputs.begin();
                while (it != n->outputs.end()) {
                  if (*it == out_node) {
                    n->outputs.erase(it);
                    break;
                  }
                  ++it;
                }
                if (n->outputs.size() == 0 && !n->is_output && !n->is_input && n->type != NODE_DFF) {
                  orphan_nodes.insert(n);
                }
              }
              out_node->inputs.clear();
              out_node->type = NODE_ZERO;
              break;
            case NODE_NOT:
            case NODE_NAND:
              for (long m = 0; m < out_node->inputs.size(); m++) {
                Node *n = out_node->inputs[m];
                if (n == cur_node)
                  continue;
                it = n->outputs.begin();
                while (it != n->outputs.end()) {
                  if (*it == out_node) {
                    n->outputs.erase(it);
                    break;
                  }
                  ++it;
                }
                if (n->outputs.size() == 0 && !n->is_output && !n->is_input && n->type != NODE_DFF) {
                  orphan_nodes.insert(n);
                }
              }
              out_node->inputs.clear();
              out_node->type = NODE_ONE;
              break;
            case NODE_OR:
            case NODE_XOR: // TODO: only 2-input XOR
              it = out_node->inputs.begin();
              while (it != out_node->inputs.end()) {
                if (*it == cur_node) {
                  out_node->inputs.erase(it);
                  break;
                }
                ++it;
              }
              if (out_node->inputs.size() == 1)
                out_node->type = NODE_BUF;
              break;
            case NODE_NOR:
            case NODE_XNOR: // TODO: only 2-input XNOR
              it = out_node->inputs.begin();
              while (it != out_node->inputs.end()) {
                if (*it == cur_node) {
                  out_node->inputs.erase(it);
                  break;
                }
                ++it;
              }
              if (out_node->inputs.size() == 1)
                out_node->type = NODE_NOT;
              break;
            default:
              if ((out_node->type == NODE_AND2_NP && out_node->inputs[1] == cur_node) || (out_node->type == NODE_AND2_PN && out_node->inputs[0] == cur_node)) {
                for (long m = 0; m < out_node->inputs.size(); m++) {
                  Node *n = out_node->inputs[m];
                  if (n == cur_node)
                    continue;
                  it = n->outputs.begin();
                  while (it != n->outputs.end()) {
                    if (*it == out_node) {
                      n->outputs.erase(it);
                      break;
                    }
                    ++it;
                  }
                  if (n->outputs.size() == 0 && !n->is_output && !n->is_input && n->type != NODE_DFF) {
                    orphan_nodes.insert(n);
                  }
                }
                out_node->inputs.clear();
                out_node->type = NODE_ZERO;
              }
              else if ((out_node->type == NODE_AND2_NP && out_node->inputs[0] == cur_node) || (out_node->type == NODE_AND2_PN && out_node->inputs[1] == cur_node)) {
                it = out_node->inputs.begin();
                while (it != out_node->inputs.end()) {
                  if (*it == cur_node) {
                    out_node->inputs.erase(it);
                    break;
                  }
                  ++it;
                }
                if (out_node->inputs.size() == 1)
                  out_node->type = NODE_BUF;
              }
              else if ((out_node->type == NODE_NAND2_NP && out_node->inputs[1] == cur_node) || (out_node->type == NODE_NAND2_PN && out_node->inputs[0] == cur_node)) {
                for (long m = 0; m < out_node->inputs.size(); m++) {
                  Node *n = out_node->inputs[m];
                  if (n == cur_node)
                    continue;
                  it = n->outputs.begin();
                  while (it != n->outputs.end()) {
                    if (*it == out_node) {
                      n->outputs.erase(it);
                      break;
                    }
                    ++it;
                  }
                  if (n->outputs.size() == 0 && !n->is_output && !n->is_input && n->type != NODE_DFF) {
                    orphan_nodes.insert(n);
                  }
                }
                out_node->inputs.clear();
                out_node->type = NODE_ONE;
              }
              else if ((out_node->type == NODE_NAND2_NP && out_node->inputs[0] == cur_node) || (out_node->type == NODE_NAND2_PN && out_node->inputs[1] == cur_node)) {
                it = out_node->inputs.begin();
                while (it != out_node->inputs.end()) {
                  if (*it == cur_node) {
                    out_node->inputs.erase(it);
                    break;
                  }
                  ++it;
                }
                if (out_node->inputs.size() == 1)
                  out_node->type = NODE_NOT;
              }
          }
        }
      }
      else if (cur_node->type == NODE_ONE) {
        for (long k = 0; k < cur_node->outputs.size(); k++) {
          Node *out_node = cur_node->outputs[k];
          switch (out_node->type) {
            case NODE_BUF:
            case NODE_OR:
              for (long m = 0; m < out_node->inputs.size(); m++) {
                Node *n = out_node->inputs[m];
                if (n == cur_node)
                  continue;
                it = n->outputs.begin();
                while (it != n->outputs.end()) {
                  if (*it == out_node) {
                    n->outputs.erase(it);
                    break;
                  }
                  ++it;
                }
                if (n->outputs.size() == 0 && !n->is_output && !n->is_input && n->type != NODE_DFF) {
                  orphan_nodes.insert(n);
                }
              }
              out_node->inputs.clear();
              out_node->type = NODE_ONE;
              break;
            case NODE_NOT:
            case NODE_NOR:
              for (long m = 0; m < out_node->inputs.size(); m++) {
                Node *n = out_node->inputs[m];
                if (n == cur_node)
                  continue;
                it = n->outputs.begin();
                while (it != n->outputs.end()) {
                  if (*it == out_node) {
                    n->outputs.erase(it);
                    break;
                  }
                  ++it;
                }
                if (n->outputs.size() == 0 && !n->is_output && !n->is_input && n->type != NODE_DFF) {
                  orphan_nodes.insert(n);
                }
              }
              out_node->inputs.clear();
              out_node->type = NODE_ZERO;
              break;
            case NODE_AND:
            case NODE_XNOR: // TODO: only 2-input XNOR
              it = out_node->inputs.begin();
              while (it != out_node->inputs.end()) {
                if (*it == cur_node) {
                  out_node->inputs.erase(it);
                  break;
                }
                ++it;
              }
              if (out_node->inputs.size() == 1)
                out_node->type = NODE_BUF;
              break;
            case NODE_NAND:
            case NODE_XOR: // TODO: only 2-input XOR
              it = out_node->inputs.begin();
              while (it != out_node->inputs.end()) {
                if (*it == cur_node) {
                  out_node->inputs.erase(it);
                  break;
                }
                ++it;
              }
              if (out_node->inputs.size() == 1)
                out_node->type = NODE_NOT;
              break;
            default:
              if ((out_node->type == NODE_AND2_NP && out_node->inputs[0] == cur_node) || (out_node->type == NODE_AND2_PN && out_node->inputs[1] == cur_node)) {
                for (long m = 0; m < out_node->inputs.size(); m++) {
                  Node *n = out_node->inputs[m];
                  if (n == cur_node)
                    continue;
                  it = n->outputs.begin();
                  while (it != n->outputs.end()) {
                    if (*it == out_node) {
                      n->outputs.erase(it);
                      break;
                    }
                    ++it;
                  }
                  if (n->outputs.size() == 0 && !n->is_output && !n->is_input && n->type != NODE_DFF) {
                    orphan_nodes.insert(n);
                  }
                }
                out_node->inputs.clear();
                out_node->type = NODE_ZERO;
              }
              else if ((out_node->type == NODE_AND2_NP && out_node->inputs[1] == cur_node) || (out_node->type == NODE_AND2_PN && out_node->inputs[0] == cur_node)) {
                it = out_node->inputs.begin();
                while (it != out_node->inputs.end()) {
                  if (*it == cur_node) {
                    out_node->inputs.erase(it);
                    break;
                  }
                  ++it;
                }
                if (out_node->inputs.size() == 1)
                  out_node->type = NODE_NOT;
              }
              else if ((out_node->type == NODE_NAND2_NP && out_node->inputs[0] == cur_node) || (out_node->type == NODE_NAND2_PN && out_node->inputs[1] == cur_node)) {
                for (long m = 0; m < out_node->inputs.size(); m++) {
                  Node *n = out_node->inputs[m];
                  if (n == cur_node)
                    continue;
                  it = n->outputs.begin();
                  while (it != n->outputs.end()) {
                    if (*it == out_node) {
                      n->outputs.erase(it);
                      break;
                    }
                    ++it;
                  }
                  if (n->outputs.size() == 0 && !n->is_output && !n->is_input && n->type != NODE_DFF) {
                    orphan_nodes.insert(n);
                  }
                }
                out_node->inputs.clear();
                out_node->type = NODE_ONE;
              }
              else if ((out_node->type == NODE_NAND2_NP && out_node->inputs[1] == cur_node) || (out_node->type == NODE_NAND2_PN && out_node->inputs[0] == cur_node)) {
                it = out_node->inputs.begin();
                while (it != out_node->inputs.end()) {
                  if (*it == cur_node) {
                    out_node->inputs.erase(it);
                    break;
                  }
                  ++it;
                }
                if (out_node->inputs.size() == 1)
                  out_node->type = NODE_BUF;
              }
          }
        }
      }
      else if (cur_node->outputs.size() == 0 && !cur_node->is_output && !cur_node->is_input && cur_node->type != NODE_DFF) {
        orphan_nodes.insert(cur_node);
      }
    }

    while (orphan_nodes.size() > 0) {
      Node *onode = *orphan_nodes.begin();
      orphan_nodes.erase(orphan_nodes.begin());
      for (long m = 0; m < onode->inputs.size(); m++) {
        Node *n = onode->inputs[m];
        if (n->type == NODE_ZERO || n->type == NODE_ONE)
          continue;
        it = n->outputs.begin();
        while (it != n->outputs.end()) {
          if (*it == onode) {
            n->outputs.erase(it);
            break;
          }
          ++it;
        }
        if (n->outputs.size() == 0 && !n->is_output && !n->is_input && n->type != NODE_DFF) {
          orphan_nodes.insert(n);
        }
      }
    }

    int copy_index = 0, cur_index = 0;
    while (cur_index < all_nodes.size()) {
      Node *cur_node = all_nodes[cur_index];
      if (cur_node->is_input || cur_node->is_output) {
        if (copy_index < cur_index)
          all_nodes[copy_index] = all_nodes[cur_index];
        cur_index++;
        copy_index++;
      }
      else if (cur_node->type == NODE_ZERO || cur_node->type == NODE_ONE || cur_node->outputs.size() == 0) {
        all_nodes_map.erase(cur_node->name);
        cur_node->inputs.clear();
        cur_node->outputs.clear();
        delete cur_node;
        cur_index++;
      }
      else {
        if (copy_index < cur_index)
          all_nodes[copy_index] = all_nodes[cur_index];
        cur_index++;
        copy_index++;
      }
    }
    all_nodes.resize(copy_index);
    all_nodes.shrink_to_fit();

    return 0;
  }

  int Circuit::Simplify(NodeVector& target_nodes) {
    if (target_nodes.size() == 0)
      return -1;

    // find affected outputs
    NodeSet fanout_cone;
    NodeSet fanout_outputs;
    NodeSet process_list;
    for (int i = 0; i < target_nodes.size(); i++)
      process_list.insert(target_nodes[i]);
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

    if (fanout_outputs.size() == outputs.size())
      return 0;

    int org_inp_size = inputs.size();
    int org_outp_size = outputs.size();
    int org_node_size = all_nodes.size();

    // remove not affected outputs
    NodeVector::iterator iter = outputs.begin();
    while (iter != outputs.end()) {
      if (fanout_outputs.find(*iter) == fanout_outputs.end()) {
        (*iter)->is_output = false;
        iter = outputs.erase(iter);
      }
      else
        ++iter;
    }

    Simplify();

    //remove not affecting inputs
    iter = inputs.begin();
    while (iter != inputs.end()) {
      if ((*iter)->outputs.size() == 0) {
        (*iter)->is_input = false;
        iter = inputs.erase(iter);
      }
      else
        ++iter;
    }

    cout << "Nodes/PIs/POs after(before) simplification based on:";
    for (int i = 0; i < target_nodes.size(); i++)
      cout << " " << target_nodes[i]->name;
    cout << " -> " << all_nodes.size() << "(" << org_node_size << ")/" << inputs.size() << "(" << org_inp_size << ")/" << outputs.size() << "(" << org_outp_size << ")" << endl;

    return 0;
  }

  int Circuit::Simulate(ValVector &input_vals, ValVector &output_vals) {
    if (ffs.size() > 0)
      return 1;
    ValVector ff_vals;
    ValVector all_vals;
    all_vals.resize(all_nodes.size());
    return Simulate(input_vals, output_vals, ff_vals, all_vals);
  }

  int Circuit::Simulate(ValVector &input_vals, ValVector &output_vals, ValVector &ff_vals) {
    ValVector all_vals;
    all_vals.resize(all_nodes.size());
    return Simulate(input_vals, output_vals, ff_vals, all_vals);
  }

  int Circuit::Simulate(ValVector &input_vals, ValVector &output_vals, ValVector &ff_vals, ValVector &all_vals) {
    if (input_vals.size() != inputs.size()) {
      ERR("Simulation Error: values are not assigned for all the inputs!");
      return 1;
    }
    if (ff_vals.size() != ffs.size()) {
      ERR("Simulation Error: values are not assigned for all the FFs!");
      return 1;
    }
    output_vals.resize(outputs.size());
    all_vals.resize(all_nodes.size());
    for (long i = 0; i < inputs.size(); i++)
      all_vals[inputs[i]->index] = input_vals[i];
    for (long k = 0; k < ffs.size(); k++)
      all_vals[ffs[k]->index] = ff_vals[k];
    NodeVector::iterator node_it = all_nodes.begin();
    while ((*node_it)->inputs.size() == 0) {
      if ((*node_it)->type == NODE_ONE)
        all_vals[(*node_it)->index] = true;
      else if ((*node_it)->type == NODE_ZERO)
        all_vals[(*node_it)->index] = false;
      // TODO: check other cases should not happen
      ++node_it;
    }
    long i = 0;
    while (node_it != all_nodes.end()) {
      // do simulate!
      bool res;
      switch ((*node_it)->type) {
        case NODE_ZERO:
          all_vals[(*node_it)->index] = false;
          break;
        case NODE_ONE:
          all_vals[(*node_it)->index] = true;
          break;
        case NODE_BUF:
          all_vals[(*node_it)->index] = all_vals[(*node_it)->inputs[0]->index];
          break;
        case NODE_NOT:
          all_vals[(*node_it)->index] = !all_vals[(*node_it)->inputs[0]->index];
          break;
        case NODE_AND:
          res = all_vals[(*node_it)->inputs[0]->index];
          for (i = 1; res && i < (*node_it)->inputs.size(); i++)
            res &= all_vals[(*node_it)->inputs[i]->index];
          all_vals[(*node_it)->index] = res;
          break;
        case NODE_NAND:
          res = all_vals[(*node_it)->inputs[0]->index];
          for (i = 1; res && i < (*node_it)->inputs.size(); i++)
            res &= all_vals[(*node_it)->inputs[i]->index];
          all_vals[(*node_it)->index] = !res;
          break;
        case NODE_OR:
          res = all_vals[(*node_it)->inputs[0]->index];
          for (i = 1; !res && i < (*node_it)->inputs.size(); i++)
            res |= all_vals[(*node_it)->inputs[i]->index];
          all_vals[(*node_it)->index] = res;
          break;
        case NODE_NOR:
          res = all_vals[(*node_it)->inputs[0]->index];
          for (i = 1; !res && i < (*node_it)->inputs.size(); i++)
            res |= all_vals[(*node_it)->inputs[i]->index];
          all_vals[(*node_it)->index] = !res;
          break;
        case NODE_XOR:
          // TODO: currently only for 2 input gate
          res = all_vals[(*node_it)->inputs[0]->index] ^ all_vals[(*node_it)->inputs[1]->index];
          all_vals[(*node_it)->index] = res;
          break;
        case NODE_XNOR:
          // TODO: currently only for 2 input gate
          res = all_vals[(*node_it)->inputs[0]->index] ^ all_vals[(*node_it)->inputs[1]->index];
          all_vals[(*node_it)->index] = !res;
          break;
        case NODE_AND2_NP:
          res = (!all_vals[(*node_it)->inputs[0]->index]) & all_vals[(*node_it)->inputs[1]->index];
          all_vals[(*node_it)->index] = res;
          break;
        case NODE_NAND2_NP:
          res = (!all_vals[(*node_it)->inputs[0]->index]) & all_vals[(*node_it)->inputs[1]->index];
          all_vals[(*node_it)->index] = !res;
          break;
        case NODE_AND2_PN:
          res = all_vals[(*node_it)->inputs[0]->index] & (!all_vals[(*node_it)->inputs[1]->index]);
          all_vals[(*node_it)->index] = res;
          break;
        case NODE_NAND2_PN:
          res = all_vals[(*node_it)->inputs[0]->index] & (!all_vals[(*node_it)->inputs[1]->index]);
          all_vals[(*node_it)->index] = !res;
          break;
        case NODE_BLIF:
          ERR("general blif gate is not yet supported! " + (*node_it)->name);
          return 1;
          break;
        case NODE_DFF:
          // update at the last stage!
          break;
        default:
          ERR("unexpected gate in simulation! " + (*node_it)->name);
          return 1;
      }
      ++node_it;
    }

    for (long i = 0; i < outputs.size(); i++)
      output_vals[i] = all_vals[outputs[i]->index];
    for (long k = 0; k < ffs.size(); k++) {
      Node* ffnode = ffs[k];
      bool val = all_vals[ffnode->inputs[0]->index];
      all_vals[ffnode->index] = val;
      ff_vals[k] = val;
    }
    return 0;
  }

  int Circuit::Simulate(Val64Vector &input_vals, Val64Vector &output_vals) {
    Val64Vector all_vals;
    return Simulate(input_vals,output_vals, all_vals);
  }

  int Circuit::Simulate(Val64Vector &input_vals, Val64Vector &output_vals, Val64Vector &all_vals) {
    if (input_vals.size() != inputs.size()) {
      ERR("Simulation Error: values are not assigned for all the inputs!");
      return 1;
    }
    if (ffs.size() > 0) {
      ERR("Simulation Error: sequential circuit cannot be simulated this way!");
      return 1;
    }
    output_vals.resize(outputs.size());
    all_vals.resize(all_nodes.size());
    for (long i = 0; i < inputs.size(); i++)
      all_vals[inputs[i]->index] = input_vals[i];
    NodeVector::iterator node_it = all_nodes.begin();
    while (node_it != all_nodes.end() && (*node_it)->inputs.size() == 0) {
      if ((*node_it)->type == NODE_ONE)
        all_vals[(*node_it)->index] = ONE64;
      else if ((*node_it)->type == NODE_ZERO)
        all_vals[(*node_it)->index] = ZERO64;
      // TODO: check other cases should not happen
      ++node_it;
    }
    long i = 0;
    while (node_it != all_nodes.end()) {
      // do simulate!
      bit64 res;
      switch ((*node_it)->type) {
        case NODE_ZERO:
          all_vals[(*node_it)->index] = ZERO64;
          break;
        case NODE_ONE:
          all_vals[(*node_it)->index] = ONE64;
          break;
        case NODE_BUF:
          all_vals[(*node_it)->index] = all_vals[(*node_it)->inputs[0]->index];
          break;
        case NODE_NOT:
          all_vals[(*node_it)->index] = ~all_vals[(*node_it)->inputs[0]->index];
          break;
        case NODE_AND:
          res = all_vals[(*node_it)->inputs[0]->index];
          for (i = 1; i < (*node_it)->inputs.size(); i++)
            res &= all_vals[(*node_it)->inputs[i]->index];
          all_vals[(*node_it)->index] = res;
          break;
        case NODE_NAND:
          res = all_vals[(*node_it)->inputs[0]->index];
          for (i = 1; i < (*node_it)->inputs.size(); i++)
            res &= all_vals[(*node_it)->inputs[i]->index];
          all_vals[(*node_it)->index] = ~res;
          break;
        case NODE_OR:
          res = all_vals[(*node_it)->inputs[0]->index];
          for (i = 1; i < (*node_it)->inputs.size(); i++)
            res |= all_vals[(*node_it)->inputs[i]->index];
          all_vals[(*node_it)->index] = res;
          break;
        case NODE_NOR:
          res = all_vals[(*node_it)->inputs[0]->index];
          for (i = 1; i < (*node_it)->inputs.size(); i++)
            res |= all_vals[(*node_it)->inputs[i]->index];
          all_vals[(*node_it)->index] = ~res;
          break;
        case NODE_XOR:
          // TODO: currently only for 2 input gate
          res = all_vals[(*node_it)->inputs[0]->index] ^ all_vals[(*node_it)->inputs[1]->index];
          all_vals[(*node_it)->index] = res;
          break;
        case NODE_XNOR:
          // TODO: currently only for 2 input gate
          res = all_vals[(*node_it)->inputs[0]->index] ^ all_vals[(*node_it)->inputs[1]->index];
          all_vals[(*node_it)->index] = ~res;
          break;
        case NODE_AND2_NP:
          res = (~all_vals[(*node_it)->inputs[0]->index]) & all_vals[(*node_it)->inputs[1]->index];
          all_vals[(*node_it)->index] = res;
          break;
        case NODE_NAND2_NP:
          res = (~all_vals[(*node_it)->inputs[0]->index]) & all_vals[(*node_it)->inputs[1]->index];
          all_vals[(*node_it)->index] = ~res;
          break;
        case NODE_AND2_PN:
          res = all_vals[(*node_it)->inputs[0]->index] & (~all_vals[(*node_it)->inputs[1]->index]);
          all_vals[(*node_it)->index] = res;
          break;
        case NODE_NAND2_PN:
          res = all_vals[(*node_it)->inputs[0]->index] & (~all_vals[(*node_it)->inputs[1]->index]);
          all_vals[(*node_it)->index] = ~res;
          break;
        case NODE_BLIF:
          ERR("general blif gate is not yet supported! " + (*node_it)->name);
          return 1;
          break;
        case NODE_DFF:
          ERR("sequential circuit is not yet supported! " + (*node_it)->name);
          return 1;
          break;
        default:
          ERR("unexpected gate in simulation! " + (*node_it)->name);
          return 1;
      }
      ++node_it;
    }

    for (long i = 0; i < outputs.size(); i++)
      output_vals[i] = all_vals[outputs[i]->index];
    return 0;
  }


  int Circuit::LevelizeSortTopological(bool ignore_1input) {
    // sort from inputs to outputs
    Levelize(false);
    SortCurrentLevel();
    if (ignore_1input) {
      Levelize(true);
      SortCurrentLevel();
    }
    return 0;
  }

  int Circuit::SortCurrentLevel() {
    vector<int> level_count;
    vector<int> level_start_index;
    level_count.reserve(128);
    level_start_index.reserve(128);
    for (int i = 0; i < all_nodes.size(); i++) {
      int level = all_nodes[i]->level;
      if (level < 0)
        continue;
      if (level_count.size() < level+1) {
        level_count.resize(level + 1, 0);
        level_start_index.resize(level+1, 0);
      }
      int cur_count = level_count[level];
      if (cur_count == 0) {
        level_count[level] = 1;
        level_start_index[level] = i;
      }
      else {
        level_count[level] = cur_count+1;
      }
    }

    NodeVector cp_all_nodes(all_nodes);
    int index = 0;
    for (int level = 0; level < level_count.size(); level++) {
      int count = level_count[level];
      int cp_index = level_start_index[level];
      //cout << "level " << level << ": " << count << " starting from: " << cp_index << endl;
      for (int j = 0; j < count; ) {
        Node* node = cp_all_nodes[cp_index];
        if (node->level == level) {
          all_nodes[index] = node;
          index++;
          j++;
        }
        cp_index++;
      }
    }

    if (index < all_nodes.size()) {
      MSG("nodes with level=-1 removed! "+std::to_string(all_nodes.size())+ " < "+std::to_string(index));
      all_nodes.resize(index);
    }

    SetIndexes();

    return 0;

    // OLD sort! very slow
    NodeVector::iterator iter1, iter2, iter3;
    iter3 = all_nodes.begin();
    iter1 = iter3;
    ++iter3;
    iter2 = iter3;
    while (iter3 != all_nodes.end()) {
      while (iter1 != all_nodes.begin() && (*iter2)->level < (*iter1)->level) {
        Node *temp = *iter1;
        *iter1 = *iter2;
        *iter2 = temp;
        //MSG("swapping "+(*iter1)->name+" "+(*iter2)->name);
        iter1--;
        iter2--;
      }
      /*while (iter1 != all_nodes.begin() && (*iter2)->level == (*iter1)->level && (*iter2)->inputs.size() > (*iter1)->inputs.size()) {
        Node* temp = *iter1;
        *iter1 = *iter2;
        *iter2 = temp;
        MSG("swapping "+(*iter1)->name+" "+(*iter2)->name);
        iter1--;
        iter2--;
      }
      while (iter1 != all_nodes.begin() && (*iter2)->level == (*iter1)->level && (*iter2)->inputs.size() == (*iter1)->inputs.size() && (*iter2)->type < (*iter1)->type) {
        Node* temp = *iter1;
        *iter1 = *iter2;
        *iter2 = temp;
        MSG("swapping "+(*iter1)->name+" "+(*iter2)->name);
        iter1--;
        iter2--;
      }*/
      iter1 = iter3;
      ++iter3;
      iter2 = iter3;
    }

    SetIndexes();
    return 0;
  }

  int Circuit::Levelize(bool ignore_1input) {
    long i, j, k;

    queue<Node *> node_process_queue;
    for (i = 0; i < all_nodes.size(); i++) {
      Node *node = all_nodes[i];
      if (node->inputs.size() == 0 || node->type == NODE_DFF) {
        node->level = 0;
        for (j = 0; j < node->outputs.size(); j++)
          node_process_queue.push(node->outputs[j]);
      }
      else {
        node->level = -1; // reset the level
      }
    }
    while (node_process_queue.size() > 0) {
      Node *node = node_process_queue.front();
      node_process_queue.pop();
      if (node->level < 0) {
        int max_level = 0;
        for (k = 0; k < node->inputs.size() && max_level >= 0; k++) {
          int lev = node->inputs[k]->level;
          if (lev < 0)
            max_level = -1;
          else if (lev > max_level)
            max_level = lev;
        }
        if (max_level >= 0) {
          if (node->inputs.size() > 1)
            node->level = max_level + 1;
          else if (ignore_1input)
            node->level = max_level;
          else
            node->level = max_level + 1;
          for (k = 0; k < node->outputs.size(); k++) {
            Node *outnode = node->outputs[k];
            if (outnode->level < 0)
              node_process_queue.push(outnode);
          }
        }
      }
    }

    for (int k = 0; k < ffs.size(); k++) {
      Node* ffnode = ffs[k];
      ffnode->level = ffnode->inputs[0]->level + 1;
    }

    return 0;
  }

  int Circuit::SetIndexes() {
    for (long i = 0; i < all_nodes.size(); i++)
      all_nodes[i]->index = i;
    return 0;
  }

  int Circuit::ResetLevels() {
    for (long i = 0; i < all_nodes.size(); i++)
      all_nodes[i]->level = -1;
    return 0;
  }

  int Circuit::CreateCodedBlifValues() {
    for (long i = 0; i < all_nodes.size(); i++) {
      Node *node = all_nodes[i];
      if (node->type == NODE_BLIF)
        ((BlifNode *) node)->CreateCodedValues();
    }

    return 0;
  }

  int Circuit::ConvertEquivalentTypeBlif() {
    for (long i = 0; i < all_nodes.size(); i++) {
      Node *node = all_nodes[i];
      if (node->type == NODE_BLIF) {
        NodeType newtype = ((BlifNode *) node)->GetEquivalentType();
        if (newtype != node->type) {
          //MSG("node type changed: " + node->name);
          //cout << " -> " << node->type << " :: " << newtype << endl;
          node->type = newtype;
        }
      }
    }

    return 0;
  }

  int Circuit::RemoveBufNot() {
    long last_index = all_nodes.size() - 1;
    long cur_index = last_index;
    while (cur_index >= 0) {
      Node *cur_node = all_nodes[cur_index];
      if (cur_node->outputs.size() > 0 && cur_node->inputs.size() == 1 && !cur_node->is_output && !cur_node->is_input && cur_node->type != NODE_DFF) {
        Node *in_node = cur_node->inputs[0];
        if (in_node->outputs.size() == 1 && in_node->type < NODE_BLIF && in_node->type > NODE_UNKNOWN) {
          if (cur_node->type == NODE_NOT) {
            switch (in_node->type) {
              case NODE_ZERO:
                in_node->type = NODE_ONE;
                break;
              case NODE_ONE:
                in_node->type = NODE_ZERO;
                break;
              case NODE_BUF:
                in_node->type = NODE_NOT;
                break;
              case NODE_NOT:
                in_node->type = NODE_BUF;
                break;
              case NODE_AND:
                in_node->type = NODE_NAND;
                break;
              case NODE_NAND:
                in_node->type = NODE_AND;
                break;
              case NODE_OR:
                in_node->type = NODE_NOR;
                break;
              case NODE_NOR:
                in_node->type = NODE_OR;
                break;
              case NODE_XOR:
                in_node->type = NODE_XNOR;
                break;
              case NODE_XNOR:
                in_node->type = NODE_XOR;
                break;
              case NODE_AND2_NP:
                in_node->type = NODE_NAND2_NP;
                break;
              case NODE_NAND2_NP:
                in_node->type = NODE_AND2_NP;
                break;
              case NODE_AND2_PN:
                in_node->type = NODE_NAND2_PN;
                break;
              case NODE_NAND2_PN:
                in_node->type = NODE_AND2_PN;
                break;
            }
          }
          // else it is a buffer!
          in_node->outputs.clear();
          in_node->outputs = cur_node->outputs;
          for (long k = 0; k < cur_node->outputs.size(); k++) {
            Node *out_node = cur_node->outputs[k];
            for (int j = 0; j < out_node->inputs.size(); j++)
              if (out_node->inputs[j] == cur_node) {
                out_node->inputs[j] = in_node;
                break; // exit the loop
              }
          }
          all_nodes_map.erase(cur_node->name);
          //for (long i = cur_index; i < last_index; i++)
          //  all_nodes[i] = all_nodes[i + 1];
          //last_index--;
          delete cur_node;
          all_nodes[cur_index] = NULL; // delete later!
        }
      }
      cur_index--;
    }
    // delete removed nodes
    last_index = 0;
    while (last_index < all_nodes.size() && all_nodes[last_index] != NULL)
      last_index++;
    cur_index = last_index;
    while (cur_index < all_nodes.size() && all_nodes[cur_index] == NULL)
      cur_index++;
    if (last_index < all_nodes.size()) {
      while (cur_index < all_nodes.size()) {
        if (all_nodes[cur_index] != NULL) {
          all_nodes[last_index] = all_nodes[cur_index];
          last_index++;
        }
        cur_index++;
      }
      all_nodes.resize(last_index);
      all_nodes.shrink_to_fit();
    }

    return 0;
  }

  int Circuit::WriteVerilog(string filename) {
    ofstream out_file_stream(filename);
    if (!out_file_stream.is_open())
      return 0;

    long i, j, k;

    out_file_stream << "module " << name << " (" << endl;
    for (i = 0; i < outputs.size(); i++)
      out_file_stream << "  " << outputs[i]->name << "," << endl;
    for (i = 0; i < inputs.size() - 1; i++)
      out_file_stream << "  " << inputs[i]->name << "," << endl;
    out_file_stream << "  " << inputs[i]->name << ");" << endl;


    out_file_stream << "input" << endl;
    for (i = 0; i < inputs.size() - 1; i++)
      out_file_stream << "  " << inputs[i]->name << "," << endl;
    out_file_stream << "  " << inputs[i]->name << ";" << endl;

    out_file_stream << "output" << endl;
    for (i = 0; i < outputs.size() - 1; i++)
      out_file_stream << "  " << outputs[i]->name << "," << endl;
    out_file_stream << "  " << outputs[i]->name << ";" << endl;

    out_file_stream << "wire" << endl;
    bool first = true;
    for (i = 0; i < all_nodes.size(); i++) {
      Node *node = all_nodes[i];
      if (!node->is_input && !node->is_output) {
        if (!first)
          out_file_stream << "," << endl;
        out_file_stream << "  " << all_nodes[i]->name;
        first = false;
      }
    }
    out_file_stream << ";" << endl;

    for (i = 0; i < all_nodes.size(); i++) {
      Node *node = all_nodes[i];
      if (node->type == NODE_ONE) {
        out_file_stream << "  buf ( " << node->name << " , 1'b1 );" << endl;
        continue;
      }
      if (node->type == NODE_ZERO) {
        out_file_stream << "  buf ( " << node->name << " , 1'b0 );" << endl;
        continue;
      }
      if (node->is_input)
        continue;
      if (node->inputs.size() == 0)
        continue;
      switch (node->type) {
        case NODE_BUF:
          out_file_stream << "  buf (";
          break;
        case NODE_NOT:
          out_file_stream << "  not (";
          break;
        case NODE_AND:
          out_file_stream << "  and (";
          break;
        case NODE_NAND:
          out_file_stream << "  nand (";
          break;
        case NODE_OR:
          out_file_stream << "  or (";
          break;
        case NODE_NOR:
          out_file_stream << "  nor (";
          break;
        case NODE_XOR:
          out_file_stream << "  xor (";
          break;
        case NODE_XNOR:
          out_file_stream << "  xnor (";
          break;
        case NODE_BLIF:
          out_file_stream << "lut (";
          break;
        case NODE_DFF:
          out_file_stream << "dff (";
          break;
        default:
          out_file_stream << "GATE_" << node->type << " (";
          // TODO: check why this case has happened!
      }
      out_file_stream << " " << node->name;
      for (j = 0; j < node->inputs.size(); j++)
        out_file_stream << " , " << node->inputs[j]->name;
      out_file_stream << " );" << endl;
    }
    out_file_stream << "endmodule" << endl;

    out_file_stream.close();

    return 0;
  }

  int Circuit::WriteBlif(string filename) {
    ofstream out_file_stream(filename);
    if (!out_file_stream.is_open())
      return 0;

    long i, j, k;

    out_file_stream << ".model " << name << endl;
    out_file_stream << ".inputs ";
    for (i = 0; i < inputs.size(); i++)
      out_file_stream << " \\" << endl << " " << inputs[i]->name;
    out_file_stream << endl;
    out_file_stream << ".outputs ";
    for (i = 0; i < outputs.size(); i++)
      out_file_stream << " \\" << endl << " " << outputs[i]->name;
    out_file_stream << endl;
    for (i = 0; i < ffs.size(); i++)
      out_file_stream << ".latch " << ffs[i]->inputs[0]->name << " " << ffs[i]->name << endl;
    for (i = 0; i < all_nodes.size(); i++) {
      Node *node = all_nodes[i];
      if (node->type == NODE_ONE) {
        out_file_stream << ".names " << node->name << endl << "1" << endl;
        continue;
      }
      if (node->type == NODE_ZERO) {
        out_file_stream << ".names " << node->name << endl << "0" << endl;
        continue;
      }
      if (node->is_input)
        continue;
      if (node->type == NODE_DFF)
        continue;
      out_file_stream << ".names";
      for (j = 0; j < node->inputs.size(); j++)
        out_file_stream << " " << node->inputs[j]->name;
      out_file_stream << " " << node->name << endl;
      switch (node->type) {
        case NODE_BUF:
          out_file_stream << "1 1" << endl;
          break;
        case NODE_NOT:
          out_file_stream << "1 0" << endl;
          break;
        case NODE_AND:
          for (k = 0; k < node->inputs.size(); k++)
            out_file_stream << "1";
          out_file_stream << " 1" << endl;
          break;
        case NODE_NAND:
          for (k = 0; k < node->inputs.size(); k++)
            out_file_stream << "1";
          out_file_stream << " 0" << endl;
          break;
        case NODE_OR:
          for (k = 0; k < node->inputs.size(); k++)
            out_file_stream << "0";
          out_file_stream << " 0" << endl;
          break;
        case NODE_NOR:
          for (k = 0; k < node->inputs.size(); k++)
            out_file_stream << "0";
          out_file_stream << " 1" << endl;
          break;
        case NODE_XOR:
          switch (node->inputs.size()) {
            case 2:
              out_file_stream << "10 1" << endl;
              out_file_stream << "01 1" << endl;
              break;
            case 3:
              out_file_stream << "100 1" << endl;
              out_file_stream << "010 1" << endl;
              out_file_stream << "001 1" << endl;
              out_file_stream << "111 1" << endl;
              break;
              // TODO: for larger XORs
          }
          break;
        case NODE_XNOR:
          switch (node->inputs.size()) {
            case 2:
              out_file_stream << "10 0" << endl;
              out_file_stream << "01 0" << endl;
              break;
            case 3:
              out_file_stream << "100 0" << endl;
              out_file_stream << "010 0" << endl;
              out_file_stream << "001 0" << endl;
              out_file_stream << "111 0" << endl;
              break;
              // TODO: for larger XNORs
          }
          break;
        case NODE_AND2_NP:
          out_file_stream << "01 1" << endl;
          break;
        case NODE_NAND2_NP:
          out_file_stream << "01 0" << endl;
          break;
        case NODE_AND2_PN:
          out_file_stream << "10 1" << endl;
          break;
        case NODE_NAND2_PN:
          out_file_stream << "10 0" << endl;
          break;
        case NODE_BLIF:
          for (k = 0; k < ((BlifNode *) node)->str_values.size(); k++)
            out_file_stream << ((BlifNode *) node)->str_values[k] << endl;
          //for (set<uint64_t>::iterator iter = ((BlifNode*)node)->coded_values.begin(); iter != ((BlifNode*)node)->coded_values.end() ; ++iter)
          //  out_file_stream << hex << *iter << endl;
          break;
        default:
          for (k = 0; k < node->inputs.size(); k++)
            out_file_stream << "-";
          out_file_stream << " 0" << endl;
          out_file_stream << "# err " << node->type << endl;
          // this case should not happen
          // TODO: check why this case has happened
      }
    }
    out_file_stream << ".end" << endl;
    out_file_stream << endl;

    out_file_stream.close();
    return 0;
  }

  Circuit *Circuit::GetDuplicate(string input_prefix, string output_prefix, string internal_prefix) {
    // TODO: different prefix for ffs?
    Circuit *new_cir = new Circuit;
    new_cir->name = name;
    new_cir->inputs.reserve(inputs.size());
    new_cir->outputs.reserve(outputs.size());
    new_cir->ffs.reserve(ffs.size());
    new_cir->all_nodes.reserve(all_nodes.size());

    for (long i = 0; i < all_nodes.size(); i++) {
      Node *cur_node = all_nodes[i];
      Node *new_node = NULL;
      if (cur_node->type == NODE_BLIF)
        new_node = new BlifNode;
      else
        new_node = new Node(cur_node->type);
      new_node->is_input = cur_node->is_input;
      new_node->is_output = cur_node->is_output;
      if (cur_node->is_input)
        new_node->name = input_prefix + cur_node->name;
      else if (cur_node->is_output)
        new_node->name = output_prefix + cur_node->name;
      else
        new_node->name = internal_prefix + cur_node->name;
      new_node->inputs.reserve(cur_node->inputs.size());
      new_node->outputs.reserve(cur_node->outputs.size());
      new_node->level = cur_node->level;
      new_node->index = i;
      cur_node->index = i;
      if (cur_node->type == NODE_BLIF) {
        ((BlifNode *) new_node)->result_is_one = ((BlifNode *) cur_node)->result_is_one;
        ((BlifNode *) new_node)->str_values = ((BlifNode *) cur_node)->str_values;
        ((BlifNode *) new_node)->coded_values = ((BlifNode *) cur_node)->coded_values;
      }
      new_cir->all_nodes.push_back(new_node);
      new_cir->all_nodes_map[new_node->name] = new_node;
    }
    for (long i = 0; i < all_nodes.size(); i++) {
      Node *cur_node = all_nodes[i];
      Node *new_node = new_cir->all_nodes[i];
      for (long j = 0; j < cur_node->inputs.size(); j++)
        new_node->inputs.push_back(new_cir->all_nodes[cur_node->inputs[j]->index]);
      for (long k = 0; k < cur_node->outputs.size(); k++)
        new_node->outputs.push_back(new_cir->all_nodes[cur_node->outputs[k]->index]);
    }
    for (long j = 0; j < inputs.size(); j++)
      new_cir->inputs.push_back(new_cir->all_nodes[inputs[j]->index]);
    for (long k = 0; k < outputs.size(); k++)
      new_cir->outputs.push_back(new_cir->all_nodes[outputs[k]->index]);
    for (long m = 0; m < ffs.size(); m++)
      new_cir->ffs.push_back(new_cir->all_nodes[ffs[m]->index]);

    return new_cir;
  }

  int Circuit::ReadVerilogG(string filename, bool simplify, bool levelize_sort, bool ignore_1input) {
    Clear();

    FILE *infile = fopen(filename.c_str(), "r");
    if (!infile) {
      ERR("cannot open input file: " + filename);
      return 0;
    }
    ofstream parser_log("parser.log");
    veriloggateparser::VG_LEXER ilexer(infile, parser_log);
    ilexer.circuit = this;
    ilexer.vg_lex();
    fclose(infile);
    parser_log.close();

    system("rm parser.log");

    if (simplify) {
      LevelizeSortTopological(false);
      RemoveBufNot();
      Simplify();
      LevelizeSortTopological(false);
    }
    else if (levelize_sort) {
      LevelizeSortTopological(ignore_1input);
    }

    SetIndexes();

    return 1;
  }

  int Circuit::ReadVerilogGL(string filename, bool simplify, bool levelize_sort, bool ignore_1input) {
    Clear();

    FILE *infile = fopen(filename.c_str(), "r");
    if (!infile) {
      ERR("cannot open input file: " + filename);
      return 0;
    }
    ofstream parser_log("parser.log");
    verilogglibparser::VGL_LEXER ilexer(infile, parser_log);
    ilexer.circuit = this;
    ilexer.vgl_lex();
    fclose(infile);
    parser_log.close();

    system("rm parser.log");

    if (simplify) {
      LevelizeSortTopological(false);
      RemoveBufNot();
      Simplify();
      LevelizeSortTopological(false);
    }
    else if (levelize_sort) {
      LevelizeSortTopological(ignore_1input);
    }

    SetIndexes();

    return 1;
  }

  int Circuit::ReadBlif(string filename, bool simplify, bool levelize_sort, bool ignore_1input) {
    Clear();

    FILE *infile = fopen(filename.c_str(), "r");
    if (!infile) {
      ERR("cannot open input file: " + filename);
      return 0;
    }
    ofstream parser_log("parser.log");
    blifparser::BL_LEXER ilexer(infile, parser_log);
    ilexer.circuit = this;
    ilexer.bl_lex();
    fclose(infile);
    parser_log.close();
    system("rm parser.log");

    CreateCodedBlifValues();
    ConvertEquivalentTypeBlif();

    if (simplify) {
      LevelizeSortTopological(false);
      RemoveBufNot();
      Simplify();
      LevelizeSortTopological(false);
    }
    else if (levelize_sort) {
      LevelizeSortTopological(ignore_1input);
    }

    SetIndexes();

    return 1;
  }


}
