// simple verilog structural gate library parser, without error recover, and without comment detection

%{
#include <cstdio>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
using namespace std;
%}

%top{
#include "circuit.h"
using namespace nodecircuit;
#include <string>
%}

%class{
Node* current_node;
Node* temp_node;
int min_range;
int max_range;
int mode;
std::string name;
std::string dff_d;
std::string dff_q;
std::string dff_qn;
virtual int wrap() { return 1; }
public:
Circuit* circuit;
%}

%init{
circuit = NULL;
current_node = NULL;
min_range = max_range = 0;
mode = 255;
// 128:start 129:output 130:input 131:wire
// 140:gateinput 145:assign
// 150:dffreset 151:dffqn 152:dffq 153:dffd 154:dffclock
// 254:nothing 255:finish others:gate
%}

newline  \r?\n
space [ \t\r\n]

index \[[0-9][0-9]*\]
range \[[0-9][0-9]*:[0-9][0-9]*\]

idindex [a-zA-Z][A-Za-z0-9_]*{index}
id  [a-zA-Z][A-Za-z0-9_]*
dotidempty  "."{id}"()"
dotid  "."{id}
paridindex  "("{idindex}")"
parid  "("{id}")"
suff1  [2-9]_[A-Za-z0-9_]*
suff2  _[A-Za-z0-9_]*
zero 1'b0
one 1'b1
parzero "("{zero}")"
parone "("{one}")"


%option fast
%option namespace=verilogglibparser
%option lexer=VGL_LEXER
%option prefix=vgl_

%%

"endmodule"    mode = 255;
"module"       mode = 128;
"output"       {mode = 129; min_range = 0; max_range = 0;}
"input"        {mode = 130; min_range = 0; max_range = 0;}
"wire"         {mode = 131; min_range = 0; max_range = 0;}
"assign"       mode = 145;
"BUFF"{suff2}         mode = NODE_BUF;
"BUF"{suff2}          mode = NODE_BUF;
"INV"{suff2}          mode = NODE_NOT;
"NAND"{suff1}         mode = NODE_NAND;
"AND"{suff1}          mode = NODE_AND;
"XNOR2"{suff2}        mode = NODE_XNOR;
"XOR2"{suff2}         mode = NODE_XOR;
"NOR"{suff1}          mode = NODE_NOR;
"OR"{suff1}           mode = NODE_OR;
"DFFR"{suff2}         mode = NODE_DFF;


{range} { // mode is not changed!
          string range_str = text();
          int idx = 1;
          max_range = 0;
          while (range_str[idx] != ':') {
            max_range *= 10;
            max_range += range_str[idx]-'0';
            idx++;
          }
          idx++;
          min_range = 0;
          while (range_str[idx] != ']') {
            min_range *= 10;
            min_range += range_str[idx]-'0';
            idx++;
          }
          if (min_range > max_range) {
            int temp = min_range;
            min_range = max_range;
            max_range = temp;
          }
        }

{idindex}  if (mode == 145) { 
             current_node = circuit->GetNode(string(text()));
           }

{id}  switch (mode) {
        case 128:
          circuit->name = text();
          mode = 254;
          break;
        case 129:
          for (int i = min_range; i <= max_range; i++) {
            current_node = new Node;
            current_node->name = text();
            if (min_range != max_range)
              current_node->name += "["+std::to_string(i)+"]";
            current_node->is_output = true;
            circuit->outputs.push_back(current_node);
            circuit->all_nodes.push_back(current_node);
            circuit->all_nodes_map[current_node->name] = current_node;
          }
          break;
        case 130:
          for (int i = min_range; i <= max_range; i++) {
            current_node = new Node;
            current_node->name = text();
            if (min_range != max_range)
              current_node->name += "["+std::to_string(i)+"]";
            current_node->is_input = true;
            circuit->inputs.push_back(current_node);
            circuit->all_nodes.push_back(current_node);
            circuit->all_nodes_map[current_node->name] = current_node;
          }
          break;
        case 131:
          for (int i = min_range; i <= max_range; i++) {
            current_node = new Node;
            current_node->name = text();
            if (min_range != max_range)
              current_node->name += "["+std::to_string(i)+"]";
            circuit->all_nodes.push_back(current_node);
            circuit->all_nodes_map[current_node->name] = current_node;
          }
          break;
        case 145:
          current_node = circuit->GetNode(string(text()));
          break;
        case NODE_DFF:
          mode = 150;
          dff_d.clear();
          dff_q.clear();
          dff_qn.clear();
          break;
      }

{dotidempty}  mode++;

{dotid}       mode=mode;

{paridindex}  |
{parid} switch (mode) {
          case NODE_BUF:
          case NODE_NOT:
          case NODE_AND:
          case NODE_NAND:
          case NODE_OR:
          case NODE_NOR:
          case NODE_XOR:
          case NODE_XNOR:
            name = text();
            name.erase(name.begin());
            name.erase(name.end()-1);
            current_node = circuit->GetNode(name);
            current_node->type = (NodeType)mode;
            mode = 140;
            break;
          case 140:
            name = text();
            name.erase(name.begin());
            name.erase(name.end()-1);
            temp_node = circuit->GetNode(name);
            current_node->inputs.push_back(temp_node);
            temp_node->outputs.push_back(current_node);
            break;
          case 150:
            mode = 151;
            break;
          case 151:
            name = text();
            name.erase(name.begin());
            name.erase(name.end()-1);
            dff_qn = name;
            mode = 152;
            break;
          case 152:
            name = text();
            name.erase(name.begin());
            name.erase(name.end()-1);
            dff_q = name;
            mode = 153;
            break;
          case 153:
            name = text();
            name.erase(name.begin());
            name.erase(name.end()-1);
            dff_d = name;
            mode = 154;
            break;
          case 154:
            if (dff_q.size() == 0) {
              dff_q = dff_qn+"_not";
              current_node = new Node;
              current_node->name = dff_q;
              circuit->all_nodes.push_back(current_node);
              circuit->all_nodes_map[current_node->name] = current_node;
            }
            else {
              current_node = circuit->GetNode(dff_q);
            }
            current_node->type = NODE_DFF;
            circuit->ffs.push_back(current_node);
            temp_node = circuit->GetNode(dff_d);
            current_node->inputs.push_back(temp_node);
            temp_node->outputs.push_back(current_node);
            if (dff_qn.size() > 0) {
              temp_node = circuit->GetNode(dff_qn);
              temp_node->type = NODE_NOT;
              current_node->outputs.push_back(temp_node);
              temp_node->inputs.push_back(current_node);
            }
            mode = 254;
            break;
        }

{zero} if (mode == 145) {
         current_node->type = NODE_ZERO;
       }

{one}  if (mode == 145) {
          current_node->type = NODE_ONE;
       }

{parzero}  if (mode == 140) {
             temp_node = circuit->GetNode("zero_zero");
             if (!temp_node) {
               temp_node = new Node(NODE_ZERO);
               temp_node->name = "zero_zero";
               circuit->all_nodes.push_back(temp_node);
               circuit->all_nodes_map[temp_node->name] = temp_node;
             }
             current_node->inputs.push_back(temp_node);
             temp_node->outputs.push_back(current_node);
           }

{parone}   if (mode == 140) {
             temp_node = circuit->GetNode("one_one");
             if (!temp_node) {
               temp_node = new Node(NODE_ONE);
               temp_node->name = "one_one";
               circuit->all_nodes.push_back(temp_node);
               circuit->all_nodes_map[temp_node->name] = temp_node;
             }
             current_node->inputs.push_back(temp_node);
             temp_node->outputs.push_back(current_node);
           }

.  mode=mode;

%%
