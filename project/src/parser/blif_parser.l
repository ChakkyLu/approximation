// simple blif parser, without error recover, and without comment detection

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
%}

%class{
Node* current_node;
NodeVector cur_nodes;
int mode;
virtual int wrap() { return 1; }
public:
Circuit* circuit;
%}

%init{
circuit = NULL;
current_node = NULL;
mode = 255; // 128:start 129:output 130:input 131:gate 132:value 254:nothing 255:finish
%}

newline  \r?\n
space [ ]
id  [a-zA-Z][A-Za-z0-9_]*("["[0-9]+"]")*
vals [01-][01-]+
bits [01]

val {bits}{newline}|{vals}{space}{bits}{newline}|{bits}{space}{bits}{newline}

%option fast
%option namespace=blifparser
%option lexer=BL_LEXER
%option prefix=bl_

%%

".end"      mode = 255;
".model"    mode = 128;
".outputs"  mode = 129;
".inputs"   mode = 130;
".names"    {
              mode = 131;
              cur_nodes.clear();
            }


{id}	switch (mode) {
          case 128:
            circuit->name = text();
            mode = 254;
            break;
          case 129:
            current_node = new BlifNode;
            current_node->name = text();
            current_node->is_output = true;
            circuit->outputs.push_back(current_node);
            circuit->all_nodes.push_back(current_node);
            circuit->all_nodes_map[current_node->name] = current_node;
            break;
          case 130:
            current_node = new BlifNode;
            current_node->name = text();
            current_node->is_input = true;
            circuit->inputs.push_back(current_node);
            circuit->all_nodes.push_back(current_node);
            circuit->all_nodes_map[current_node->name] = current_node;
            break;
          case 131:
            current_node = circuit->GetNode(string(text()));
            if (!current_node) {
              current_node = new BlifNode;
              current_node->name = text();
              circuit->all_nodes.push_back(current_node);
              circuit->all_nodes_map[current_node->name] = current_node;
            }
            cur_nodes.push_back(current_node);
            break;
        }

{val}   if (cur_nodes.size() == 0)
          mode = 254; 
        else {
          if (mode == 131) {
            current_node = cur_nodes[cur_nodes.size()-1];
            for (int i = 0; i < cur_nodes.size()-1; i++) {
              Node* node = cur_nodes[i];
              current_node->inputs.push_back(node);
              node->outputs.push_back(current_node);
            }
            mode = 132;
          }
          if (mode = 132) {
            string str(text());
            while (str[str.size()-1] == '\n' || str[str.size()-1] == '\r')
              str.resize(str.size()-1);
            ((BlifNode*)current_node)->str_values.push_back(str);
          }
        }

.	mode=mode;

%%
