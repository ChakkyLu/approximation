#ifndef _NODE_H_INCLUDED
#define _NODE_H_INCLUDED

#include "global.h"

namespace nodecircuit {

  enum NodeType {
    NODE_UNKNOWN,
    NODE_ZERO,//1 0
    NODE_ONE,//2 1
    NODE_BUF,//3 2,3
    NODE_NOT,//4 4,5
    NODE_AND,//5 6
    NODE_NAND,//6 7
    NODE_OR,//7 8
    NODE_NOR,//8 9
    NODE_XOR,//9 10
    NODE_XNOR,//10 11
    NODE_AND2_NP,//11 12
    NODE_NAND2_NP,//12 13
    NODE_AND2_PN,//13 14
    NODE_NAND2_PN,//14 15
    NODE_BLIF,
    NODE_DFF
  };

  class Node;

  typedef std::vector<Node *> NodeVector;
  typedef std::set<Node *> NodeSet;

  class Node {
  public:
    Node(NodeType _type = NODE_UNKNOWN) {
      type = _type;
      is_input = false;
      is_output = false;
      level = -1; // not levelized yet!
      index = -1;
    }

    virtual ~Node() {}

  public:
    NodeType type;
    std::string name;
    bool is_input;
    bool is_output;
    NodeVector inputs;  // fanins  of the node
    NodeVector outputs; // fanouts of the node
    long level; // topological level
    long index; // index in the array of all nodes in the circuit

  };

#define CODE_ZERO 1
#define CODE_ONE  2
#define CODE_DC   3

  class BlifNode : public Node {
  public:
    BlifNode() : Node(NODE_BLIF) {
      result_is_one = false;
    }

    virtual ~BlifNode() {}

  public:
    std::vector<std::string> str_values;
    bool result_is_one; // or it is zero!
    // 2-bit codes: '0' -> 01 , '1' -> 10 , '-' -> 11
    std::set<uint64_t> coded_values;

  public:
    // convert the strings to the coded values
    // maximum number of inputs is 32 (32*2 = 64 bits)
    int CreateCodedValues();

    // check if the Blif corresponds to a standard gate type, like and, nor, xor,...
    NodeType GetEquivalentType();
  };

}

#endif // _NODE_H_INCLUDED
