#ifndef _NODE_H_INCLUDED
#define _NODE_H_INCLUDED

#include "global.h"

namespace nodecircuit {

  enum NodeType {
    NODE_UNKNOWN,
    NODE_ZERO,
    NODE_ONE,
    NODE_BUF,
    NODE_NOT,
    NODE_AND,
    NODE_NAND,
    NODE_OR,
    NODE_NOR,
    NODE_XOR,
    NODE_XNOR,
    NODE_AND2_NP,
    NODE_NAND2_NP,
    NODE_AND2_PN,
    NODE_NAND2_PN,
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
