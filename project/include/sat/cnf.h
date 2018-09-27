#ifndef _CNF_H_INCLUDED
#define _CNF_H_INCLUDED

#include "global.h"
#include "circuit.h"

// ABC headers
#include "base/abc/abc.h"
#include "base/main/main.h"
#include "base/main/mainInt.h"

namespace satcnf {

  typedef long LiteralType;
  typedef std::vector<LiteralType> LiteralVector;
  typedef std::vector<LiteralVector> ClauseVector;

  class Cnf {
  public:
    Cnf() { org_circuit = NULL; }

    virtual ~Cnf() {};

  public:
    int Convert2Cnf(nodecircuit::Circuit *circuit, bool optimize = true);
    int Convert2Cnf(nodecircuit::Circuit *circuit, nodecircuit::NodeVector *ex_ins, nodecircuit::NodeVector *ex_outs, bool optimize = true);
    void WriteDimacs(std::string filename);

  protected:
    int Circuit2Abc(ABC::Abc_Ntk_t *pNtk);
    int Circuit2Abc(ABC::Abc_Ntk_t *pNtk, nodecircuit::NodeVector *ex_ins, nodecircuit::NodeVector *ex_outs);
    int Abc2Cnf(ABC::Abc_Ntk_t *pNtk);

  public:
    nodecircuit::Circuit *org_circuit;
    std::map<nodecircuit::Node *, LiteralType> literal_node_map;
    std::unordered_map<std::string, LiteralType> literal_str_map;
    ClauseVector clauses;
    LiteralType GetLiteral(std::string name) {
      std::unordered_map<std::string, LiteralType>::iterator it = literal_str_map.find(name);
      if (it == literal_str_map.end())
        return 0;
      return it->second;
    }
    LiteralType GetLiteral(nodecircuit::Node *node) {
      std::map<nodecircuit::Node *, LiteralType>::iterator it = literal_node_map.find(node);
      if (it == literal_node_map.end())
        return 0;
      return it->second;
    }
  };

}

#endif // _CNF_H_INCLUDED