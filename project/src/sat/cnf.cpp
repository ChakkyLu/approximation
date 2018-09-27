#include "cnf.h"

namespace satcnf {

  using namespace std;
  using namespace nodecircuit;

  // TODO: sequential circuits not supported -> error message?


  int Cnf::Convert2Cnf(Circuit *circuit, NodeVector *ex_ins, NodeVector *ex_outs, bool optimize) {
    // TODO: currently only optimized version is implemented!
    org_circuit = circuit;
    clauses.clear();
    literal_str_map.clear();
    literal_node_map.clear();

    ABC::Abc_Ntk_t *pNtk, *pTemp;
    ABC::Abc_Frame_t *pAbc;

    pNtk = ABC::Abc_NtkStartRead(strdup(circuit->name.c_str()));
    Circuit2Abc(pNtk, ex_ins, ex_outs);

    pAbc = ABC::Abc_FrameGetGlobalFrame();
    Abc_FrameReplaceCurrentNetwork(pAbc, pNtk);
    ABC::Abc_NtkFinalizeRead(pNtk);
    if (ABC::Abc_NtkCheckRead(pNtk) == 0) {
      ERR("Convert2Cnf:: Circuit is wrong");
    }
    pNtk = ABC::Abc_NtkToLogic(pTemp = pNtk);

    ABC::Abc_FrameReplaceCurrentNetwork(pAbc, pNtk);
    ABC::Abc_FrameClearVerifStatus(pAbc);
    std::string Command;
    if (optimize)
      Command = "balance -l; rewrite -l; refactor -l; balance -l; rewrite -l; rewrite -z -l; balance -l; refactor -z -l; rewrite -z -l; balance -l;";
    else
      Command = "strash";
    if (ABC::Cmd_CommandExecute(pAbc, Command.c_str())) {
      ERR(string("ABC:: Command execution error: ") + Command);
      return 1;
    }

    pTemp = Abc_NtkToNetlist(pAbc->pNtkCur);
    if (!ABC::Abc_NtkHasSop(pTemp) && !ABC::Abc_NtkHasMapping(pTemp))
      ABC::Abc_NtkToSop(pTemp, -1, ABC_INFINITY);

    Abc2Cnf(pTemp);

    return 0;
  }

  int Cnf::Convert2Cnf(Circuit *circuit, bool optimize) {
    // TODO: currently only optimized version is implemented!
    org_circuit = circuit;
    clauses.clear();
    literal_str_map.clear();
    literal_node_map.clear();

    ABC::Abc_Ntk_t *pNtk, *pTemp;
    ABC::Abc_Frame_t *pAbc;

    pNtk = ABC::Abc_NtkStartRead(strdup(circuit->name.c_str()));
    Circuit2Abc(pNtk);

    pAbc = ABC::Abc_FrameGetGlobalFrame();
    Abc_FrameReplaceCurrentNetwork(pAbc, pNtk);
    ABC::Abc_NtkFinalizeRead(pNtk);
    if (ABC::Abc_NtkCheckRead(pNtk) == 0) {
      ERR("Convert2Cnf:: Circuit is wrong");
    }
    pNtk = ABC::Abc_NtkToLogic(pTemp = pNtk);

    ABC::Abc_FrameReplaceCurrentNetwork(pAbc, pNtk);
    ABC::Abc_FrameClearVerifStatus(pAbc);
    char Command[300] = "balance -l; rewrite -l; refactor -l; balance -l; rewrite -l; rewrite -z -l; balance -l; refactor -z -l; rewrite -z -l; balance -l;";
    //char Command[100] = "strash;dc2;dc2;";
    if (ABC::Cmd_CommandExecute(pAbc, Command)) {
      ERR(string("ABC:: Command execution error: ") + Command);
      return 1;
    }

    pTemp = Abc_NtkToNetlist(pAbc->pNtkCur);
    if (!ABC::Abc_NtkHasSop(pTemp) && !ABC::Abc_NtkHasMapping(pTemp))
      ABC::Abc_NtkToSop(pTemp, -1, ABC_INFINITY);

    Abc2Cnf(pTemp);

    return 0;
  }

  int Cnf::Circuit2Abc(ABC::Abc_Ntk_t *pNtk, NodeVector *ex_ins, NodeVector *ex_outs) {

    for (auto input = begin(org_circuit->inputs); input != end(org_circuit->inputs); input++) {
      Io_ReadCreatePi(pNtk, strdup(((*input)->name).c_str()));
    }
    for (auto output = begin(org_circuit->outputs); output != end(org_circuit->outputs); output++) {
      Io_ReadCreatePo(pNtk, strdup(((*output)->name).c_str()));
    }

    if (ex_ins)
      for (auto tmp = (*ex_ins).begin(); tmp != (*ex_ins).end(); tmp++) {
        if ((*tmp)->is_input)
          continue;
        Io_ReadCreatePi(pNtk, strdup(((*tmp)->name).c_str()));
      }
    if (ex_outs)
      for (auto tmp = (*ex_outs).begin(); tmp != (*ex_outs).end(); tmp++) {
        if ((*tmp)->is_output)
          continue;
        Io_ReadCreatePo(pNtk, strdup(((*tmp)->name).c_str()));
      }

    for (auto node = begin(org_circuit->all_nodes); node != end(org_circuit->all_nodes); node++) {
      bool ignoreFlag = false;
      if (ex_ins)
        for (auto exin : *ex_ins) {
          if ((*node)->name == exin->name && (*node)->inputs.size() != 0) {
            ignoreFlag = true;
          }
        }
      for (auto in : org_circuit->inputs) {
        if ((*node)->name == in->name && (*node)->inputs.size() != 0) {
          ignoreFlag = true;
        }
      }
      if (ignoreFlag) {
        continue;
      }

      if ((*node)->type == NODE_UNKNOWN ||
          (*node)->type >= NODE_BLIF) // TODO: BLIF node is not supported yet!
        continue;

      if ((*node)->type == NODE_ZERO) {
        ABC::Abc_Obj_t *pNode = Io_ReadCreateNode(pNtk, strdup((*node)->name.c_str()), NULL, 0);
        ABC::Abc_ObjSetData(pNode, ABC::Abc_SopCreateConst0((ABC::Mem_Flex_t *) pNtk->pManFunc));
        continue;
      }

      if ((*node)->type == NODE_ONE) {
        ABC::Abc_Obj_t *pNode = Io_ReadCreateNode(pNtk, strdup((*node)->name.c_str()), NULL, 0);
        ABC::Abc_ObjSetData(pNode, ABC::Abc_SopCreateConst1((ABC::Mem_Flex_t *) pNtk->pManFunc));
        continue;
      }

      bool FLAG = true;
      if ((*node)->is_input) {
        continue;
      }

      long nNames = (*node)->inputs.size();

      char **pNames = (char **) malloc(sizeof(char *) * nNames);
      for (long i = 0; i < nNames; i++) {
        pNames[i] = strdup((*node)->inputs[i]->name.c_str());
      }

      if (ex_ins)
        for (auto tmp = (*ex_ins).begin(); FLAG && tmp != (*ex_ins).end(); tmp++) {
          FLAG = ((*tmp) != (*node));
        }

      if (FLAG) {
        ABC::Abc_Obj_t *pNode = NULL;
        if ((*node)->type == nodecircuit::NODE_NAND2_NP || (*node)->type == nodecircuit::NODE_NAND2_PN)
          pNode = Io_ReadCreateNode(pNtk, strdup(string((*node)->name+"temp").c_str()), pNames, nNames);
        else
          pNode = Io_ReadCreateNode(pNtk, strdup((*node)->name.c_str()), pNames, nNames);

        switch ((*node)->type) {
          case nodecircuit::NODE_AND:
            Abc_ObjSetData(pNode, ABC::Abc_SopCreateAnd((ABC::Mem_Flex_t *) pNtk->pManFunc, nNames, NULL));
            break;
          case nodecircuit::NODE_OR:
            Abc_ObjSetData(pNode, ABC::Abc_SopCreateOr((ABC::Mem_Flex_t *) pNtk->pManFunc, nNames, NULL));
            break;
          case nodecircuit::NODE_NAND:
            Abc_ObjSetData(pNode, ABC::Abc_SopCreateNand((ABC::Mem_Flex_t *) pNtk->pManFunc, nNames));
            break;
          case nodecircuit::NODE_NOR:
            Abc_ObjSetData(pNode, ABC::Abc_SopCreateNor((ABC::Mem_Flex_t *) pNtk->pManFunc, nNames));
            break;
          case nodecircuit::NODE_XOR:
            Abc_ObjSetData(pNode, ABC::Abc_SopCreateXor((ABC::Mem_Flex_t *) pNtk->pManFunc, nNames));
            break;
          case nodecircuit::NODE_XNOR:
            Abc_ObjSetData(pNode, ABC::Abc_SopCreateNxor((ABC::Mem_Flex_t *) pNtk->pManFunc, nNames));
            break;
          case nodecircuit::NODE_BUF:
            Abc_ObjSetData(pNode, ABC::Abc_SopCreateBuf((ABC::Mem_Flex_t *) pNtk->pManFunc));
            break;
          case nodecircuit::NODE_NOT:
            Abc_ObjSetData(pNode, ABC::Abc_SopCreateInv((ABC::Mem_Flex_t *) pNtk->pManFunc));
            break;
          case nodecircuit::NODE_AND2_NP:
            Abc_ObjSetData(pNode, ABC::Abc_SopCreateAnd2((ABC::Mem_Flex_t *) pNtk->pManFunc, 1, 0));
            break;
          case nodecircuit::NODE_AND2_PN:
            Abc_ObjSetData(pNode, ABC::Abc_SopCreateAnd2((ABC::Mem_Flex_t *) pNtk->pManFunc, 0, 1));
            break;
          case nodecircuit::NODE_NAND2_NP:
            Abc_ObjSetData(pNode, ABC::Abc_SopCreateAnd2((ABC::Mem_Flex_t *) pNtk->pManFunc, 1, 0));
            {
              char **pNames2 = (char **) malloc(sizeof(char *) * 1);
              pNames2[0] = strdup(string((*node)->name + "temp").c_str());
              ABC::Abc_Obj_t *pNode2 = Io_ReadCreateNode(pNtk, strdup((*node)->name.c_str()), pNames2, 1);
              Abc_ObjSetData(pNode2, ABC::Abc_SopCreateInv((ABC::Mem_Flex_t *) pNtk->pManFunc));
            }
            break;
          case nodecircuit::NODE_NAND2_PN:
            Abc_ObjSetData(pNode, ABC::Abc_SopCreateAnd2((ABC::Mem_Flex_t *) pNtk->pManFunc, 0, 1));
            {
              char **pNames2 = (char **) malloc(sizeof(char *) * 1);
              pNames2[0] = strdup(string((*node)->name + "temp").c_str());
              ABC::Abc_Obj_t *pNode2 = Io_ReadCreateNode(pNtk, strdup((*node)->name.c_str()), pNames2, 1);
              Abc_ObjSetData(pNode2, ABC::Abc_SopCreateInv((ABC::Mem_Flex_t *) pNtk->pManFunc));
            }
            break;
          default:
            ERR("Circuit2Abc:: Unexpected Node! "+(*node)->name);
            return 1;
        }
      }
    }
    return 0;
  }

  int Cnf::Circuit2Abc(ABC::Abc_Ntk_t *pNtk) {

    for (auto input = begin(org_circuit->inputs); input != end(org_circuit->inputs); input++) {
      Io_ReadCreatePi(pNtk, strdup(((*input)->name).c_str()));
    }
    for (auto output = begin(org_circuit->outputs); output != end(org_circuit->outputs); output++) {
      Io_ReadCreatePo(pNtk, strdup(((*output)->name).c_str()));
    }

    for (auto node = begin(org_circuit->all_nodes); node != end(org_circuit->all_nodes); node++) {
      bool ignoreFlag = false;
      for (auto in : org_circuit->inputs) {
        if ((*node)->name == in->name && (*node)->inputs.size() != 0) {
          ignoreFlag = true;
        }
      }
      if (ignoreFlag) {
        continue;
      }

      if ((*node)->type == NODE_UNKNOWN ||
          (*node)->type >= NODE_BLIF) // TODO: BLIF node is not supported yet!
        continue;

      if ((*node)->is_input) {
        continue;
      }

      if ((*node)->type == NODE_ONE) {
        ABC::Abc_Obj_t *pNode = Io_ReadCreateNode(pNtk, strdup((*node)->name.c_str()), NULL, 0);
        ABC::Abc_ObjSetData(pNode, ABC::Abc_SopCreateConst1((ABC::Mem_Flex_t *) pNtk->pManFunc));
        continue;
      }

      if ((*node)->type == NODE_ZERO) {
        ABC::Abc_Obj_t *pNode = Io_ReadCreateNode(pNtk, strdup((*node)->name.c_str()), NULL, 0);
        ABC::Abc_ObjSetData(pNode, ABC::Abc_SopCreateConst0((ABC::Mem_Flex_t *) pNtk->pManFunc));
        continue;
      }

      long nNames = (*node)->inputs.size();

      char **pNames = (char **) malloc(sizeof(char *) * nNames);
      for (long i = 0; i < nNames; i++) {
        pNames[i] = strdup((*node)->inputs[i]->name.c_str());
      }

      ABC::Abc_Obj_t *pNode = NULL;
      if ((*node)->type == nodecircuit::NODE_NAND2_NP || (*node)->type == nodecircuit::NODE_NAND2_PN)
        pNode = Io_ReadCreateNode(pNtk, strdup(string((*node)->name+"temp").c_str()), pNames, nNames);
      else
        pNode = Io_ReadCreateNode(pNtk, strdup((*node)->name.c_str()), pNames, nNames);

      switch ((*node)->type) {
        case nodecircuit::NODE_AND:
          Abc_ObjSetData(pNode, ABC::Abc_SopCreateAnd((ABC::Mem_Flex_t *) pNtk->pManFunc, nNames, NULL));
          break;
        case nodecircuit::NODE_OR:
          Abc_ObjSetData(pNode, ABC::Abc_SopCreateOr((ABC::Mem_Flex_t *) pNtk->pManFunc, nNames, NULL));
          break;
        case nodecircuit::NODE_NAND:
          Abc_ObjSetData(pNode, ABC::Abc_SopCreateNand((ABC::Mem_Flex_t *) pNtk->pManFunc, nNames));
          break;
        case nodecircuit::NODE_NOR:
          Abc_ObjSetData(pNode, ABC::Abc_SopCreateNor((ABC::Mem_Flex_t *) pNtk->pManFunc, nNames));
          break;
        case nodecircuit::NODE_XOR:
          Abc_ObjSetData(pNode, ABC::Abc_SopCreateXor((ABC::Mem_Flex_t *) pNtk->pManFunc, nNames));
          break;
        case nodecircuit::NODE_XNOR:
          Abc_ObjSetData(pNode, ABC::Abc_SopCreateNxor((ABC::Mem_Flex_t *) pNtk->pManFunc, nNames));
          break;
        case nodecircuit::NODE_BUF:
          Abc_ObjSetData(pNode, ABC::Abc_SopCreateBuf((ABC::Mem_Flex_t *) pNtk->pManFunc));
          break;
        case nodecircuit::NODE_NOT:
          Abc_ObjSetData(pNode, ABC::Abc_SopCreateInv((ABC::Mem_Flex_t *) pNtk->pManFunc));
          break;
        case nodecircuit::NODE_AND2_NP:
          Abc_ObjSetData(pNode, ABC::Abc_SopCreateAnd2((ABC::Mem_Flex_t *) pNtk->pManFunc, 1, 0));
          break;
        case nodecircuit::NODE_AND2_PN:
          Abc_ObjSetData(pNode, ABC::Abc_SopCreateAnd2((ABC::Mem_Flex_t *) pNtk->pManFunc, 0, 1));
          break;
        case nodecircuit::NODE_NAND2_NP:
          Abc_ObjSetData(pNode, ABC::Abc_SopCreateAnd2((ABC::Mem_Flex_t *) pNtk->pManFunc, 1, 0));
          {
            char **pNames2 = (char **) malloc(sizeof(char *) * 1);
            pNames2[0] = strdup(string((*node)->name + "temp").c_str());
            ABC::Abc_Obj_t *pNode2 = Io_ReadCreateNode(pNtk, strdup((*node)->name.c_str()), pNames2, 1);
            Abc_ObjSetData(pNode2, ABC::Abc_SopCreateInv((ABC::Mem_Flex_t *) pNtk->pManFunc));
          }
          break;
        case nodecircuit::NODE_NAND2_PN:
          Abc_ObjSetData(pNode, ABC::Abc_SopCreateAnd2((ABC::Mem_Flex_t *) pNtk->pManFunc, 0, 1));
          {
            char **pNames2 = (char **) malloc(sizeof(char *) * 1);
            pNames2[0] = strdup(string((*node)->name + "temp").c_str());
            ABC::Abc_Obj_t *pNode2 = Io_ReadCreateNode(pNtk, strdup((*node)->name.c_str()), pNames2, 1);
            Abc_ObjSetData(pNode2, ABC::Abc_SopCreateInv((ABC::Mem_Flex_t *) pNtk->pManFunc));
          }
          break;
        default:
          ERR("Circuit2Abc:: Unexpected Node! "+(*node)->name);
          return 1;
      }
    }

    return 0;
  }

  int Cnf::Abc2Cnf(ABC::Abc_Ntk_t *pNtk) {
    LiteralType index = 1;

    string name;
    ABC::Abc_Obj_t *pNode, *pNet;
    long i, j;
    Abc_NtkForEachPi(pNtk, pNode, i) {
      name = Abc_ObjName(Abc_ObjFanout0(pNode));
      if (literal_str_map.find(name) == literal_str_map.end()) {
        literal_node_map[org_circuit->all_nodes_map[name]] = index;
        literal_str_map[name] = index;
        index++;
      }
    }
    Abc_NtkForEachPo(pNtk, pNode, i) {
      name = Abc_ObjName(Abc_ObjFanin0(pNode));
      if (literal_str_map.find(name) == literal_str_map.end()) {
        literal_node_map[org_circuit->all_nodes_map[name]] = index;
        literal_str_map[name] = index;
        index++;
      }
    }
    Abc_NtkForEachLatch(pNtk, pNode, i)ERR("Abc2Cnf:: Latch isn't supported");

    Abc_NtkForEachNode(pNtk, pNode, i) {
        char *pName;
        long nFanins = ABC::Abc_ObjFaninNum(pNode);
        LiteralType iOut;
        LiteralVector iIn;

        Abc_ObjForEachFanin(pNode, pNet, j) {
          pName = Abc_ObjName(pNet);
          name = pName;
          if (literal_str_map.find(name) == literal_str_map.end()) {
            literal_str_map[name] = index;
            iIn.push_back(index);
            index++;
          }
          else {
            iIn.push_back(literal_str_map[name]);
          }
        }
        pName = Abc_ObjName(Abc_ObjFanout0(pNode));
        name = pName;
        if (literal_str_map.find(name) == literal_str_map.end()) {
          literal_str_map[name] = index;
          iOut = index;
          index++;
        }
        else {
          iOut = literal_str_map[name];
        }

        LiteralVector lits;
        if (nFanins == 0) {
          if (Abc_NodeIsConst1(pNode)) {
            lits.push_back(iOut);
          }
          else {
            lits.push_back(-iOut);
          }
          clauses.push_back(lits);
          lits.clear();
        }
        else if (nFanins == 1) {
          if (Abc_NodeIsBuf(pNode)) {
            //BUFF;
            lits.push_back(iOut);
            lits.push_back(-iIn[0]);
            clauses.push_back(lits);
            lits.clear();
            lits.push_back(-iOut);
            lits.push_back(iIn[0]);
            clauses.push_back(lits);
            lits.clear();
          }
          else {
            //NOT;
            lits.push_back(-iOut);
            lits.push_back(-iIn[0]);
            clauses.push_back(lits);
            lits.clear();
            lits.push_back(iOut);
            lits.push_back(iIn[0]);
            clauses.push_back(lits);
            lits.clear();
          }
        }
        else if (nFanins == 2) {
          vector<bool> logic;
          string buffer = (char *) Abc_ObjData(pNode);

          logic.push_back(buffer[0] - '0');
          logic.push_back(buffer[1] - '0');
          logic.push_back(buffer[3] - '0');

          if (!logic[0])
            iIn[0] = -iIn[0];
          if (!logic[1])
            iIn[1] = -iIn[1];
          if (!logic[2])
            iOut = -iOut;

          lits.push_back(iOut);
          lits.push_back(-iIn[0]);
          lits.push_back(-iIn[1]);
          clauses.push_back(lits);
          lits.clear();
          lits.push_back(-iOut);
          lits.push_back(iIn[0]);
          clauses.push_back(lits);
          lits.clear();
          lits.push_back(-iOut);
          lits.push_back(iIn[1]);
          clauses.push_back(lits);
          lits.clear();
        }
        else {
          ERR("Abc2Cnf:: Not supported gate!");
        }

      }

    return 0;
  }

  void Cnf::WriteDimacs(string filename) {
    ofstream dimacs;
    dimacs.open(filename);

    dimacs << "p cnf " << literal_str_map.size() << " " << clauses.size() << endl;

    for (auto sig : literal_str_map)
      dimacs << "c " << sig.first << " " << sig.second << endl;

    for (auto itr : clauses) {
      for (auto lit: itr)
        dimacs << lit << " ";
      dimacs << "0" << endl;
    }
    dimacs.close();
  }
}