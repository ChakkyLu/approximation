#include "glucose_interface.h"

#include <chrono>

using namespace Glucose;
using namespace std;


namespace satcnf {

  GInterface::GInterface() {
    //solver = new SimpSolver();
    solver = &mainsolver;
    InitSolver();
  }

  GInterface::GInterface(const GInterface &old) : mainsolver(*old.solver) {
    //solver = new SimpSolver(*old.solver);
    solver = &mainsolver;
    //InitSolver();
  }

  void GInterface::InitSolver() {
#if defined(__linux__)
    fpu_control_t oldcw, newcw;
    _FPU_GETCW(oldcw);
    newcw = (oldcw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
    _FPU_SETCW(newcw);
#endif
    solver->eliminate(false);
    solver->setIncrementalMode();
    solver->parsing = 1;
    solver->verbosity = 0;
    return;
  }

  Lit GInterface::GenLit(LiteralType x, long offset) {
    LiteralType abs_x = 0;
    if (x < 0)
      abs_x = -x + offset;
    else
      abs_x = x + offset;

    while (abs_x >= solver->nVars())
      solver->newVar();

    if (x > 0)
      return mkLit(abs_x);
    else
      return ~mkLit(abs_x);
  }

  void GInterface::AddCnf(ClauseVector &clauses, long offset) {
    for (auto cl = begin(clauses); cl != end(clauses); cl++) {
      vec<Lit> tmp;
      for (auto lit : (*cl)) {
        tmp.push(GenLit(lit, offset));
      }
      solver->addClause_(tmp);
    }
  }

  void GInterface::AddDimacs(const char *dimacsFile, long offset) {
    solver->parsing = 1;
    if (dimacsFile == NULL) {
      ERR("No input dimacs!");
      return;
    }
    gzFile fp = gzopen(dimacsFile, "rb");
    StreamBuffer in(fp);
    while (1) {
      skipWhitespace(in);
      if (*in == EOF)
        break;
      else if (*in == 'c' || *in == 'p')
        skipLine(in);
      else {
        vec<Lit> lits;
        ReadClause(in, lits, offset);
        solver->addClause_(lits);
      }
    }
    gzclose(fp);
  }

  void GInterface::SetEqual(LiteralVector &values1, LiteralVector &values2) {
    if (values1.size() != values2.size())
      return;
    for (int i = 0; i < values1.size(); i++) {
      SetEqual(values1[i], values2[i]);
    }
  }

  void GInterface::SetEqual(LiteralType a, LiteralType b) {
    if (a == b)
      return;
    vec<Lit> lits1, lits2;
    lits1.push(GenLit(a));
    lits1.push(GenLit(-b));
    lits2.push(GenLit(b));
    lits2.push(GenLit(-a));
    solver->addClause_(lits1);
    solver->addClause_(lits2);
  }


  void GInterface::SetNEqual(LiteralVector &values1, LiteralVector &values2) {
    if (values1.size() != values2.size())
      return;

    int max = solver->nVars();

    for (int i = 0; i < values1.size(); i++) {
      SetNEqual(values1[i], values2[i]);
    }

    vec<Lit> tmp;
    for (int i = 0; i < values1.size(); i++)
      tmp.push(GenLit(max + i));
    solver->addClause_(tmp);
  }

  void GInterface::SetNEqual(LiteralType a, LiteralType b) {
    if (a == b)
      return;

    int max = solver->nVars();

    vec<Lit> lits1, lits2, lits3, lits4;
    lits1.push(GenLit(-a));
    lits1.push(GenLit(-b));
    lits1.push(GenLit(-max));
    lits2.push(GenLit(-a));
    lits2.push(GenLit(b));
    lits2.push(GenLit(max));
    lits3.push(GenLit(a));
    lits3.push(GenLit(-b));
    lits3.push(GenLit(max));
    lits4.push(GenLit(a));
    lits4.push(GenLit(b));
    lits4.push(GenLit(-max));
    solver->addClause_(lits1);
    solver->addClause_(lits2);
    solver->addClause_(lits3);
    solver->addClause_(lits4);
  }

  void GInterface::SetEqualIfNotFlag(LiteralVector &values1, LiteralVector &values2, LiteralType flag) {
    if (values1.size() != values2.size())
      return;

    GenLit(flag);

    for (int i = 0; i < values1.size(); i++) {
      SetEqualIfNotFlag(values1[i], values2[i], flag);
    }
  }

  void GInterface::SetEqualIfNotFlag(LiteralType a, LiteralType b, LiteralType flag) {
    if (a == b)
      return;

    GenLit(flag);

    vec<Lit> lits1, lits2;
    lits1.push(GenLit(a));
    lits1.push(GenLit(-b));
    lits2.push(GenLit(b));
    lits2.push(GenLit(-a));

    lits1.push(GenLit(flag));
    lits2.push(GenLit(flag));

    solver->addClause_(lits1);
    solver->addClause_(lits2);
  }


  void GInterface::SetNEqualIfNotFlag(LiteralVector &values1, LiteralVector &values2, LiteralType flag) {
    if (values1.size() != values2.size())
      return;

    GenLit(flag);

    int max = solver->nVars();

    for (int i = 0; i < values1.size(); i++) {
      SetNEqualIfNotFlag(values1[i], values2[i], flag);
    }

    vec<Lit> tmp;
    for (int i = 0; i < values1.size(); i++)
      tmp.push(GenLit(max + i));
    tmp.push(GenLit(flag));
    solver->addClause_(tmp);
  }

  void GInterface::SetNEqualIfNotFlag(LiteralType a, LiteralType b, LiteralType flag) {
    if (a == b)
      return;

    GenLit(flag);

    int max = solver->nVars();

    vec<Lit> lits1, lits2, lits3, lits4;
    lits1.push(GenLit(-a));
    lits1.push(GenLit(-b));
    lits1.push(GenLit(-max));
    lits2.push(GenLit(-a));
    lits2.push(GenLit(b));
    lits2.push(GenLit(max));
    lits3.push(GenLit(a));
    lits3.push(GenLit(-b));
    lits3.push(GenLit(max));
    lits4.push(GenLit(a));
    lits4.push(GenLit(b));
    lits4.push(GenLit(-max));

    lits1.push(GenLit(flag));
    lits2.push(GenLit(flag));
    lits3.push(GenLit(flag));
    lits4.push(GenLit(flag));

    solver->addClause_(lits1);
    solver->addClause_(lits2);
    solver->addClause_(lits3);
    solver->addClause_(lits4);
  }

  void GInterface::IfTrueThenEQElseNEQ(LiteralVector &values1, LiteralVector &values2, LiteralType flag) {
    IfTrueThenEQ(values1, values2, flag);
    IfFalseThenNEQ(values1, values2, flag);
  }

  void GInterface::IfTrueThenNEQElseEQ(LiteralVector &values1, LiteralVector &values2, LiteralType flag) {
    IfTrueThenNEQ(values1, values2, flag);
    IfFalseThenEQ(values1, values2, flag);
  }

  void GInterface::IfTrueThenNEQ(LiteralVector &values1, LiteralVector &values2, LiteralType flag) {
    return SetNEqualIfNotFlag(values1, values2, -flag);
  }

  void GInterface::IfFalseThenNEQ(LiteralVector &values1, LiteralVector &values2, LiteralType flag) {
    return SetNEqualIfNotFlag(values1, values2, flag);
  }

  void GInterface::IfFalseThenEQ(LiteralVector &values1, LiteralVector &values2, LiteralType flag) {
    return SetEqualIfNotFlag(values1, values2, flag);
  }

  void GInterface::IfTrueThenEQ(LiteralVector &values1, LiteralVector &values2, LiteralType flag) {
    return SetEqualIfNotFlag(values1, values2, -flag);
  }

  lbool GInterface::ModelNum(int num) {
    return solver->model[num];
  }

  void GInterface::ReadClause(StreamBuffer &in, vec<Lit> &lits, long offset) {
    lits.clear();
    int parsed_lit = parseInt(in);
    while (parsed_lit != 0) {
      lits.push(GenLit(parsed_lit, offset));
      parsed_lit = parseInt(in);
    }
  }

  bool GInterface::GetAnswer(const bool VERBOSE) {
    vec<Lit> dummy;
    std::chrono::high_resolution_clock::time_point start, end;
    solver->parsing = 0;
    solver->verbosity = 0;

    if (VERBOSE)
      start = std::chrono::high_resolution_clock::now();
    lbool res = solver->solveLimited(dummy);
    if (VERBOSE) {
      end = std::chrono::high_resolution_clock::now();
      double time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    }

    solver->parsing = 1;

    if (VERBOSE)
      std::cout << "#SOLVER vars:" << solver->nVars() << ", clauses: "
                << solver->nClauses() << ", time: " << time << "ms" << std::endl;

    return (res == l__True) ? true : false;
  }

  bool GInterface::GetAnswer(const vec<Lit> &assump, const bool VERBOSE) {
    std::chrono::high_resolution_clock::time_point start, end;
    solver->parsing = 0;
    solver->verbosity = 0;

    if (VERBOSE) {
      start = std::chrono::high_resolution_clock::now();
    }
    lbool res = solver->solveLimited(assump);
    if (VERBOSE) {
      end = std::chrono::high_resolution_clock::now();
      double time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    }
    solver->parsing = 1;

    if (VERBOSE)
      std::cout << "#SOLVER vars:" << solver->nVars() << ", clauses: "
                << solver->nClauses() << ", time: " << time << "ms" << std::endl;

    return (res == l__True) ? true : false;
  }

}