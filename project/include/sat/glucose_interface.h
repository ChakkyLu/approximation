#ifndef GLUCOSE_INTERFACE_H_
#define GLUCOSE_INTERFACE_H_

#include "global.h"
#include "cnf.h"

#include "utils/System.h"
#include "core/Dimacs.h"
#include "simp/SimpSolver.h"
#include "core/SolverTypes.h"
#include "mtl/Vec.h"


namespace satcnf {

  class GInterface {
  public:
    GInterface();
    GInterface(const GInterface &obj);

    virtual ~GInterface() {};

  public:
    void AddDimacs(const char *dimacsFile, long offset);
    void AddCnf(ClauseVector &clauses, long offset);

    bool GetAnswer(bool VERBOSE = false);
    bool GetAnswer(const Glucose::vec<Glucose::Lit> &assump, bool VERBOSE = false);

    void SetEqual(LiteralType a, LiteralType b);
    void SetNEqual(LiteralType a, LiteralType b);
    void SetEqual(LiteralVector &values1, LiteralVector &values2);
    void SetNEqual(LiteralVector &values1, LiteralVector &values2);

    void SetEqualIfNotFlag(LiteralType a, LiteralType b, LiteralType flag);
    void SetNEqualIfNotFlag(LiteralType a, LiteralType b, LiteralType flag);
    void SetEqualIfNotFlag(LiteralVector &values1, LiteralVector &values2, LiteralType flag);
    void SetNEqualIfNotFlag(LiteralVector &values1, LiteralVector &values2, LiteralType flag);

    void IfTrueThenEQElseNEQ(LiteralVector &values1, LiteralVector &values2, LiteralType flag);
    void IfTrueThenNEQElseEQ(LiteralVector &values1, LiteralVector &values2, LiteralType flag);

    void IfTrueThenEQ(LiteralVector &values1, LiteralVector &values2, LiteralType flag);
    void IfFalseThenEQ(LiteralVector &values1, LiteralVector &values2, LiteralType flag);
    void IfTrueThenNEQ(LiteralVector &values1, LiteralVector &values2, LiteralType flag);
    void IfFalseThenNEQ(LiteralVector &values1, LiteralVector &values2, LiteralType flag);

    Glucose::lbool ModelNum(int num);
    Glucose::Lit GenLit(LiteralType x, long offset = 0);

  public:
    Glucose::SimpSolver *solver;
    Glucose::SimpSolver mainsolver;

  private:
    void InitSolver();
    void ReadClause(Glucose::StreamBuffer &in, Glucose::vec<Glucose::Lit> &lits, long offset);

  };

}

#endif
