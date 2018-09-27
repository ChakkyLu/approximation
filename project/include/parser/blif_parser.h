// ./include/parser/blif_parser.h generated by reflex 1.0.9 from ./src/parser/blif_parser.l

#ifndef bl_REFLEX_HEADER_H
#define bl_REFLEX_HEADER_H
#define bl_IN_HEADER 1

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  OPTIONS USED                                                              //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#define bl_REFLEX_OPTION_fast                true
#define bl_REFLEX_OPTION_header_file         ./include/parser/blif_parser.h
#define bl_REFLEX_OPTION_lex                 bl_lex
#define bl_REFLEX_OPTION_lexer               BL_LEXER
#define bl_REFLEX_OPTION_namespace           blifparser
#define bl_REFLEX_OPTION_outfile             ./src/parser/blif_parser.cpp
#define bl_REFLEX_OPTION_prefix              bl_

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  SECTION 1: %top{ user code %}                                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#line 12 "./src/parser/blif_parser.l"

#include "circuit.h"
using namespace nodecircuit;


////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  REGEX MATCHER                                                             //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <reflex/matcher.h>

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  ABSTRACT LEXER CLASS                                                      //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

#include <reflex/abslexer.h>

////////////////////////////////////////////////////////////////////////////////
//                                                                            //
//  LEXER CLASS                                                               //
//                                                                            //
////////////////////////////////////////////////////////////////////////////////

namespace blifparser {

class BL_LEXER : public reflex::AbstractLexer<reflex::Matcher> {
#line 17 "./src/parser/blif_parser.l"

Node* current_node;
NodeVector cur_nodes;
int mode;
virtual int wrap() { return 1; }
public:
Circuit* circuit;

 public:
  typedef reflex::AbstractLexer<reflex::Matcher> AbstractBaseLexer;
  BL_LEXER(
      const reflex::Input& input = reflex::Input(),
      std::ostream&        os    = std::cout)
    :
      AbstractBaseLexer(input, os)
  {
#line 26 "./src/parser/blif_parser.l"

circuit = NULL;
current_node = NULL;
mode = 255; // 128:start 129:output 130:input 131:gate 132:value 254:nothing 255:finish

  }
  static const int INITIAL = 0;
  virtual int bl_lex();
  int bl_lex(
      const reflex::Input& input,
      std::ostream        *os = NULL)
  {
    in(input);
    if (os)
      out(*os);
    return bl_lex();
  }
};

} // namespace blifparser

#endif
