#include "node.h"

using namespace std;

namespace nodecircuit {

  int BlifNode::CreateCodedValues() {
    int res = 0;
    if (is_input) {
      return res;
    }
    if (inputs.size() == 0) {
      if (str_values.size() != 1) {
        ERR("Expecting only one line for value! "+name);
        res = 1;
      }
      else {
        string cur_str = str_values[0];
        if (str_values.size() != 1) {
          ERR("expecting only one value! "+name);
          res = 1;
        }
        else if (str_values[0] == "0") {
          type = NODE_ZERO;
        }
        else if (str_values[0] == "1") {
          type = NODE_ONE;
        }
        else {
          ERR("unknown output value! "+name);
          res = 1;
        }
      }
      return res;
    }
    for (long i = 0; i < str_values.size(); i++) {
      string cur_str = str_values[i];
      uint64_t cur_val = 0;
      if (cur_str.size() != inputs.size() + 2) {
        ERR("size of string values and the number of inputs mismatch! "+name);
        res = 1;
      }
      for (long j = 0; j < cur_str.size() - 2; j++) {
        cur_val <<= 2;
        switch (cur_str[j]) {
          case '0':
            cur_val |= CODE_ZERO;
            break;
          case '1':
            cur_val |= CODE_ONE;
            break;
          case '-':
            cur_val |= CODE_DC;
            break;
          default:
            ERR("unknown value in blif string! "+name);
        }
      }
      // TODO: find redundant values
      coded_values.insert(cur_val);
      //cout << "coded value: " << name << " -> " << hex << cur_val << endl;
      switch (cur_str[cur_str.size() - 1]) {
        case '0':
          if (i > 0) {
            if (result_is_one) {
              ERR("output values mismatch in blif string! "+name);
              res = 1;
            }
          }
          else
            result_is_one = false;
          break;
        case '1':
          if (i > 0) {
            if (!result_is_one) {
              ERR("output values mismatch in blif string! "+name);
              res = 1;
            }
          }
          else
            result_is_one = true;
          break;
        default:
          ERR("unknown output value in blif string! "+name);
          res = 1;
      }
    }
    return res;
  }

  NodeType BlifNode::GetEquivalentType() {
    switch (inputs.size()) {
      case 0:
        if (coded_values.size() != 1)
          return NODE_UNKNOWN;
        else if (result_is_one)
          return NODE_ONE;
        else
          return NODE_ZERO;
        break;
      case 1:
        if (coded_values.size() == 1) {
          switch (*coded_values.begin()) {
            case CODE_DC:
              if (result_is_one)
                return NODE_ONE;
              else
                return NODE_ZERO;
              break;
            case CODE_ONE:
              if (result_is_one)
                return NODE_BUF;
              else
                return NODE_NOT;
              break;
            case CODE_ZERO:
              if (result_is_one)
                return NODE_NOT;
              else
                return NODE_BUF;
              break;
          }
        }
        break;
      case 2:
        switch (coded_values.size()) {
          case 1:
            switch (*coded_values.begin()) {
              case (CODE_ONE | (CODE_ONE << 2)):
                if (result_is_one)
                  return NODE_AND;
                else
                  return NODE_NAND;
                break;
              case (CODE_ZERO | (CODE_ZERO << 2)):
                if (result_is_one)
                  return NODE_NOR;
                else
                  return NODE_OR;
                break;
              case (CODE_ZERO | (CODE_ONE << 2)):
                if (result_is_one)
                  return NODE_AND2_PN;
                else
                  return NODE_NAND2_PN;
                break;
              case (CODE_ONE | (CODE_ZERO << 2)):
                if (result_is_one)
                  return NODE_AND2_NP;
                else
                  return NODE_NAND2_NP;
                break;
            }
            break;
          case 2:
            if ((*coded_values.begin() == (CODE_ONE | (CODE_ZERO << 2))) && (*next(coded_values.begin()) == (CODE_ZERO | (CODE_ONE << 2)))) {
              if (result_is_one)
                return NODE_XOR;
              else
                return NODE_XNOR;
            }
            else if ((*coded_values.begin() == (CODE_ZERO | (CODE_ZERO << 2))) && (*next(coded_values.begin()) == (CODE_ONE | (CODE_ONE << 2)))) {
              if (result_is_one)
                return NODE_XNOR;
              else
                return NODE_XOR;
            }
            else if ((*coded_values.begin() == (CODE_DC | (CODE_ZERO << 2))) && (*next(coded_values.begin()) == (CODE_ZERO | (CODE_DC << 2)))) {
              if (result_is_one)
                return NODE_NAND;
              else
                return NODE_AND;
            }
            else if ((*coded_values.begin() == (CODE_DC | (CODE_ONE << 2))) && (*next(coded_values.begin()) == (CODE_ONE | (CODE_DC << 2)))) {
              if (result_is_one)
                return NODE_OR;
              else
                return NODE_NOR;
            }
            break;
          case 3:
            if ((*coded_values.begin() == (CODE_ONE | (CODE_ZERO << 2))) && (*next(coded_values.begin()) == (CODE_ZERO | (CODE_ONE << 2))) && (*next(next(coded_values.begin())) == (CODE_ONE | (CODE_ONE << 2)))) {
              if (result_is_one)
                return NODE_OR;
              else
                return NODE_NOR;
            }
            else if ((*coded_values.begin() == (CODE_ZERO | (CODE_ZERO << 2))) && (*next(coded_values.begin()) == (CODE_ONE | (CODE_ZERO << 2))) && (*next(next(coded_values.begin())) == (CODE_ZERO | (CODE_ONE << 2)))) {
              if (result_is_one)
                return NODE_NAND;
              else
                return NODE_AND;
            }
            break;
        }
        break;
      default:
        uint64_t all_one_val = 0;
        uint64_t all_zero_val = 0;
        for (long i = 0; i < inputs.size(); i++) {
          all_one_val <<= 2;
          all_zero_val <<= 2;
          all_one_val |= CODE_ONE;
          all_zero_val |= CODE_ZERO;
        }
        if (coded_values.size() == 1) {
          if (*coded_values.begin() == all_one_val) {
            if (result_is_one)
              return NODE_AND;
            else
              return NODE_NAND;
          }
          else if (*coded_values.begin() == all_zero_val) {
            if (result_is_one)
              return NODE_NOR;
            else
              return NODE_OR;
          }
        }
        else if (coded_values.size() == inputs.size()) {
          uint64_t one_dc_val = CODE_ONE;
          uint64_t zero_dc_val = CODE_ZERO;
          uint64_t last_dc_val = CODE_DC << 2 * (inputs.size() - 1);
          for (long k = 0; k < inputs.size() - 1; k++) {
            one_dc_val <<= 2;
            zero_dc_val <<= 2;
            one_dc_val |= CODE_DC;
            zero_dc_val |= CODE_DC;
          }
          set<uint64_t>::iterator iter = coded_values.begin();
          if (*iter == one_dc_val) {
            ++iter;
            while (iter != coded_values.end()) {
              one_dc_val >>= 2;
              one_dc_val |= last_dc_val;
              if (*iter != one_dc_val)
                return type;
              ++iter;
            }
            if (result_is_one)
              return NODE_OR;
            else
              return NODE_NOR;
          }
          else if (*iter == zero_dc_val) {
            ++iter;
            while (iter != coded_values.end()) {
              zero_dc_val >>= 2;
              zero_dc_val |= last_dc_val;
              if (*iter != zero_dc_val)
                return type;
              ++iter;
            }
            if (result_is_one)
              return NODE_NAND;
            else
              return NODE_AND;
          }
        }

    }
    return type; // no change!
  }

}
