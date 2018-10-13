#include "global.h"
#include <chrono>

using namespace std;

#include "circuit.h"

using namespace nodecircuit;

#include "verilog_gate_parser.h"
#include "blif_parser.h"

// these are duplicate definitions! in two cpp files..
const char default_spec_name[] = "new-spec-final.blif";
const char default_impl_name[] = "new-impl-final.blif";

const char default_exsim_out_name[] = "exsim_res.out";

extern int PrepareDebugFiles(string filename, vector<string> &targets, string out_node_name, int num_min_change);
extern int DoDebug(string spec_filename, string impl_filename, int num_init_patterns, string pattern_filename);
extern int DoExSim(string spec_filename, string output_filename);
extern int DoSimEq(string spec_filename, string impl_filename);
extern int RndSimEq(string spec_filename, string impl_filename, int num_patterns);
extern int DoSimEqMod(string spec_filename, string target_node_name);
extern int CreateNodeMiter(string impl_filename, string node_name);
extern int SimActivity(string spec_filename);
extern int SatIncActivity(string spec_filename, int high_act_level, int low_act_level, int num_cycles);
extern int IncActivity(string spec_filename, string allregions_filename, string highregion_filename, int high_act_level, int num_cycles, bool rndsim);
extern int DecActivity(string spec_filename, string allregions_filename, string lowregion_filename, int low_act_level, int num_cycles, bool rndsim);
extern int SandActivity(string spec_filename, string allregions_filename, string sandregion_filename, int low_act_level, int high_act_level, int num_cycles_low, int num_cycles_high);


extern int TestSimEq(string spec_filename, int method);
extern int DoSimEqLim(string spec_filename, int eqmode, int amount, int submode);
extern int DoSimEqConnect(string spec_filename, int mode, int amount);
extern int DoSimEqModTwo(string spec_filename_1, string spec_filename_2, int mode, int amount);


extern int MyTest(string filename);

const char sim_debug_str[] = "simdebug";
const char add_lut_str[] = "addlut";
const char minterm_change_str[] = "mintermchange";
const char ex_sim_str[] = "exsim";
const char sim_eq_str[] = "simeq";
const char rnd_sim_eq_str[] = "rndsimeq";
const char sim_eq_mod_str[] = "simeqmod";
const char node_miter_str[] = "nodemiter";
const char sim_act_str[] = "simact";
const char sat_act_inc_str[] = "satactinc";
const char inc_act_str[] = "incact";
const char dec_act_str[] = "decact";
const char sand_act_str[] = "sandact";


const char sim_eq_lim_str[] = "simeqlim";
const char sim_eq_con_str[] = "simeqcon";
const char test_sim_eq[] = "testeq";
const char sim_eq_mod_two_str[] = "simeqmodtwo";


inline void msg_sim_debug_cmd () {
  MSG(string("  ")+sim_debug_str+" [NumInitialPatterns] or");
  MSG(string("  ")+sim_debug_str+" [InputPatternFile]");
  MSG("      ->  Debugging with the previously prepared spec/impl by simulating spec");
  MSG("      ->  Optionally number of initial random patterns can be specified");
  MSG("      ->    (default=1) or");
  MSG("      ->  Optionally file name with initial patterns can be specified");
}
inline void msg_add_lut_cmd () {
  MSG(string("  ")+add_lut_str+" SpecFile Node1 [ Node2 ... ]");
  MSG("      ->  Inserting LUT for the given Nodes");
}
inline void msg_minterm_change_cmd () {
  MSG(string("  ")+minterm_change_str+" SpecFile PO NumMinterms Node1 [ Node2 ... ]");
  MSG("      ->  Changing NumMinterms minterms of PO randomly ");
  MSG("      ->  as well as inserting LUT for the given Nodes");
}
inline void msg_ex_sim_cmd () {
  MSG(string("  ")+ex_sim_str+" SpecFile");
  MSG("      ->  Exhaustive simulation of SpecFile");
}
inline void msg_sim_eq_cmd () {
  MSG(string("  ")+sim_eq_str+" SpecFile ImplFile");
  MSG("      ->  Equivalency checking between SpecFile and ImplFile by simulation");
}
inline void msg_rnd_sim_eq_cmd () {
  MSG(string("  ")+rnd_sim_eq_str+" SpecFile ImplFile NumPatterns");
  MSG("      ->  Equivalency checking between SpecFile and ImplFile");
  MSG("      ->  by random simulating ~NumPatterns");
}
inline void msg_sim_eq_mod_cmd () {
  MSG(string("  ")+sim_eq_mod_str+" SpecFile Node");
  MSG("      ->  Equivalency checking between SpecFile and the Node modified version by simulation");
}
inline void msg_node_miter_cmd () {
  MSG(string("  ")+node_miter_str+" ImplFile Node");
  MSG("      ->  Create a miter by removing Node and assigning it to both 0 and 1");
}
inline void msg_sim_act_cmd () {
  MSG(string("  ")+sim_act_str+" SpecFile");
  MSG("      ->  Measure signal activity by random simulation of SpecFile");
}
inline void msg_sat_act_inc_cmd () {
  MSG(string("  ")+sat_act_inc_str+" SpecFile HighSA LowSA [NumCycles]");
  MSG("      ->  Increase signal activity of SpecFile to target HighSA and LowSA (1..99)");
  MSG("      ->  for NumCycles (default=1) by SAT starting with a random input");
}
inline void msg_inc_act_cmd () {
  MSG(string("  ")+inc_act_str+" SpecFile AllRegions HighRegion HighSA [NumCycles]");
  MSG("      ->  Increase signal activity of SpecFile to target HighSA(1..99)");
  MSG("      ->  for NumCycles (default=1) by SAT starting with a random input");
  MSG("      ->  SpecFile is a Verilog file synthesized with specific library");
  MSG("      ->  AllRegions have the name of the files having signals names in that region");
  MSG("      ->  Highregion have the name of the signals names in the high region");
}
inline void msg_dec_act_cmd () {
  MSG(string("  ")+dec_act_str+" SpecFile AllRegions LowRegion LowSA [NumCycles]");
  MSG("      ->  Decrease signal activity of SpecFile to target LowSA(1..99)");
  MSG("      ->  for NumCycles (default=1) by SAT starting with a random input");
  MSG("      ->  SpecFile is a Verilog file synthesized with specific library");
  MSG("      ->  AllRegions have the name of the files having signals names in that region");
  MSG("      ->  Lowregion have the name of the signals names in the low region");
}
inline void msg_sand_act_cmd () {
  MSG(string("  ")+sand_act_str+" SpecFile AllRegions SandwitchedRegion LowSA HighSA NumCyclesLow NumCyclesHigh");
  MSG("      ->  Decrease signal activity of SpecFile to target LowSA(1..99) for NumCyclesLow");
  MSG("      ->  Then increase signal activity to target HighSA(1..99) for NumCyclesHigh");
  MSG("      ->  Then again decrease signal activity to target LowSA(1..99) for NumCyclesLow");
  MSG("      ->  SpecFile is a Verilog file synthesized with specific library");
  MSG("      ->  AllRegions have the name of the files having signals names in that region");
  MSG("      ->  SandwitchedRegion have the name of the signals names in the low region");
}

inline void msg_sim_eq_lim_cmd () {
  MSG(string("  ")+sim_eq_lim_str+" SpecFile "+ " Run Mode " + "Case Number");
  MSG("      ->  Equivalency checking of SpecFile under Gate Type Change By Different Run Mode");
}
inline void msg_sim_eq_con_cmd () {
  MSG(string("  ")+sim_eq_con_str+" SpecFile "+" Case Number");
  MSG("      ->  Equivalency checking between SpecFile and the Connection modified version by simulation through (parameters) cases");
}
inline void sim_eq_mod_two_cmd () {
  MSG(string("") + sim_eq_mod_two_str + "SpecFile 1 " + "SpecFile 2" + "Run Mode " + "Case Number ");
  MSG("       -> Equvalency checking to SpecFile 1 and SpecFile 2 in terms of gate type change");
}



// TODO: check the circuit is combinational for the commands that expect the circuit to be combinational!

int main(int argc, char *argv[]) {

  if (argc == 1) {
    MSG("Use one of the following commands: ");
    msg_sim_debug_cmd();
    msg_add_lut_cmd();
    msg_minterm_change_cmd();
    msg_ex_sim_cmd();
    msg_sim_eq_cmd();
    msg_rnd_sim_eq_cmd();
    msg_sim_eq_mod_cmd();

    msg_node_miter_cmd();
    return 0;
  }

  string cmd_str(argv[1]);

  if (cmd_str == sim_debug_str) {
    if (argc > 3) {
      MSG("Too many parameters! Use as following:");
      msg_sim_debug_cmd();
      return 0;
    }
    int num_init_patterns = 1;
    string pat_str;
    if (argc == 3) {
      pat_str = argv[2];
      if (pat_str[0] < '0' || pat_str[0] > '9')
        // input parameter is not a number!
        num_init_patterns = -1;
      else {
        size_t p;
        num_init_patterns = stoi(pat_str, &p);
        if (p != pat_str.size())
          // input parameter is not a number!
          num_init_patterns = -1;
      }
    }
    MSG("Debugging with prepared sepc/impl files...");
    DoDebug(default_spec_name, default_impl_name, num_init_patterns, pat_str);
  }
  else if (cmd_str == add_lut_str) {
    if (argc < 4) {
      MSG("Too few parameters! Use as following:");
      msg_add_lut_cmd();
      return 0;
    }
    MSG("Preparing spec/impl files for debug by adding LUTs...");
    vector<string> targets;
    string out_node_min_change;
    int min_change_num = 0;
    for (int i = 3; i < argc; i++) {
      targets.push_back(argv[i]);
    }
    PrepareDebugFiles(argv[2], targets, out_node_min_change, min_change_num);
    MSG("Run again with \"simdebug\" command to do debug!");
  }
  else if (cmd_str == minterm_change_str) {
    if (argc < 5) {
      MSG("Too few parameters! Use as following:");
      msg_minterm_change_cmd();
      return 0;
    }
    MSG("Preparing spec/impl files for debug by changing minterms and adding LUTs...");
    vector<string> targets;
    string out_node_min_change = argv[3];
    int min_change_num = stoi(string(argv[4]));
    for (int i = 5; i < argc; i++) {
      targets.push_back(argv[i]);
    }
    PrepareDebugFiles(argv[2], targets, out_node_min_change, min_change_num);
    MSG("Run again with \"simdebug\" command to do debug!");
  }
  else if (cmd_str == ex_sim_str) {
    if (argc != 3) {
      MSG("Incorrect number of parameters! Use as following:");
      msg_ex_sim_cmd();
      return 0;
    }
    MSG("Exhaustive simulating spec file...");
    DoExSim(argv[2], default_exsim_out_name);
  }
  else if (cmd_str == sim_eq_str) {
    if (argc != 4) {
      MSG("Incorrect number of parameters! Use as following:");
      msg_sim_eq_cmd();
      return 0;
    }
    MSG("Equivalency checking by simulation...");
    DoSimEq(argv[2], argv[3]);
  }
  else if (cmd_str == rnd_sim_eq_str) {
    if (argc != 5) {
      MSG("Incorrect number of parameters! Use as following:");
      msg_rnd_sim_eq_cmd();
      return 0;
    }
    MSG("Equivalency checking by random simulation...");
    RndSimEq(argv[2], argv[3], stoi(string(argv[4])));
  }
  else if (cmd_str == sim_eq_mod_str) {
    if (argc != 4) {
      MSG("Incorrect number of parameters! Use as following:");
      msg_sim_eq_mod_cmd();
      return 0;
    }
    MSG("Modifying and equivalency checking by simulation...");
    DoSimEqMod(argv[2], argv[3]);
  }
  else if (cmd_str == node_miter_str) {
    if (argc != 4) {
      MSG("Incorrect number of parameters! Use as following:");
      msg_node_miter_cmd();
      return 0;
    }
    MSG("Creating miter from node...");
    CreateNodeMiter(argv[2], argv[3]);
  }
  else if (cmd_str == sim_act_str) {
    if (argc != 3) {
      MSG("Incorrect number of parameters! Use as following:");
      msg_sim_act_cmd();
      return 0;
    }
    MSG("Simulating and measuring activities...");
    SimActivity(argv[2]);
  }
  else if (cmd_str == sat_act_inc_str) {
    if (argc != 5 && argc != 6) {
      MSG("Incorrect number of parameters! Use as following:");
      msg_sat_act_inc_cmd();
      return 0;
    }
    MSG("Trying to increasing signal activity by SAT...");
    SatIncActivity(argv[2], stoi(string(argv[3])), stoi(string(argv[4])), argc == 6 ? stoi(string(argv[5])) : 1 );
  }
  else if (cmd_str == inc_act_str) {
    if (argc != 6 && argc != 7) {
      MSG("Incorrect number of parameters! Use as following:");
      msg_inc_act_cmd();
      return 0;
    }
    MSG("Trying to increase signal activity of region(s) by SAT...");
    IncActivity(argv[2], argv[3], argv[4], stoi(string(argv[5])), argc == 7 ? stoi(string(argv[6])) : 1, false);
  }
  else if (cmd_str == dec_act_str) {
    if (argc != 6 && argc != 7) {
      MSG("Incorrect number of parameters! Use as following:");
      msg_dec_act_cmd();
      return 0;
    }
    MSG("Trying to decrease signal activity of region(s) by SAT...");
    DecActivity(argv[2], argv[3], argv[4], stoi(string(argv[5])), argc == 7 ? stoi(string(argv[6])) : 1, false);
  }
  else if (cmd_str == sand_act_str) {
    if (argc != 9) {
      MSG("Incorrect number of parameters! Use as following:");
      msg_sand_act_cmd();
      return 0;
    }
    MSG("Trying to decrease/increase signal activity of sandwitched region(s) by SAT/Sim...");
    SandActivity(argv[2], argv[3], argv[4], stoi(string(argv[5])), stoi(string(argv[6])), stoi(string(argv[7])), stoi(string(argv[8])));
  }
  else if (cmd_str == "mytest") {
    MSG("mytest command is JUST for my test!");
    MyTest(argv[2]);
  }
  else if (cmd_str == test_sim_eq) {
    if (argc != 4) {
      MSG("Incorrect number of parameters! Use as following:");
      msg_sat_act_inc_cmd();
      return 0;
    }
    TestSimEq(argv[2], stoi(string(argv[3])));
  }
  else if (cmd_str == sim_eq_lim_str) {
    if (argc != 5 && argc != 6) {
      MSG("Incorrect number of parameters! Use as following:");
      msg_sim_eq_lim_cmd();
      return 0;
    }

    int mode = stoi(string(argv[3]));
    int amount = stoi(string(argv[4]));

    int submode = 0;
    if (argc == 6) {
      submode = stoi(string(argv[5]));
    }

    if(amount!=-1 && amount<0){
      MSG("You must select more than 1 gates or input -1 to select all gates ");
    }

    if(amount==-1)  MSG("You have chose all gates");
    else MSG("You select "+string(argv[4])+" gates randomly");
    if(mode==0)  MSG("This mode is to do exhaustive equivalency checking of selected gates");
    if(mode==1) MSG("This mode is to do exhaustive equivalency checking after checking its wrong cases after random simulation");
    if(mode==2) MSG("This mode is to do random simulation of selected gate by 10000 times");
    if(mode==3) MSG("This mode is to do random simulation by difference times(1000/5000/10000/20000/50000)");

    DoSimEqLim(argv[2], mode, amount, submode);
  }
  else if (cmd_str == sim_eq_con_str) {
    if (argc != 5) {
      MSG("Incorrect number of parameters! Use as following:");
      msg_sim_eq_con_cmd();
      return 0;
    }
    int amount = stoi(string(argv[4]));
    int mode = stoi(string(argv[3]));
    if(amount==-1)  MSG("You have chose all gates");
    else MSG("You select "+string(argv[4])+" gates randomly");
    if(mode==0) MSG("This mode is to do exhaustive simulation of selected gates from connection change");
    if(mode==1) MSG("This mode is to do exhaustive equivalency checking after checking its wrong cases after random simulation");
    if(mode==2) MSG("This mode is to do random simulation of selected gate by 10000 times");
    if(mode==3) MSG("This mode is to do random simulation by difference times(1000/5000/10000/20000/50000)");
    DoSimEqConnect(argv[2], mode , amount);
  }
  else if (cmd_str == sim_eq_mod_two_str) {
    if (argc != 6) {
      MSG("Incorrect number of parameters! Use as following:");
      sim_eq_mod_two_cmd();
      return 0;
    }
    int amount = stoi(string(argv[5]));
    int mode = stoi(string(argv[4]));
    if(amount==-1)  MSG("You have chose all gates");
    else MSG("You select "+string(argv[5])+" gates randomly");
    DoSimEqModTwo(argv[2], argv[3], mode , amount);
  }
  else {
    MSG("Use one of the following commands: ");
    msg_sim_debug_cmd();
    msg_add_lut_cmd();
    msg_minterm_change_cmd();
    msg_ex_sim_cmd();
    msg_sim_eq_cmd();
    msg_rnd_sim_eq_cmd();
    msg_sim_eq_mod_cmd();

    msg_node_miter_cmd();
    return 0;
  }

  return 0;
}
