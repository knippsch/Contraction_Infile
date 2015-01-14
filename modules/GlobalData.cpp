/*
 * GlobalData.cpp
 *
 *  Created on: Mar 28, 2013
 *      Author: knippsch
 */

#include "GlobalData.h"
#include "GlobalData_tests.hpp"

namespace po = boost::program_options;

GlobalData* GlobalData::instance_ = 0;

GlobalData* GlobalData::Instance () {

  if(instance_ == 0) instance_ = new GlobalData;

  return instance_;
}
// *****************************************************************************
/// @brief Convenience function for when a 'store_to' value is being provided
///        to typed_value.
///
/// @param store_to The variable that will hold the parsed value upon notify.
///
/// @return Pointer to a type_value.
template<typename T>
boost::program_options::typed_value<T>* make_value (T* store_to) {
  return boost::program_options::value<T>(store_to);
}
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
bool compare_quantum_numbers_of_pdg(const pdg& in1, const pdg& in2){

  if( (in1.p3 == in2.p3) && 
      (in1.dis3 == in2.dis3) && 
      (in1.gamma == in2.gamma))
    return true;
  else
    return false;

}
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
bool compare_mom_dis_of_pdg(const pdg& in1, const pdg& in2){

  if( (in1.p3 == in2.p3) && 
      (in1.dis3 == in2.dis3))
    return true;
  else
    return false;

}
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
static void copy_quantum_numbers(const pdg& in, std::array<int, 6>& out){
  out[0] = in.dis3[0];
  out[1] = in.dis3[1];
  out[2] = in.dis3[2];
  out[3] = in.p3[0];
  out[4] = in.p3[1]; 
  out[5] = in.p3[2];
}
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
void GlobalData::set_Corr(){

  // first number is the operator id, the next three is the displacement vector
  // and the last three the momentum vector
  std::vector<std::array<int, 6> > rvdaggervr_qu_nb;
  std::vector<std::array<int, 6> > vdaggerv_qu_nb;
  size_t counter_rvdvr = 0;
  size_t counter_vdv = 0;
  for(auto& op : op_Corr){
    std::array<int, 6> write;
    if(op.id != 0){
      copy_quantum_numbers(op, write);
      // ######################################################################
      // check if quantum numbers are already stored in rvdaggervr_qu_nb
      bool is_known_rvdvr = false;
      size_t fast_counter_rvdvr = 0;// this gives the Op id if QN are duplicate
      for(const auto& rvdvr : rvdaggervr_qu_nb){
        if(rvdvr == write){
          is_known_rvdvr = true;
          break;
        }
        fast_counter_rvdvr++;
      }
      if(!is_known_rvdvr){ // setting the unknown quantum numbers
        op.id_rVdaggerVr = counter_rvdvr;
        counter_rvdvr++;
        rvdaggervr_qu_nb.push_back(write);
      }
      else
        op.id_rVdaggerVr = fast_counter_rvdvr;
      // ######################################################################
      // check if quantum numbers are already stored in vdaggerv_qu_nb
      bool is_known_vdv = false;
      size_t fast_counter_vdv = 0;// this gives the Op id if QN are duplicate
      // first check for duplicate quantum numbers
      for(const auto& vdv : vdaggerv_qu_nb){
        if(vdv == write){
          is_known_vdv = true;
          break;
        }
        fast_counter_vdv++;
      }
      if(!is_known_vdv){ // second check for complex conjugate momenta
        fast_counter_vdv = 0;
        for(size_t i = 3; i < 6; i++)
          write[i] *= -1;
        for(const auto& vdv : vdaggerv_qu_nb){
          if(vdv == write){
            is_known_vdv = true;
            break;
          }
          fast_counter_vdv++;
        }
        if(!is_known_vdv){
          op.id_VdaggerV = counter_vdv;
          vdaggerv_qu_nb.push_back(write);
          counter_vdv++;
        }
        else{
          op.flag_VdaggerV = -1;
          op.id_VdaggerV = fast_counter_vdv;
        }
      }
      else{
        op.flag_VdaggerV = 1;
        op.id_VdaggerV = fast_counter_vdv;
      }
    }
    else{ // setting the very first entry
      copy_quantum_numbers(op, write);
      rvdaggervr_qu_nb.push_back(write);
      vdaggerv_qu_nb.push_back(write);
      op.id_VdaggerV = counter_vdv;
      op.id_rVdaggerVr = counter_rvdvr;
      counter_rvdvr++;
      counter_vdv++;
    }
  }

  // setting the lookuptables to be able to reconstruct the quantum numbers
  // when computing VdaggerV and rVdaggerVr
  op_VdaggerV.resize(vdaggerv_qu_nb.size());
  op_rVdaggerVr.resize(rvdaggervr_qu_nb.size());

  size_t index = 0;
  for(auto& op_vdv : op_VdaggerV){
    op_vdv.id = index;
    for(const auto& op : op_Corr){
      if(index == op.id_VdaggerV)
        op_vdv.index = op.id;
    }
    index++;
  }
  index = 0;
  for(auto& op_rvdvr : op_rVdaggerVr){
    op_rvdvr.id = index;
    for(const auto& op : op_Corr){
      if(index == op.id_VdaggerV){
        op_rvdvr.index = op.id;
        if(op.flag_VdaggerV == 1)
          op_rvdvr.adjoint = false;
        else
          op_rvdvr.adjoint = true;
      }
    }
    index++;
  }

}
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
void GlobalData::init_from_infile() {

  // extracting all operators which are used in correlations functions
  std::vector<int> used_operators;
  for(const auto& corr_list : correlator_list)
    used_operators.insert(used_operators.end(), 
                          corr_list.operator_numbers.begin(), 
                          corr_list.operator_numbers.end());
  
  sort(used_operators.begin(), used_operators.end());
  used_operators.erase(std::unique(used_operators.begin(), 
                                   used_operators.end()),
                       used_operators.end());
  // write quantum number in op_Corr
  for(const auto& op_entry : used_operators){
    for(const auto& individual_operator : operator_list[op_entry]){
      pdg write;
      write.gamma = individual_operator.gammas;
      write.dis3 = individual_operator.dil_vec;
      for(auto mom : individual_operator.mom_vec){
        write.p3 = mom;
        op_Corr.push_back(write);
      }
    }
  }
  // doubly counted op_Corr entries are deleted
  auto it = op_Corr.begin();
  while(it != op_Corr.end()) {
    auto it2 = it;
    it2++;
    while(it2 != op_Corr.end()) {
      if(compare_quantum_numbers_of_pdg(*it, *it2))
        op_Corr.erase(it2);
      else
        it2++;
    }
    it++;
  }
  // sorting op_Corr for equal momentum and displacement vectors - makes it
  // easier to run over it with auto loops
  std::vector<pdg> dump_write;
  while(op_Corr.size() != 0){

    it = op_Corr.begin();
    dump_write.push_back(*it);
    op_Corr.erase(it);

    auto it2 = dump_write.end()-1;
    while(it != op_Corr.end()) {
      if(compare_mom_dis_of_pdg(*it, *it2)){
        dump_write.push_back(*it);
        op_Corr.erase(it);
      }
      else
        it++;
    }
  }
  for(auto a : dump_write){
    std::cout << a.gamma;
    for(auto b : a.dis3)
      std::cout << " " << b;
    for(auto b : a.p3)
      std::cout << " " << b;
    std::cout << std::endl;
  }
  op_Corr.swap(dump_write);

  // setting the identification numbers of op_Corr
  size_t counter = 0;
  for(auto& op : op_Corr)
    op.id = counter++;

  // final setting lookuptables for vdaggerv and so on
  set_Corr();

}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void GlobalData::read_parameters (int ac, char* av[]) {

  try{
    std::string input_file;
    std::string output_file;
    // Variables that will store parsed values for quarks.
    std::vector<std::string> quark_configs;
    // Variables that will store parsed values for operators.
    std::vector<std::string> operator_list_configs;
    // Variables that will store parsed values for correlators.
    std::vector<std::string> correlator_list_configs;

    // Declare a group of options that will be allowed only on command line
    po::options_description generic("Command line options");
    generic.add_options()("help,h", "produce help message")("version,v",
        "print version string")("verbose",
        "does additional tests and prints more details")("input,i",
        po::value<std::string>(&input_file)->default_value("LapHs.in"),
        "name of input file.")("output,o",
        po::value<std::string>(&output_file)->default_value("LapHs.out"),
        "name of output file.");

    // Declare a group of options that will be
    // allowed both on command line and in input file
    po::options_description config("Input file options");
    // lattice options
    config.add_options()("output_path",
        po::value<std::string>(&path_output)->
        default_value("../../contractions"),
        "path for output")
        ("config_path",
        po::value<std::string>(&path_config)->default_value("../../configs"),
        "path for configurations")
        ("lattice", po::value<std::string>(&name_lattice)->
        default_value("lattice"),"Codename of the lattice")("Lt", 
        po::value<int>(&Lt)->default_value(0),
        "Lt: temporal lattice extend")("Lx",
        po::value<int>(&Lx)->default_value(0),
        "Lx: lattice extend in x direction")("Ly",
        po::value<int>(&Ly)->default_value(0),
        "Ly: lattice extend in y direction")("Lz",
        po::value<int>(&Lz)->default_value(0),
        "Lz: lattice extend in z direction");
    // eigenvector options
    config.add_options()("number_of_eigen_vec",
        po::value<int>(&number_of_eigen_vec)->default_value(0),
        "Number of eigen vectors")("path_eigenvectors",
        po::value<std::string>(&path_eigenvectors)->default_value("."),
        "directory of eigenvectors")("name_eigenvectors",
        po::value<std::string>(&name_eigenvectors)->default_value(
            "eigenvector"),
        "name of eigenvectors\nThe full name is internally created to:\n"
            "\"name_of_eigenvectors.eigenvector\n."
            "time slice.configuration\"");
    // perambulator options
    config.add_options()("path_perambulators",
        po::value<std::string>(&path_perambulators)->default_value("."),
        "directory of perambulators")("name_perambulators",
        po::value<std::string>(&name_perambulators)->default_value(
            "perambulator"),
        "name of perambulators\nThe full name is internally created to:\n"
            "\"rather long\n."
            "t_sink.configuration\"");
    // quark options
    config.add_options()("quarks.quark", make_value(&quark_configs),
        "quark input, must be of type:\n"
            "quark = \n type:number of rnd. vec.:\n"
            " dil type time:number of dil time:\n"
            " dil type ev:number of dil ev:\n"
            " dil type Dirac:number of dil Dirac");
    // operator list options
    config.add_options()("operator_lists.operator_list", 
        make_value(&operator_list_configs),
        "operator input is rather complicated - see documentation!!");
    // correlator list options
    config.add_options()("correlator_lists.correlator_list", 
        make_value(&correlator_list_configs),
        "correlator input is rather complicated - see documentation!!");
    // configuration options
    config.add_options()("start_config",
        po::value<int>(&start_config)->
        default_value(-1), "First configuration")(
        "end_config", po::value<int>(&end_config)->default_value(0),
        "Last configuration")("delta_config",
        po::value<int>(&delta_config)->default_value(0),
        "Stepsize between two configurations");

    po::options_description cmdline_options;
    cmdline_options.add(generic).add(config);

    po::options_description input_file_options;
    input_file_options.add(config);

    po::options_description visible("Allowed options");
    visible.add(generic).add(config);
    po::positional_options_description p;
    p.add("input-file", -1);

    po::variables_map vm;
    po::store(
        po::command_line_parser(ac, av).options(cmdline_options).
                                        positional(p).run(), vm);
    po::notify(vm);
    // *************************************************************************
    // command line options ****************************************************
    if(vm.count("help")){
      std::cout << visible << "\n";
      exit(0);
    }
    if(vm.count("verbose")){
      verbose = 1;
    }
    else verbose = 0;
    if(vm.count("version")){
      std::cout << "stochastic LapH code, version under construction \n";
      exit(0);
    }
    std::ifstream ifs(input_file.c_str());
    if(!ifs){
      std::cout << "CANNOT open input file: " << input_file << "\n";
      exit(0);
    }
    else{
      po::store(parse_config_file(ifs, input_file_options), vm);
      po::notify(vm);
    }
    ifs.close();

    // *************************************************************************
    // reading input file options **********************************************
    //
    lattice_input_data_handling(path_output, name_lattice, path_config, 
                                                                Lt, Lx, Ly, Lz);
    //
    eigenvec_perambulator_input_data_handling(number_of_eigen_vec,
        path_eigenvectors, name_eigenvectors, path_perambulators,
        name_perambulators);
    //
    quark_input_data_handling(quark_configs);
    //
    operator_input_data_handling(operator_list_configs);
    //
    correlator_input_data_handling(correlator_list_configs);
    //
    config_input_data_handling(start_config, end_config, delta_config);


    // *************************************************************************
    // setting all lookup tables
    init_from_infile();    



    // computing some global variables depending on the input values ***********
    dim_row = Lx * Ly * Lz * 3;

    //needed for config_utils.h
    //4 is number of directions, 3 number of colors and 2 factor
    //for memory requirement of complex numbers
    V_TS = dim_row * 4 * 3 * 2;
    V_for_lime = V_TS * Lt;
  }
  catch(std::exception& e){
    std::cout << e.what() << "\n";
    exit(0);
  }
}
