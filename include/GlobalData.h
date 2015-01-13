/*
 * GlobalData.h
 *
 *  Created on: Mar 28, 2013
 *      Author: knippsch
 */

#ifndef GLOBALDATA_H_
#define GLOBALDATA_H_

#include <array>
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>

#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SparseCore>

#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/program_options.hpp>

#include "Operators.h"
#include "typedefs.h"

class GlobalData {

private:
  //! A pointer on the class itself
  static GlobalData* instance_;
  //! globally accessible data
  int Lx, Ly, Lz, Lt;
  int dim_row, V_TS, V_for_lime;
  int number_of_eigen_vec;
  int number_of_rnd_vec;
  int number_of_inversions;
  int start_config, end_config, delta_config;
  int verbose;
  std::string path_eigenvectors;
  std::string name_eigenvectors;
  std::string path_perambulators;
  std::string name_perambulators;
  std::string name_lattice;
  std::string path_output;
  std::string path_config;
  std::vector<quark> quarks;
  std::vector<Operator_list> operator_list;
  Correlator_list correlator_list;
  void quark_input_data_handling(
                        const std::vector<std::string> quark_configs);
  void operator_input_data_handling(
                      const std::vector<std::string> operator_list_configs);
  void correlator_input_data_handling(
                         const std::vector<std::string>& correlator_string);

  // new stuff
  void init_from_infile();
  void set_Corr();
  vec_pdg_Corr op_Corr;
  vec_pd_VdaggerV op_VdaggerV;
  vec_pd_rVdaggerVr op_rVdaggerVr;
  //


public:
  static GlobalData* Instance ();

  void read_parameters(int ac, char* av[]);

  inline std::string get_name_lattice() {
    return name_lattice;
  }
  inline std::string get_output_path() {
    return path_output;
  }
  inline std::string get_config_path() {
    return path_config;
  }
  inline int get_Lx () {
    return Lx;
  }
  inline int get_Ly () {
    return Ly;
  }
  inline int get_Lz () {
    return Lz;
  }
  inline int get_Lt () {
    return Lt;
  }
  inline int get_dim_row () {
    return dim_row;
  }
  inline int get_V_TS () {
    return V_TS;
  }
  inline int get_V_for_lime () {
    return V_for_lime;
  }
  inline int get_number_of_inversions () {
    return number_of_inversions;
  }
  inline int get_number_of_rnd_vec () {
    return number_of_rnd_vec;
  }
  inline int get_start_config () {
    return start_config;
  }
  inline int get_end_config () {
    return end_config;
  }
  inline int get_delta_config () {
    return delta_config;
  }
  inline int get_number_of_eigen_vec() {
    return number_of_eigen_vec;
  }
  inline int get_verbose() {
    return verbose;
  }
  inline std::string get_path_eigenvectors() {
    return path_eigenvectors;
  }
  inline std::string get_name_eigenvectors() {
    return name_eigenvectors;
  }
  inline std::string get_path_perambulators() {
    return path_perambulators;
  }
  inline std::string get_name_perambulators() {
    return name_perambulators;
  }
  inline std::vector<quark> get_quarks() {
    return quarks;
  }

  //! All con/de-structors are protected to assure that only one instance exists
  //! at once. DO NOT CHANGE!!
protected:
  GlobalData () {
  }
  GlobalData (const GlobalData& other) {
  }
  virtual ~GlobalData () {
  }

};

#endif /* GLOBALDATA_H_ */
