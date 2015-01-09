
// TODO: Put these things in a name space !!!!


// *****************************************************************************
// A helper function to simplify the main part.
template<class T>
std::ostream& operator<< (std::ostream& os, const std::vector<T>& v) {
  std::copy(v.begin(), v.end(), std::ostream_iterator<T>(os, " "));
  return os;
}
// *****************************************************************************
/// @brief Stream insertion operator for slave.
///
/// @param stream The stream into which quark is being inserted.
/// @param q The quark object.
///
/// @return Reference to the ostream.
std::ostream& operator<< (std::ostream& stream, const quark& quark) {
  return stream << "\tQUARK type: ****  " << quark.type
      << "  ****\n\t number of random vectors: " << quark.number_of_rnd_vec
      << "\n\t dilution scheme in time: " << quark.dilution_T
      << quark.number_of_dilution_T << "\n\t dilution scheme in ev space: "
      << quark.dilution_E << quark.number_of_dilution_E
      << "\n\t dilution scheme in Dirac space: " << quark.dilution_D
      << quark.number_of_dilution_D << "\n";
}
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
// simplifies and cleans read_parameters function
void lattice_input_data_handling (const std::string path_output,
    const std::string name_lattice, const std::string path_config, int Lt,
    int Lx, int Ly, int Lz) {
  try{
    if(Lt < 1){
      std::cout << "\ninput file error:\n" << "\toption \"Lt\""
          << " is mendatory and its value must be an integer greater than 0!"
          << "\n\n";
      exit(0);
    }
    else std::cout << "\n\ttemporal lattice extend .................. " << Lt
        << "\n";
    //
    if(Lx < 1){
      std::cout << "\ninput file error:\n" << "\toption \"Lx\""
          << " is mandatory and its value must be an integer greater than 0!"
          << "\n\n";
      exit(0);
    }
    else std::cout << "\tspatial lattice extend in x direction .... " << Lx
        << "\n";
    //
    if(Ly < 1){
      std::cout << "\ninput file error:\n" << "\toption \"Ly\""
          << " is mandatory and its value must be an integer greater than 0!"
          << "\n\n";
      exit(0);
    }
    else std::cout << "\tspatial lattice extend in y direction .... " << Ly
        << "\n";
    //
    if(Lz < 1){
      std::cout << "\ninput file error:\n" << "\toption \"Lz\""
          << " is mandatory and its value must be an integer greater than 0!"
          << "\n\n\n";
      exit(0);
    }
    else std::cout << "\tspatial lattice extend in z direction .... " << Lz
        << "\n\n";
    std::cout << "\tEnsemble ...................................... " <<
      name_lattice << std::endl;
    std::cout << "\tResults will be saved to path:\n\t\t"
        << path_output << "/" << std::endl;
    std::cout << "\tConfigurations will be read from:\n\t\t"
        << path_config << "/" << std::endl;
  }
  catch(std::exception& e){
    std::cout << e.what() << "\n";
    exit(0);
  }
}
// *****************************************************************************
// simplifies and cleans read_parameters function
void eigenvec_perambulator_input_data_handling (
    const int number_of_eigen_vec, const std::string path_eigenvectors,
    const std::string name_eigenvectors, const std::string path_perambulators,
    const std::string name_perambulators) {

  try{
    if(number_of_eigen_vec < 1){
      std::cout << "\ninput file error:\n" << "\toption \"number_of_eigen_vec\""
          << " is mandatory and its value must be an integer greater than 0!"
          << "\n\n";
      exit(0);
    }
    else{
      std::cout << "\tnumber of eigen vectors .................. "
          << number_of_eigen_vec << "\n";
    }
    std::cout << "\tEigenvectors will be read from files:\n\t\t"
        << path_eigenvectors << "/" << name_eigenvectors
        << "\".eigenvector.t.config\"\n";
    std::cout << "\tPerambulators will be read from files:\n\t\t"
        << path_perambulators << "/" << name_perambulators
        << "\".rnd_vec.scheme.t_sink.config\"\n\n";
  }
  catch(std::exception& e){
    std::cout << e.what() << "\n";
    exit(0);
  }
}
// *****************************************************************************
// simplifies and cleans read_parameters function
void config_input_data_handling (const int start_config,
    const int end_config, const int delta_config) {

  try{
    if(start_config < 0){
      std::cout << "\ninput file error:\n" << "\toption \"start config\""
          << " is mandatory and its value must be an integer greater or equal 0!"
          << "\n\n";
      exit(0);
    }
    else if(end_config < 1 || end_config < start_config){
      std::cout << "\ninput file error:\n" << "\toption \"end_config\""
          << " is mandatory, its value must be an integer greater than 0,"
          << " and it must be larger than start config!" << "\n\n";
      exit(0);
    }
    else if(delta_config < 1){
      std::cout << "\ninput file error:\n" << "\toption \"delta_config\""
          << " is mandatory and its value must be an integer greater than 0!"
          << "\n\n";
      exit(0);
    }
    else std::cout << "\tprocessing configurations " << start_config << " to "
        << end_config << " in steps of " << delta_config << "\n\n";
  }
  catch(std::exception& e){
    std::cout << e.what() << "\n";
    exit(0);
  }
}
// *****************************************************************************
// simplifies and cleans read_parameters function
void quark_check (quark quarks) {

  try{
    if(quarks.type != "u" && quarks.type != "d" && quarks.type != "s"
        && quarks.type != "c"){
      std::cout << "quarks.quark.type must be u, d, s or c" << std::endl;
      exit(0);
    }
    else if(quarks.number_of_rnd_vec < 1){
      std::cout << "quarks.quark.number_of_rnd_vec must be greater than 0"
          << std::endl;
      exit(0);
    }
    else if(quarks.dilution_T != "TI" && quarks.dilution_T != "TB"){
      std::cout << "quarks.quark.dilutione_T must be TI or TB" << std::endl;
      exit(0);
    }
    else if(quarks.number_of_dilution_T < 1){
      std::cout << "quarks.quark.number_of_dilution_T must be greater than 0 "
          "and smaller than the temporal extend" << std::endl;
      exit(0);
    }
    else if(quarks.dilution_E != "EI" && quarks.dilution_E != "EB"){
      std::cout << "quarks.quark.dilutione_E must be EI or EB" << std::endl;
      exit(0);
    }
    else if(quarks.number_of_dilution_E < 1){
      std::cout << "quarks.quark.number_of_dilution_E must be greater than 0 "
          "and smaller than number of eigen vectors" << std::endl;
      exit(0);
    }
    else if(quarks.dilution_D != "DI" && quarks.dilution_D != "DI"){
      std::cout << "quarks.quark.dilutione_D must be DI or DB" << std::endl;
      exit(0);
    }
    else if(quarks.number_of_dilution_D < 1 || quarks.number_of_dilution_D > 4){
      std::cout << "quarks.quark.number_of_dilution_D must be greater than 0 "
          "and smaller than 5" << std::endl;
      exit(0);
    }
    else std::cout << quarks << std::endl;
  }
  catch(std::exception& e){
    std::cout << e.what() << "\n";
    exit(0);
  }

}
