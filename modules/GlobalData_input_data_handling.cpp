#include "GlobalData.h"


// *****************************************************************************
/// @brief Makes a quark object from a string
static quark make_quark (const std::string& quark_string) {
  // Tokenize the string on the ":" delimiter.
  std::vector<std::string> tokens;
  boost::split(tokens, quark_string, boost::is_any_of(":"));

  // If the split did not result in exactly 8 tokens, then the value
  // is formatted wrong.
  if(8 != tokens.size()){
    using boost::program_options::validation_error;
    throw validation_error(validation_error::invalid_option_value,
        "quarks.quark", quark_string);
  }

  // Create a quark from the token values.
  return quark(tokens[0], boost::lexical_cast<int>(tokens[1]), tokens[2],
      boost::lexical_cast<int>(tokens[3]), tokens[4],
      boost::lexical_cast<int>(tokens[5]), tokens[6],
      boost::lexical_cast<int>(tokens[7]));
}
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
static std::array<int, 3> create_3darray_from_string(std::string in) { 

  std::array<int, 3> out;
  std::vector<std::string> tokens;
  // erasing the brakets at the beginning and the end
  in.erase(0,2);
  in.erase(in.end()-1);

  boost::split(tokens, in, boost::is_any_of(","));

  return {{boost::lexical_cast<int>(tokens[0]),
          boost::lexical_cast<int>(tokens[1]),
          boost::lexical_cast<int>(tokens[2]) }};

}
// *****************************************************************************
static void create_all_momentum_combinations(const std::vector<int>& in, 
                                        std::vector<std::array<int, 3> >& out) {
  // creating all momentum combinations possible and needed
  int max_p = sqrt(*std::max_element(in.begin(), in.end()));
  std::vector<std::array<int, 3> > all_p;
  for(int p1 = -max_p; p1 < max_p+1; p1++)
    for(int p2 = -max_p; p2 < max_p+1; p2++)
      for(int p3 = -max_p; p3 < max_p+1; p3++)
        all_p.push_back({{p1, p2, p3}});
  // copying wanted combinations into out array
  for(const auto& p : in)
    for(const auto& all : all_p)
      if(p == all[0]*all[0] + all[1]*all[1] + all[2]*all[2])
        out.push_back(all);

}
// *****************************************************************************
static void create_mom_array_from_string(std::string in, 
                                        std::vector<std::array<int, 3> >& out) {
  // erase the p (first entry)
  in.erase(0,1);
  std::vector<std::string> tokens;
  boost::split(tokens, in, boost::is_any_of(","));
  std::vector<int> p;
  for(const auto& t : tokens)
    p.push_back(boost::lexical_cast<int>(t));

  create_all_momentum_combinations(p, out);

}
// *****************************************************************************
/// @brief Makes an operator list object from a string
Operator_list make_operator_list(const std::string& operator_string) {

  Operator_list op_list; // return object

  // Two steps are necessary: 1. Getting all operators in one list which are 
  //                             separated by ";"
  //                          2. Separating the individual operators into its
  //                             smaller bits, which are separated by "."
  // Tokenize the string on the ";" delimiter -> Individual operators
  std::vector<std::string> operator_tokens;
  boost::split(operator_tokens, operator_string, boost::is_any_of(":"));

  // running over opeator tokens and split them further (Step 2):
  for (const auto& op_t : operator_tokens){
    std::vector<std::string> tokens;
    boost::split(tokens, op_t, boost::is_any_of("."));
    std::vector<int> gammas;
    std::array<int, 3> dil_vec;
    std::vector<std::array<int, 3> > mom_vec;
    for (auto str : tokens){
      // getting the gamma structure
      if(str.compare(0,1,"g") == 0)
        gammas.push_back(boost::lexical_cast<int>(str.erase(0,1)));
      // getting the displacement indices
      else if (str.compare(0,1,"d") == 0) {
        if(str.compare(1,1,"0") == 0)
          dil_vec = {0, 0, 0};
        else if (str.compare(1,1,"(") == 0)
          dil_vec = create_3darray_from_string(str);
        else {
         std::cout << "Something wrong with the displacement in the operator" \
                      " definition" << std::endl;
         exit(0);
        }
      }
      // getting the momenta
      else if (str.compare(0,1,"p") == 0) {
        if(str.compare(1,1,"(") == 0)
          mom_vec.push_back(create_3darray_from_string(str));
        else 
          create_mom_array_from_string(str, mom_vec);
      }
      // catching wrong entries
      else {
        std::cout << "there is something wrong with the operators" << std::endl;
        exit(0);
      }
    }
    op_list.push_back(Operators(gammas, dil_vec, mom_vec));
  }
  return op_list;
}
// *****************************************************************************
/// @brief Makes an operator list object from a string
void GlobalData::correlator_input_data_handling (
                             const std::vector<std::string>& correlator_string){

  for(auto str : correlator_string){
  
    std::vector<std::string> correlator_tokens;
    boost::split(correlator_tokens, str, boost::is_any_of(":"));
  
    std::string type;
    std::vector<int> quark_number;
    std::vector<int> operator_number;
    std::string GEVP;
    std::vector<int> tot_mom;
    for (auto corr_t : correlator_tokens){
      // getting the type name
      if (corr_t.compare(0,1,"C") == 0)
        type = corr_t;
      // getting quark numbers
      else if (corr_t.compare(0,1,"Q") == 0) 
        quark_number.push_back(boost::lexical_cast<int>(corr_t.erase(0,1)));
      // getting operator numbers
      else if (corr_t.compare(0,2,"Op") == 0)
        operator_number.push_back(boost::lexical_cast<int>(corr_t.erase(0,2)));
      // getting the GEVP type
      else if (corr_t.compare(0,1,"G") == 0)
        GEVP = corr_t;
      // getting total momenta for moving frames
      else if (corr_t.compare(0,1,"P") == 0) {
        corr_t.erase(0,1);
        std::vector<std::string> tokens;
        boost::split(tokens, corr_t, boost::is_any_of(","));
        for(auto t : tokens)
          tot_mom.push_back(boost::lexical_cast<int>(t));
      }
      // catching wrong entries
      else {
        std::cout << "there is something wrong with the correlators" 
                  << std::endl;
        exit(0);
      }
    }
    correlator_list.push_back(Correlators
                          (type, quark_number, operator_number, GEVP, tot_mom));
  }
  // TODO: write check for correctness of input data
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void GlobalData::quark_input_data_handling (
    const std::vector<std::string> quark_configs) {
  try{
    // Transform each configured quark into a quark via make_quark, 
    // inserting each object into the quark vector.
    std::transform(quark_configs.begin(), quark_configs.end(),
        std::back_inserter(quarks), make_quark);
    // checking the contents for correctness
    std::for_each(quarks.begin(), quarks.end(), quark_check);
  }
  catch(std::exception& e){
    std::cout << e.what() << "\n";
    exit(0);
  }
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void GlobalData::operator_input_data_handling (
    const std::vector<std::string> operator_list_configs) {
  try{
    // Transform each configured quark into a quark via make_quark,
    // inserting each object into the quark vector.
    //std::transform(operator_list_configs.begin(), operator_list_configs.end(),
    //    std::back_inserter(operator_list), make_operator_list);
    for(auto op_list : operator_list_configs)
      operator_list.push_back(make_operator_list(op_list));
    // TODO write a check for correctness of input    
  }
  catch(std::exception& e){
    std::cout << "operator_input_data_handling: " << e.what() << "\n";
    exit(0);
  }
}
