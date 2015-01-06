/*
 * operators.h
 *
 *  Created on: Jan, 2015
 *      Author: knippsch
 */

#ifndef OPERATORS_H_
#define OPERATORS_H_


// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
/// @brief operator type that contains all operator informations
struct Operators {

private:
  std::vector<int> gammas;
  std::array<int, 3> dil_vec;
  std::vector<std::array<int, 3> > mom_vec;

public:
	/// @brief Constructor.
	Operators (std::vector<int> gammas, std::array<int, 3> dil_vec, 
             std::vector<std::array<int, 3> > mom_vec) :
                           gammas(gammas), dil_vec(dil_vec), mom_vec(mom_vec) {}

};
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
/// @brief correlator type that contains all correlator informations
struct Correlators {

private:
  std::string type;
  std::vector<quark> quarks;
  std::vector<Operators> operators;

public:
	/// @brief Constructor.
	Correlators (std::string type, std::vector<quark> quarks, 
             std::vector<Operators> operators) :
                           type(type), quarks(quarks), operators(operators) {}

};
// *****************************************************************************
// *****************************************************************************
// *****************************************************************************
typedef std::vector<Operators> Operator_list;
typedef std::vector<Correlators> Correlator_list;

#endif /* OPERATORS_H_ */
