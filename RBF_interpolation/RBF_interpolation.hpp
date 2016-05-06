#include <Eigen/Dense>
#include <iostream>
using namespace Eigen;
using namespace std;


MatrixXd rbfcreate(MatrixXd x, MatrixXd y, MatrixXd RBFFunction(MatrixXd r, double const_num),
	double RBFConstant, double RBFSmooth);
MatrixXd rbfphi_linear(MatrixXd r, double const_num);
MatrixXd rbfphi_cubic(MatrixXd r, double const_num);
MatrixXd rbfphi_gaussian(MatrixXd r, double const_num);
MatrixXd rbfphi_multiquadrics(MatrixXd r, double const_num);
MatrixXd rbfphi_thinplate(MatrixXd r, double const_num);
MatrixXd rbfAssemble(MatrixXd x, MatrixXd RBFFunction(MatrixXd r, double const_num),
	double RBFConstant, double RBFSmooth);
MatrixXd rbfinterp(MatrixXd nodes, MatrixXd rbfcoeff, MatrixXd x,
	MatrixXd RBFFunction(MatrixXd r, double const_num),
	double RBFConstant);
double rbfcheck(MatrixXd x, MatrixXd y, MatrixXd rbfcoeff,
	MatrixXd RBFFunction(MatrixXd r, double const_num),
	double RBFConstant);