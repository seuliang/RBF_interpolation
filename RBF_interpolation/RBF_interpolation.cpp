#include <RBF_interpolation.hpp>

MatrixXd rbfcreate(MatrixXd x,							//采样点坐标,每一行是一个数据
	MatrixXd y,											//采样点值,函数值是一列多行
	MatrixXd RBFFunction(MatrixXd r, double const_num),	//插值函数
	double RBFConstant,									//插值函数参数
	double RBFSmooth = 0)								//平滑参数 0是不平滑
	//求RBF插值参数
{
	int nXDim = x.cols();
	int nYDim = y.cols();
	int nXCount = x.rows();
	int nYCount = y.rows();
	if (nXCount != nYCount)
		cerr << "x and y should have the same number of rows" << endl;
	if (nYDim != 1)
		cerr << "y should be n by 1 vector" << endl;
	MatrixXd A;
	A = rbfAssemble(x, RBFFunction, RBFConstant, RBFSmooth);
	MatrixXd b = MatrixXd::Zero(y.rows() + nXDim + 1, 1);
	b.topRows(y.rows()) = y;
	MatrixXd rbfcoeff;
	rbfcoeff = A.lu().solve(b);
	return rbfcoeff;
}

MatrixXd rbfAssemble(MatrixXd x,						//采样点坐标
	MatrixXd RBFFunction(MatrixXd r, double const_num),	//插值函数
	double RBFConstant,									//插值函数参数
	double RBFSmooth)									//平滑参数
	//构造系数矩阵
{
	int dim = x.cols();//维度数
	int n = x.rows();//采样点数
	MatrixXd r = MatrixXd::Zero(n, n);
	MatrixXd temp_A;
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < i; j++)
		{
			r(i, j) = (x.row(i) - x.row(j)).norm();
			r(j, i) = r(i, j);
		}
	}
	temp_A = RBFFunction(r, RBFConstant);
	//平滑:使得采样点不再精确等于采样值
	for (int i = 0; i < n; i++)
	{
		temp_A(i, i) = temp_A(i, i) - RBFSmooth;
	}
	MatrixXd P(x.rows(), x.cols() + 1);
	P.leftCols(1) = MatrixXd::Ones(n, 1);
	P.rightCols(x.cols()) = x;
	MatrixXd A = MatrixXd::Zero(temp_A.rows() + P.cols(), temp_A.cols() + P.cols());
	A.topLeftCorner(temp_A.rows(), temp_A.cols()) = temp_A;
	A.topRightCorner(P.rows(), P.cols()) = P;
	A.bottomLeftCorner(P.cols(), P.rows()) = P.transpose();
	return A;
}

MatrixXd rbfphi_linear(MatrixXd r, double const_num)
//Linear
{
	MatrixXd u(r.rows(), r.cols());
	u = r;
	return u;
}
MatrixXd rbfphi_cubic(MatrixXd r, double const_num)
//Cubic
{
	MatrixXd u(r.rows(), r.cols());
	u = r.array().pow(3);
	return u;
}
MatrixXd rbfphi_gaussian(MatrixXd r, double const_num)
//Gaussian
{
	MatrixXd u(r.rows(), r.cols());
	u = (-0.5*r.array().square() / (const_num*const_num)).array().exp();
	return u;
}
MatrixXd rbfphi_multiquadrics(MatrixXd r, double const_num)
//Multiquadrics
{
	MatrixXd u(r.rows(), r.cols());
	u = (1.0 + r.array().square() / (const_num*const_num)).array().sqrt();
	return u;
}
MatrixXd rbfphi_thinplate(MatrixXd r, double const_num)
//Thinplate
{
	MatrixXd u(r.rows(), r.cols());
	u = (r.array().square()).cwiseProduct((r.array() + 1).log());
	return u;
}

MatrixXd rbfinterp(MatrixXd nodes,						//采样点坐标
	MatrixXd rbfcoeff,									//RBF插值参数
	MatrixXd x,											//待插值点坐标
	MatrixXd RBFFunction(MatrixXd r, double const_num),	//插值函数
	double RBFConstant)									//插值函数参数
//计算待插值点的值
{
	int dim = nodes.cols();
	int n = nodes.rows();
	int dimPoints = x.cols();
	int nPoints = x.rows();
	if (dim != dimPoints)
		cerr << "x should have the same number of rows as an array used to create RBF interpolation" << endl;
	MatrixXd r = MatrixXd::Zero(nPoints, n);
	MatrixXd temp_A;
	for (int i = 0; i < nPoints; i++)
	{
		for (int j = 0; j < n; j++)
		{
			r(i, j) = (x.row(i) - nodes.row(j)).norm();
		}
	}
	temp_A = RBFFunction(r, RBFConstant);
	MatrixXd P(x.rows(), x.cols() + 1);
	P.leftCols(1) = MatrixXd::Ones(x.rows(), 1);
	P.rightCols(x.cols()) = x;
	MatrixXd A(nPoints, n + x.cols() + 1);
	A.topLeftCorner(temp_A.rows(), temp_A.cols()) = temp_A;
	A.topRightCorner(P.rows(), P.cols()) = P;
	MatrixXd f;
	f = A*rbfcoeff;
	return f;
}

double rbfcheck(MatrixXd x, MatrixXd y, MatrixXd rbfcoeff,
	MatrixXd RBFFunction(MatrixXd r, double const_num),
	double RBFConstant)
	//计算采样点的采样值和插值值的最大误差s
{
	MatrixXd S;
	S = rbfinterp(x, rbfcoeff, x, RBFFunction, RBFConstant);
	double maxdiff;
	maxdiff = (S - y).cwiseAbs().maxCoeff();
	return maxdiff;
}