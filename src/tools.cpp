#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{

    Eigen::VectorXd rmse(4);
    rmse << 0,0,0,0;
    if(estimations.size() != ground_truth.size() || estimations.size() == 0)
    {
        std::cerr << "Wrong size: Size is not equal to the ground truth or equal to zero" << endl;
        return rmse;
    }

    for(uint i=0; i < estimations.size(); ++i)
    {
        Eigen::VectorXd tmp = estimations[i] - ground_truth[i];

        tmp = tmp.array()*tmp.array();
        rmse += tmp;
    }

    rmse = rmse/estimations.size();

    rmse = rmse.array().sqrt();

    return rmse;

}
