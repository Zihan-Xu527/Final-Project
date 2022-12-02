//
// Created by Zihan Xu on 11/28/22.
//

#ifndef FINAL_PROJECT_STATISTIC_TOOL_H
#define FINAL_PROJECT_STATISTIC_TOOL_H
#include <vector>
#include <numeric>


double average(std::vector<double> const & vec);
double variance(std::vector<double> const & vec);
double std_deviation(std::vector<double> const & vec);
double std_error(std::vector<double> const & vec);
std::vector<double> confidenceInterval(std::vector<double> const & vec);
double get_b(std::vector<double> const & spot, std::vector<double> const & payoff);
std::vector<double> getYb(std::vector<double> const & spot, std::vector<double> const & payoff, double expected_spot);
double normalCDF(double x);
double expected_Payoff(double t, double spot_price, double interest_rate, double volatility, double strike, double maturity);
#endif //FINAL_PROJECT_STATISTIC_TOOL_H
