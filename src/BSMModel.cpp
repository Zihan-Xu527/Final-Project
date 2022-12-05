//
// Created by Zihan Xu on 11/29/22.
//

#include "BSMModel.h"
#include "Payoff.h"

std::vector<std::vector<double>> BSMModel::sim_path(double maturity, double S_init, double strike, unsigned length_out)
{
//    std::vector<double> path;

//
//    // trivial case
//    if (length_out == 1) {
//        return path;
//    }

    // first simulate path of log S, then take exp (numerically more stable)
//    path[0] = std::log(S_init);

//    unsigned num_period = length_out - 1;
//    double growth_rate = interest_rate - 0.5*sigma*sigma;
//    double dt = maturity/num_period;
//    double drift = growth_rate*dt;
//    double volatility = sigma*std::sqrt(dt);

//    for (unsigned i = 1; i <= num_period; ++i) {
//        double dlogS = drift + volatility*normal(rng);
//        path[i] = path[i - 1] + dlogS;

//    }


    std::vector<std::vector<double>> path(length_out, std::vector<double> (batch_size, std::log(S_init)));
    Call my_call(strike);

    // trivial case
    if (length_out > 1)
    {
        unsigned num_period = length_out - 1;
        double growth_rate = interest_rate - 0.5 * sigma * sigma;
        double dt = maturity / num_period;
        double drift = growth_rate * dt;
        double volatility = sigma * std::sqrt(dt);
//#pragma omp parallel for
        for (int i = 1; i <= num_period; ++i) {
            for (int j = 0; j<batch_size; j++){
                double dlogS = drift + volatility * normal(rng);
                path[i][j]=path[i-1][j] + dlogS;
            }
        }

    }
    for (int i=0; i< length_out; i++){
        for (int j = 0; j<batch_size; j++){
            path[i][j] = std::exp(path[i][j]);
        }
    }
    // get S = exp(log S)
//    for (unsigned i = 0; i < length_out ; ++i)
//        path[i] = std::exp(path[i]);
//
    return path;
}