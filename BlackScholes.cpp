#define _USE_MATH_DEFINES
#include<iostream>
#include<cmath>
#include<fstream>

struct Contract {
	double time_to_mature;
	double strike_price;
	double current_price;
	double volatility;
	double risk_free_int;
	bool isCallOption;

	Contract(double ttm, double sp, double current, double vol, double risk, bool is) {
		time_to_mature = ttm;
		strike_price = sp;
		current_price = current;
		volatility = vol;
		risk_free_int = risk;
		isCallOption = is; 
	}

};

struct Greeks {
	double delta;
	double gamma;
	double vega;
	double theta;
	double rho;
};

//standard norm. probability density function
double norm_pdf(const double& x) {
	return (1.0 / (pow(2 * M_PI, 0.5))) * exp(-0.5 * x * x);
}

//approx. of the normal cumulative distribution function. 
double norm_cdf(const double& x) {
	double k = 1.0 / (1.0 + 0.2316419 * x);
	double k_sum = k * (0.319381530 + k * (-0.356563782 + k * (1.781477937 + k * (-1.821255978 + 1.330274429 * k))));

	if (x >= 0.0) {
		return (1.0 - (1.0 / (pow(2 * M_PI, 0.5))) * exp(-0.5 * x * x) * k_sum);
	}
	else {
		return 1.0 - norm_cdf(-x);
	}
}

static double Getd1(Contract& params) {
	double x;
	double add1 = std::log(params.current_price / params.strike_price);
	double add2 = (params.risk_free_int + (0.5 * ((std::pow(params.volatility, 2.0)))) * params.time_to_mature);
	double divisor = params.volatility * (std::sqrt(params.time_to_mature));
	return x = (add1 + add2) / divisor;
}

static float Getd2(Contract& params) {
	float firstVar = Getd1(params);
	return firstVar - params.volatility * (std::sqrt(params.time_to_mature));
}

double PriceBlackScholesModel(Contract& params) {
	double d1;
	double d2;
	double eRT = std::pow(M_E, (-params.risk_free_int * params.time_to_mature)); // Variable used to represent (e^-rt)

	d1 = Getd1(params);
	d2 = Getd2(params);

	//If the option is a call use C = N(d1)St - N(d2)Ke^-rt ---> Call option price. 
	if (params.isCallOption) {
		float call_price = norm_cdf(d1) * params.current_price - norm_cdf(d2) * params.strike_price * eRT;
		return call_price;
	}
	else { // compute for put option
		float put_price = params.strike_price * eRT * norm_cdf(-d2) - params.current_price * norm_cdf(-d1);
		return put_price;

	}
}

Greeks calculateGreeks(Contract& params) {
	Greeks greekList;
	double d1 = Getd1(params);
	double d2 = Getd2(params);

	if (params.isCallOption) {
		greekList.delta = norm_cdf(d1);
		greekList.gamma = (norm_pdf(d1) / (params.current_price * params.volatility * std::sqrt(params.time_to_mature)));
		greekList.vega = params.current_price * norm_pdf(d1) * std::sqrt(params.time_to_mature);

		//Calulations for variables in theta. 
		double CallThetaVar = -(params.strike_price * norm_pdf(d1) * params.volatility) / (2 * std::sqrt(params.time_to_mature));
		double CallThetaVar2 = params.risk_free_int * params.strike_price * (std::pow(M_E, -params.risk_free_int * params.time_to_mature)) * norm_cdf(d2);

		greekList.theta = CallThetaVar - CallThetaVar2;

		greekList.rho = params.strike_price * params.time_to_mature * (std::pow(M_E, -params.risk_free_int * params.time_to_mature)) * norm_cdf(d2);
	}
	else {
		greekList.delta = norm_cdf(d1) - 1;
		greekList.gamma = (norm_pdf(d1) / (params.current_price * params.volatility * std::sqrt(params.time_to_mature)));
		greekList.vega = params.current_price * norm_pdf(d1) * std::sqrt(params.time_to_mature);

		double PutThetaVar = -(params.strike_price * norm_pdf(d1) * params.volatility) / std::sqrt((2 * params.time_to_mature));
		double PutThetaVar2 = params.risk_free_int * params.strike_price * (std::pow(M_E, -params.risk_free_int * params.time_to_mature)) * norm_cdf(-d2);

		greekList.theta = PutThetaVar + PutThetaVar2;

		greekList.rho = -params.strike_price * params.time_to_mature * (std::pow(M_E, -params.risk_free_int * params.time_to_mature)) * norm_cdf(-d2);
	}
	return greekList; 
}
int main() {
	std::ofstream file("black_scholes_output.csv");
	// Write header row
	double S = 100.0, K = 100.0, T = 1.0, r = 0.05, sigma = 0.2;
	Contract Callparams = Contract(S, K, T, r, sigma, true);
	return 0;
} 
