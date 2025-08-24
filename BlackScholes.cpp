#define _USE_MATH_DEFINES
#include<iostream>
#include<cmath>
#include<fstream>

struct Contract {
	float time_to_mature;
	float strike_price;
	float current_price;
	float volatility;
	float risk_free_int;

	Contract(float ttm, float sp, float current, float vol, float risk) {
		time_to_mature = ttm;
		strike_price = sp;
		current_price = current;
		volatility = vol;
		risk_free_int = risk;
	}

};

struct Greeks {
	float delta;
	float gamma;
	float vega;
	float theta;
	float rho;
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

static float Getd1(Contract& params) {
	float x;
	float add1 = std::log(params.current_price / params.strike_price);
	float add2 = (params.risk_free_int + (0.5 * ((std::pow(params.volatility, 2.0)))) * params.time_to_mature);
	float divisor = params.volatility * (std::sqrt(params.time_to_mature));
	return x = (add1 + add2) / divisor;
}

static float Getd2(Contract& params) {
	float firstVar = Getd1(params);
	return firstVar - params.volatility * (std::sqrt(params.time_to_mature));
}

double PriceBlackScholesModel(Contract& params, bool isCallOption) {
	float d1;
	float d2;
	float eRT = std::pow(M_E, (-params.risk_free_int * params.time_to_mature)); // Variable used to represent (e^-rt)

	d1 = Getd1(params);
	d2 = Getd2(params);

	//If the option is a call use C = N(d1)St - N(d2)Ke^-rt ---> Call option price. 
	if (isCallOption) {
		float call_price = norm_cdf(d1) * params.current_price - norm_cdf(d2) * params.strike_price * eRT;
		return call_price;
	}
	else { // compute for put option
		float put_price = params.strike_price * eRT * norm_cdf(-d2) - params.current_price * norm_cdf(-d1);
		return put_price;

	}
}

Greeks calculateGreeks(Contract& params, bool isCallOption) {
	Greeks greekList;
	float d1 = Getd1(params);
	float d2 = Getd2(params);

	if (isCallOption){
		greekList.delta = std::pow(M_E, (-params.risk_free_int * params.time_to_mature)) * norm_cdf(d1);
		greekList.gamma = (norm_pdf(d1) * std::pow(M_E, (-params.risk_free_int * params.time_to_mature))) /
			(params.current_price * params.volatility * std::sqrt(params.time_to_mature));
	}
	else {

	}
}
int main() {
	//implement user inputs and csv file output. 
}
