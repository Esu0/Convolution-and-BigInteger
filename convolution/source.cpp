#include<iostream>
#include<vector>
#include<ios>
#include"convolution.h"
#include<iomanip>
#include<convenience.h>
#include<fourier_trans.h>

using namespace convenient;
using namespace instant;
int main()
{
	std::vector<long long> a = { 3,2,14,1,3,4,1,4 };
	number_theoretic_transform<17, 5>(a.data(), a.size());
	std::cout << to_string(a) << std::endl;
	number_theoretic_transform<17, 5, true>(a.data(), a.size());
	std::cout << to_string(a) << std::endl;
}