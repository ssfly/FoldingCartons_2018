#pragma once
#include <string>

using namespace std;

struct Parameters
{
	double MyPi;
	string InputFilePath;
	double MyEpsilon;
	int MyInfinity;
	int ResolutionMultipe;
	double ThresholdForSeparateLines;
};



extern Parameters g_parameters;
