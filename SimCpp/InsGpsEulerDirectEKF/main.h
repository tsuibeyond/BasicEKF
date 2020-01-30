#pragma once

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include "InsGpsEulerDirectEKF.hpp"
#include "matrix/math.hpp"

typedef float DataType;

struct UavSensors
{
	DataType outcome[15];// = { 0 };
	DataType time_s;


	DataType gx_rps;
	DataType gy_rps;
	DataType gz_rps;
	DataType fx_mps2;
	DataType fy_mps2;
	DataType fz_mps2;

	int isYawUpdate;
	DataType yaw_rad;
	DataType yaw_deg;

	int isGPSUpdate;
	DataType Pn_m;
	DataType Pe_m;
	DataType h_msl_m;
	DataType Vn_mps;
	DataType Ve_mps;
	DataType Vd_mps;
};

int prepare_data(const char * datapath);

int save_estvaule(const char* datapath);

