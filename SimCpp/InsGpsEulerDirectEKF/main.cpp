// InsGpsEulerDirectEKF.cpp : ���ļ����� "main" ����������ִ�н��ڴ˴���ʼ��������
//


#include "main.h"
#include <time.h>

using namespace matrix;
static struct UavSensors uavSensors;
static std::vector<UavSensors> vUavSensors;

DataType outcome[15] = { 0 };
int outcome_dim = 15;

int main()
{
	char sensor_data_path[] = "./DATA/UavData_SquareSpiral.txt";
	prepare_data(sensor_data_path);

	clsInsGpsLooseEulerDirectEKF<DataType, 15, 1, 6> InsGpsLooseEulerDirectEKF;

	DataType pQt_data[15] = { 4.00000000000000e-06,4.00000000000000e-06,4.00000000000000e-06,2.50000000000000e-05,2.50000000000000e-05,2.50000000000000e-05,1.00000000000000e-06,1.00000000000000e-06,1.00000000000000e-06,1.00000000000000e-06,1.00000000000000e-06,1.00000000000000e-06,1.00000000000000e-06,1.00000000000000e-06,1.00000000000000e-06 };
	int Qt_dim = 15;

	DataType pR_yaw_data[1] = { 0.0200000000000000 };
	int Ryaw_dim = 1;

	DataType pR_gps_data[6] = { 1,1,1,0.00250000000000000,0.00250000000000000,0.00250000000000000 };
	int Rgps_dim = 6;

	clock_t time_start = clock();
	for (uint32_t i = 0; i < vUavSensors.size(); i++)
	{
		// �˲�ʱ����
		DataType time_s = vUavSensors.at(i).time_s;
		DataType dt;
		if (i == 0)	
			dt = vUavSensors.at(1).time_s - vUavSensors.at(0).time_s;
		else
			dt = time_s - vUavSensors.at(i-1).time_s;
		// IMU���ݸ���
		DataType gx_rps = vUavSensors.at(i).gx_rps;
		DataType gy_rps = vUavSensors.at(i).gy_rps;
		DataType gz_rps = vUavSensors.at(i).gz_rps;
		DataType fx_mps2 = vUavSensors.at(i).fx_mps2;
		DataType fy_mps2 = vUavSensors.at(i).fy_mps2;
		DataType fz_mps2 = vUavSensors.at(i).fz_mps2;
		// ��������ݸ���
		int isYawUpdate = vUavSensors.at(i).isYawUpdate;
		DataType yaw_rad = vUavSensors.at(i).yaw_rad;
		// GPS���ݸ���
		int isGPSUpdate = vUavSensors.at(i).isGPSUpdate;
		DataType Pn_m = vUavSensors.at(i).Pn_m;
		DataType Pe_m = vUavSensors.at(i).Pe_m;
		DataType h_msl_m = vUavSensors.at(i).h_msl_m;
		DataType Vn_mps = vUavSensors.at(i).Vn_mps;
		DataType Ve_mps = vUavSensors.at(i).Ve_mps;
		DataType Vd_mps = vUavSensors.at(i).Vd_mps;
		// �˲�����
		InsGpsLooseEulerDirectEKF.Run(
			outcome, outcome_dim,				// �˲����
			//vUavSensors.at(i).outcome, outcome_dim,
			dt,									// ʱ����
			gx_rps, gy_rps, gz_rps,				// ������
			fx_mps2, fy_mps2, fz_mps2,			// ���ٶȼ� 
			isYawUpdate, yaw_rad,				// ����ǲ���ֵ
			isGPSUpdate, Pn_m, Pe_m, h_msl_m, Vn_mps, Ve_mps, Vd_mps,	// GPS����ֵ
			pQt_data, Qt_dim,					// Q ģ����������
			pR_yaw_data, Ryaw_dim,				// R_yaw ����ǲ�����������
			pR_gps_data, Rgps_dim				// R_gps λ���ٶȲ�����������
		);
		// ��ӡ�˲����
		for (uint32_t j = 0; j < outcome_dim; j++)
		{
			vUavSensors.at(i).outcome[j] = outcome[j];
			//if (j == 0)	printf("[ ");
			//printf("%0.5lf ", vUavSensors.at(i).outcome[j]);
			//if (j == (outcome_dim - 1))	printf("]\n");
		}
	}
	clock_t time_end = clock();
	// �����ƽ��д�롰./DATA/EstInsGpsEulerDirectEKF.txt��
	char est_data_path[] = "./DATA/EstInsGpsEulerDirectEKF.txt";
	save_estvaule(est_data_path);
	std::cout << "Running Time (s): " << (1.0*time_end - 1.0*time_start) / (1.0*CLOCKS_PER_SEC) << std::endl;
	std::cout << "Hello World!\n" << std::endl;
	return 0;
}

int prepare_data(const char * datapath)
{
	// ׼�����ݣ�(MATLAB)��DATA�ļ���������generate_uav_sensors_measurement_for_cpp.m�ļ������� UavData_SquareSpiral.txt
	// ��ȡ����
	std::ifstream f_read(datapath);
	//std::ifstream f_read("./DATA/UavData_SquareSpiral.txt");
	if (!f_read)
	{
		std::cout << "ifstream - Cannot open ./DATA/UavData_SquareSpiral.txt" << std::endl;
		std::cout << "Exit :(" << std::endl;
		return -1;
	}
	std::string str_fread;
	while (getline(f_read, str_fread))
	{
		//cout << str_fread << endl;
		// �Ѷ��Ż��ɿո�
		//for(std::string::iterator iter = str_fread.begin(); iter!=str_fread.end();iter++)
		for (uint32_t i = 0; i < str_fread.length(); i++)
		{
			if (str_fread.at(i) == ',')	str_fread.at(i) = ' ';
			//if (*iter == ',') str_fread = ' ';
		}
		//cout << str_fread << endl;
		// ��ȡ���֣�����istringstream��ͷ�ļ�<sstream>��
		std::istringstream out(str_fread);
		std::string tmp;
		while (out >> tmp)
		{
			// ���ַ���ת��Ϊ���� stoi stof stod
			static int stage = 0;	// �趨 stage ��������
			if (stage == 0)
			{
				uavSensors.time_s = stod(tmp);	// time_s
				stage++;
				continue;
			}
			if (stage == 1)
			{
				uavSensors.gx_rps = stod(tmp);
				stage++;
				continue;
			}
			if (stage == 2)
			{
				uavSensors.gy_rps = stod(tmp);
				stage++;
				continue;
			}
			if (stage == 3)
			{
				uavSensors.gz_rps = stod(tmp);
				stage++;
				continue;
			}
			if (stage == 4)
			{
				uavSensors.fx_mps2 = stod(tmp);
				stage++;
				continue;
			}
			if (stage == 5)
			{
				uavSensors.fy_mps2 = stod(tmp);
				stage++;
				continue;
			}
			if (stage == 6)
			{
				uavSensors.fz_mps2 = stod(tmp);
				stage++;
				continue;
			}
			if (stage == 7)
			{
				uavSensors.isYawUpdate = stoi(tmp);
				stage++;
				continue;
			}
			if (stage == 8)
			{
				uavSensors.yaw_deg = stod(tmp);
				uavSensors.yaw_rad = uavSensors.yaw_deg * 3.1415926535f / 180.0f;
				stage++;
				continue;
			}
			if (stage == 9)
			{
				uavSensors.isGPSUpdate = stoi(tmp);
				stage++;
				continue;
			}
			if (stage == 10)
			{
				uavSensors.Pn_m = stod(tmp);
				stage++;
				continue;
			}
			if (stage == 11)
			{
				uavSensors.Pe_m = stod(tmp);
				stage++;
				continue;
			}
			if (stage == 12)
			{
				uavSensors.h_msl_m = stod(tmp);
				stage++;
				continue;
			}
			if (stage == 13)
			{
				uavSensors.Vn_mps = stod(tmp);
				stage++;
				continue;
			}
			if (stage == 14)
			{
				uavSensors.Ve_mps = stod(tmp);
				stage++;
				continue;
			}
			if (stage == 15)
			{
				uavSensors.Vd_mps = stod(tmp);
				stage = 0;
				// continue;
			}
			//printf("%f,%f,%f,%f,%f,%f,%f,%d,%f,%d,%f,%f,%f,%f,%f,%f\n", 
			//	uavSensors.time_s,
			//	uavSensors.gx_rps, uavSensors.gy_rps, uavSensors.gz_rps,
			//	uavSensors.fx_mps2, uavSensors.fy_mps2, uavSensors.fz_mps2,
			//	uavSensors.isYawUpdate, uavSensors.yaw_rad,
			//	uavSensors.isGPSUpdate, uavSensors.Pn_m, uavSensors.Pe_m, uavSensors.h_msl_m, uavSensors.Vn_mps, uavSensors.Ve_mps, uavSensors.Vd_mps);
			// �����ݻ�����vector��
			vUavSensors.push_back(uavSensors);
		}
	}
	f_read.close();
	std::cout << "Total count: " << vUavSensors.size() << std::endl;
	//for (int i = 0; i < vUavSensors.size(); i++)
	//{
	//	printf("%f,%f,%f,%f,%f,%f,%f,%d,%f,%d,%f,%f,%f,%f,%f,%f\n",
	//		vUavSensors.at(i).time_s,
	//		vUavSensors.at(i).gx_rps, vUavSensors.at(i).gy_rps, vUavSensors.at(i).gz_rps,
	//		vUavSensors.at(i).fx_mps2, vUavSensors.at(i).fy_mps2, vUavSensors.at(i).fz_mps2,
	//		vUavSensors.at(i).isYawUpdate, vUavSensors.at(i).yaw_rad,
	//		vUavSensors.at(i).isGPSUpdate, vUavSensors.at(i).Pn_m, vUavSensors.at(i).Pe_m, vUavSensors.at(i).h_msl_m, vUavSensors.at(i).Vn_mps, vUavSensors.at(i).Ve_mps, vUavSensors.at(i).Vd_mps);
	//}
	return 0;
}

int save_estvaule(const char* datapath)
{
	std::ofstream fout("./DATA/EstInsGpsEulerDirectEKF.txt");
	if (!fout)
	{
		std::cout << "ofsteam - Cannot create and open EstInsGpsEulerDirectEKF.txt" << std::endl;
		return -1;
	}
	for (uint32_t i = 0; i < vUavSensors.size(); i++)
	{
		fout << vUavSensors.at(i).time_s << ",";
		for (uint32_t j = 0; j < outcome_dim; j++)
		{
			if (j < outcome_dim-1)
			{
				fout << vUavSensors.at(i).outcome[j] << ",";
			}
			if (j == outcome_dim-1)	fout << vUavSensors.at(i).outcome[j] << std::endl;
		}
	}
	fout.close();
	return 0;
}
