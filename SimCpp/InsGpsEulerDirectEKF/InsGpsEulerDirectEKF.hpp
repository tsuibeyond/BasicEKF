#pragma once

#include <array>
#include <math.h>
#include "matrix/math.hpp"

using namespace matrix;

template<typename Type, size_t N, size_t M1, size_t M2>	//N表示状态向量维数，M1表示航向角测量， M2表示GPS测量值（三轴速度和三轴位置）
class clsInsGpsLooseEulerDirectEKF
{
public:
	clsInsGpsLooseEulerDirectEKF() {
		isGPSvalid = 0;
		isMagYawVaild = 0;
		Gravity(0, 0) = 0.0f;	Gravity(1, 0) = 0.0f;	Gravity(2, 0) = g_mps2;
		I_nxn.identity();
		I_2nx2n.identity();
		wx = 0.0f, wy = 0.0f, wz = 0.0f;
		fx = 0.0f, fy = 0.0f, fz = 0.0f;
	}

	~clsInsGpsLooseEulerDirectEKF() {
		// 导航&控制算法禁用动态内存
	}

	// 常值
	const Type PI = 3.141592653f;
	const Type g_mps2 = 9.81f;
	Matrix<Type, 3, 1> Gravity;
	// 
	int isGPSvalid;
	int isMagYawVaild;
	// 
	Type wx, wy, wz;
	Type fx, fy, fz;
	Matrix<Type, N, 1> xdot;	//	状态向量的导数
	Matrix<Type, N, 1> xhat;	//	状态向量
	Matrix<Type, N, N> F;		// 
	Matrix<Type, N * 2, N * 2> AA;		// 
	Matrix<Type, N * 2, N * 2> BB;		// 
	Matrix<Type, N, N> PHI;		//	一步预测矩阵
	Matrix<Type, N, N> Qt;		//	过程噪声
	Matrix<Type, N, N> Qk;		//	过程噪声
	Matrix<Type, N, N> P;		//	估计均方误差
	Matrix<Type, M1, 1> Z_yaw;		//	航向角测量值
	Matrix<Type, M1, N> H_yaw;		//	航向角测量矩阵
	Matrix<Type, M1, M1> R_yaw;		//	航向角测量方差
	Matrix<Type, N, M1> K_yaw;		//	航向角与一步预测的状态量的融合权重
	Matrix<Type, M2, 1> Z_gps;		//	GPS测量值
	Matrix<Type, M2, N> H_gps;		//	GPS测量矩阵
	Matrix<Type, M2, M2> R_gps;		//	GPS测量方差
	Matrix<Type, N, M2> K_gps;		//	GPS测量值与一步预测的状态量的融合权重
	Matrix<Type, N, N> I_nxn;			// 单位矩阵
	Matrix<Type, N * 2, N * 2> I_2nx2n;			// 单位矩阵
	Matrix<Type, 3, 1> rpyDot;		// 姿态角导数
	Matrix<Type, 3, 1> GyroRmBias;	//
	Matrix<Type, 3, 1> AccRmBias;	//
	Matrix<Type, 3, 1> PosNEDdot;	//
	Matrix<Type, 3, 1> VelNEDdot;	//
	Matrix<Type, 3, 3> C_bodyrate2eulerdot;
	Matrix<Type, 3, 3> C_ned2b;
	Matrix<Type, 3, 3> C_b2ned;
	
	inline Type sec(Type x)
	{
		return 1.0f / cos(x);
	}
	// IMU数据更新
	void imu_measure(Type gx_rps, Type gy_rps, Type gz_rps, Type fx_mps2, Type fy_mps2, Type fz_mps2)
	{
		wx = gx_rps;	wy = gy_rps;	wz = gz_rps;
		fx = fx_mps2;	fy = fy_mps2;	fz = fz_mps2;
	}
	// 航向角测量更新
	void yaw_measure(int isUpdate, Type yaw_rad)
	{
		if (isUpdate == 1)
		{
			Z_yaw(0, 0) = yaw_rad;
			//printf("%lf\n",Z_yaw(0,0));
			isMagYawVaild = 1;
		}
		else
		{
			isMagYawVaild = 0;
		}
	}
	// GPS(三轴位置、速度更新)
	void gps_measure(int isUpdate, Type Pn_m, Type Pe_m, Type h_msl_m, Type Vn_mps, Type Ve_mps, Type Vd_mps)
	{
		if (isUpdate == 1)
		{
			Z_gps(0, 0) = Pn_m;	Z_gps(1, 0) = Pe_m;	Z_gps(2, 0) = h_msl_m;
			Z_gps(3, 0) = Vn_mps;	Z_gps(4, 0) = Ve_mps;	Z_gps(5, 0) = Vd_mps;
			isGPSvalid = 1;
		}
		else
		{
			isGPSvalid = 0;
		}
	}
	// 
	void init()
	{
		init_xhat();	// 状态向量初始化
		init_P();		// 状态估计协方差初始化
	}

	inline void init_xhat()
	{
		xhat.zero();
		// 仅仿真时使用
		xhat(2, 0) = -2.72343926657247f;
		xhat(3, 0) = -0.0161376263168116f;
		xhat(4, 0) = 0.104874716016494f;
		xhat(5, 0) = 109.623435014780f;
		xhat(6, 0) = -18.1430696059217f;
		xhat(7, 0) = -3.30602783578041f;
		xhat(8, 0) = 0.0849970564432384f;
	}

	inline void init_P()
	{
		P.zero();
		// 姿态角估计协方差初始化
		P(0, 0) = 30.0f * PI / 180.0f;
		P(1, 1) = 30.0f * PI / 180.0f;
		P(2, 2) = 30.0f * PI / 180.0f;
		// 位置估计协方差初始化
		P(3, 3) = 100.0f;
		P(4, 4) = 100.0f;
		P(5, 5) = 100.0f;
		// 速度估计协方程初始化
		P(6, 6) = 10.0f;
		P(7, 7) = 10.0f;
		P(8, 8) = 10.0f;
		// 三轴陀螺仪常值偏差估计协方程初始化
		P(9, 9) = 0.01f;
		P(10, 10) = 0.01f;
		P(11, 11) = 0.01f;
		// 三轴加速度计常值偏差估计协方程初始化
		P(12, 12) = 0.05f;
		P(13, 13) = 0.05f;
		P(14, 14) = 0.05f;
		// 对其平方
		P = P * P;
	}

	void updateXdot()
	{
		updateRPYdot();
		xdot(0, 0) = rpyDot(0, 0);	xdot(1, 0) = rpyDot(1, 0);	xdot(2, 0) = rpyDot(2, 0);

		updatePosNEDdot();
		xdot(3, 0) = PosNEDdot(0, 0);	xdot(4, 0) = PosNEDdot(1, 0);	xdot(5, 0) = PosNEDdot(2, 0);

		updateVelNEDdot();
		xdot(6, 0) = VelNEDdot(0, 0);	xdot(7, 0) = VelNEDdot(1, 0);	xdot(8, 0) = VelNEDdot(2, 0);

		xdot(9, 0) = 0.0f;	xdot(10, 0) = 0.0f;	xdot(11, 0) = 0.0f;

		xdot(12, 0) = 0.0f;	xdot(13, 0) = 0.0f;	xdot(14, 0) = 0.0f;
	}
	//
	inline void updateRPYdot()
	{
		updateC_bodyrate2eulerdot();
		updateGyroRmBias();
		rpyDot = C_bodyrate2eulerdot * GyroRmBias;
	}
	inline void updateC_bodyrate2eulerdot()
	{
		Type phi = xhat(0, 0);
		Type theta = xhat(1, 0);
		Type psi = xhat(2, 0);
		C_bodyrate2eulerdot(0, 0) = 1.0f;	C_bodyrate2eulerdot(0, 1) = sin(phi) * tan(theta);	C_bodyrate2eulerdot(0, 2) = cos(phi) * tan(theta);
		C_bodyrate2eulerdot(1, 0) = 0.0f;	C_bodyrate2eulerdot(1, 1) = cos(phi);	C_bodyrate2eulerdot(1, 2) = -sin(phi);
		C_bodyrate2eulerdot(2, 0) = 0.0f;	C_bodyrate2eulerdot(2, 1) = sin(phi) * sec(theta);	C_bodyrate2eulerdot(2, 2) = cos(phi) * sec(theta);
	}
	inline void updateGyroRmBias()
	{
		GyroRmBias(0, 0) = wx - xhat(9, 0);	GyroRmBias(1, 0) = wy - xhat(10, 0);	GyroRmBias(2, 0) = wz - xhat(11, 0);
	}
	//
	inline void updatePosNEDdot()
	{
		Type Vn = xhat(6, 0);
		Type Ve = xhat(7, 0);
		Type Vd = xhat(8, 0);
		PosNEDdot(0, 0) = Vn;	PosNEDdot(1, 0) = Ve;	PosNEDdot(2, 0) = -Vd;
	}
	//
	inline void updateVelNEDdot()
	{
		updateC_ned2b();
		updateAccRmBias();
		VelNEDdot = C_b2ned * AccRmBias + Gravity;
	}
	inline void updateC_ned2b()
	{
		Type phi = xhat(0, 0);
		Type theta = xhat(1, 0);
		Type psi = xhat(2, 0);
		C_ned2b(0, 0) = cos(theta) * cos(psi);	C_ned2b(0, 1) = cos(theta) * sin(psi);	C_ned2b(0, 2) = -sin(theta);
		C_ned2b(1, 0) = sin(phi) * sin(theta) * cos(psi) - cos(phi) * sin(psi);	C_ned2b(1, 1) = sin(phi) * sin(theta) * sin(psi) + cos(phi) * cos(psi);	C_ned2b(1, 2) = sin(phi) * cos(theta);
		C_ned2b(2, 0) = cos(phi) * sin(theta) * cos(psi) + sin(phi) * sin(psi);	C_ned2b(2, 1) = cos(phi) * sin(theta) * sin(psi) - sin(phi) * cos(psi);	C_ned2b(2, 2) = cos(phi) * cos(theta);
		C_b2ned = C_ned2b.transpose();
	}
	inline void updateAccRmBias(void)
	{
		AccRmBias(0, 0) = fx - xhat(12, 0);	AccRmBias(1, 0) = fy - xhat(13, 0);	AccRmBias(2, 0) = fz - xhat(14, 0);
	}
	//
	void updateF()
	{
		Type phi = xhat(0, 0);	Type theta = xhat(1, 0);	Type psi = xhat(2, 0);
		Type bwx = xhat(9, 0);	Type bwy = xhat(10, 0);	Type bwz = xhat(11, 0);
		Type bax = xhat(12, 0);	Type bay = xhat(13, 0);	Type baz = xhat(14, 0);
		Type F_data[] = {
																					 sin(phi)* tan(theta)* (bwz - wz) - cos(phi) * tan(theta) * (bwy - wy),                                  -cos(phi) * (bwz - wz) * (tan(theta)* tan(theta) + 1) - sin(phi) * (bwy - wy) * (tan(theta)* tan(theta) + 1),                                                                                                                                                              0, 0, 0, 0, 0, 0,  0, -1, -sin(phi) * tan(theta), -cos(phi) * tan(theta),                    0,                                                  0,                                                  0,
																							   cos(phi)* (bwz - wz) + sin(phi) * (bwy - wy),                                                                                                                  0,                                                                                                                                                              0, 0, 0, 0, 0, 0,  0,  0,            -cos(phi),             sin(phi),                    0,                                                  0,                                                  0,
																	 (sin(phi) * (bwz - wz)) / cos(theta) - (cos(phi) * (bwy - wy)) / cos(theta),                    -(cos(phi) * sin(theta) * (bwz - wz)) / (cos(theta)* cos(theta)) - (sin(phi) * sin(theta) * (bwy - wy)) / (cos(theta)*cos(theta)),                                                                                                                                                              0, 0, 0, 0, 0, 0,  0,  0, -sin(phi) / cos(theta), -cos(phi) / cos(theta),                    0,                                                  0,                                                  0,
																																	   0,                                                                                                                  0,                                                                                                                                                              0, 0, 0, 0, 1, 0,  0,  0,                    0,                    0,                    0,                                                  0,                                                  0,
																																	   0,                                                                                                                  0,                                                                                                                                                              0, 0, 0, 0, 0, 1,  0,  0,                    0,                    0,                    0,                                                  0,                                                  0,
																																	   0,                                                                                                                  0,                                                                                                                                                              0, 0, 0, 0, 0, 0, -1,  0,                    0,                    0,                    0,                                                  0,                                                  0,
		 -(sin(phi) * sin(psi) + cos(phi) * cos(psi) * sin(theta)) * (bay - fy) - (cos(phi) * sin(psi) - cos(psi) * sin(phi) * sin(theta)) * (baz - fz), cos(psi)* sin(theta)* (bax - fx) - cos(phi) * cos(psi) * cos(theta) * (baz - fz) - cos(psi) * cos(theta) * sin(phi) * (bay - fy), (cos(phi) * cos(psi) + sin(phi) * sin(psi) * sin(theta))* (bay - fy) - (cos(psi) * sin(phi) - cos(phi) * sin(psi) * sin(theta)) * (baz - fz) + cos(theta) * sin(psi) * (bax - fx), 0, 0, 0, 0, 0,  0,  0,                    0,                    0, -cos(psi) * cos(theta),   cos(phi)* sin(psi) - cos(psi) * sin(phi) * sin(theta), -sin(phi) * sin(psi) - cos(phi) * cos(psi) * sin(theta),
		   (cos(psi) * sin(phi) - cos(phi) * sin(psi) * sin(theta))* (bay - fy) + (cos(phi) * cos(psi) + sin(phi) * sin(psi) * sin(theta)) * (baz - fz), sin(psi)* sin(theta)* (bax - fx) - cos(phi) * cos(theta) * sin(psi) * (baz - fz) - cos(theta) * sin(phi) * sin(psi) * (bay - fy), (cos(phi) * sin(psi) - cos(psi) * sin(phi) * sin(theta))* (bay - fy) - (sin(phi) * sin(psi) + cos(phi) * cos(psi) * sin(theta)) * (baz - fz) - cos(psi) * cos(theta) * (bax - fx), 0, 0, 0, 0, 0,  0,  0,                    0,                    0, -cos(theta) * sin(psi), -cos(phi) * cos(psi) - sin(phi) * sin(psi) * sin(theta),   cos(psi)* sin(phi) - cos(phi) * sin(psi) * sin(theta),
																		 cos(theta)* sin(phi)* (baz - fz) - cos(phi) * cos(theta) * (bay - fy),                            cos(theta)* (bax - fx) + cos(phi) * sin(theta) * (baz - fz) + sin(phi) * sin(theta) * (bay - fy),                                                                                                                                                              0, 0, 0, 0, 0, 0,  0,  0,                    0,                    0,           sin(theta),                               -cos(theta) * sin(phi),                               -cos(phi) * cos(theta),
																																	   0,                                                                                                                  0,                                                                                                                                                              0, 0, 0, 0, 0, 0,  0,  0,                    0,                    0,                    0,                                                  0,                                                  0,
																																	   0,                                                                                                                  0,                                                                                                                                                              0, 0, 0, 0, 0, 0,  0,  0,                    0,                    0,                    0,                                                  0,                                                  0,
																																	   0,                                                                                                                  0,                                                                                                                                                              0, 0, 0, 0, 0, 0,  0,  0,                    0,                    0,                    0,                                                  0,                                                  0,
																																	   0,                                                                                                                  0,                                                                                                                                                              0, 0, 0, 0, 0, 0,  0,  0,                    0,                    0,                    0,                                                  0,                                                  0,
																																	   0,                                                                                                                  0,                                                                                                                                                              0, 0, 0, 0, 0, 0,  0,  0,                    0,                    0,                    0,                                                  0,                                                  0,
																																	   0,                                                                                                                  0,                                                                                                                                                              0, 0, 0, 0, 0, 0,  0,  0,                    0,                    0,                    0,                                                  0,                                                  0
		};
		//
		for (int i = 0; i <= N - 1; i++)
		{
			for (int j = 0; j <= N - 1; j++)
			{
				F(i, j) = F_data[i*N+j];
			}
		}
	}
	//
	void updatePHIQk(Type dt)
	{
		updateF();
		// AA左上角[(0,0)-(N-1,N-1)]
		for (int i = 0; i <= N - 1; i++)
		{
			for (int j = 0; j <= N - 1; j++)
			{
				AA(i, j) = -F(i, j);
			}
		}
		// AA右上角[(0,N)-(N-1,2N-1)]
		for (int i = 0; i <= N - 1; i++)
		{
			for (int j = N; j <= 2 * N - 1; j++)
			{
				AA(i, j) = Qt(i,j-N);
			}
		}
		// AA左下角[(N,0)-(2N-1,N-1)]
		for (int i = N; i <= 2 * N - 1; i++)
		{
			for (int j = 0; j <= N - 1; j++)
			{
				AA(i, j) = 0.0f;
			}
		}
		// AA右下角[(N,N)-(2N-1,2N-1)]
		for (int i = N; i <= 2 * N - 1; i++)
		{
			for (int j = N; j <= 2 * N - 1; j++)
			{
				AA(i, j) = F(j - N, i - N);// (i-N)与(j-N)的位置调换，形成转置
			}
		}
		//
		BB = I_2nx2n + AA*dt; 
		//
		for (int i = 0; i <= N - 1; i++)
		{
			for (int j = 0; j <= N - 1; j++)
			{
				PHI(i, j) = BB(N+i,N+j);
			}
		}
		PHI = PHI.transpose();
		for (int i = 0; i <= N - 1; i++)
		{
			for (int j = 0; j <= N - 1; j++)
			{
				Qk(i, j) = BB(i, N + j);
			}
		}
		Qk = PHI * Qk;
	}

	void updateQt(Type * pQt_data, int size)
	{
		if (size == N)
		{
			// 所有元素置零
			Qt.zero();	
			// 设定对角线元素
			for (int i = 0; i <= N - 1; i++)
			{
				Qt(i, i) = pQt_data[i];
			}
		}
		else
		{
			while (1) {};
		}
	}

	void updateR_yaw(Type * pR_yaw_data, int size) 
	{
		if (size == M1)
		{
			// 所有元素置零
			R_yaw.zero();
			// 设定对角线元素
			for (int i = 0; i <= M1 - 1; i++)
			{
				R_yaw(i, i) = pR_yaw_data[i];
			}
		}
		else
		{
			while (1) {};
		}
	}

	void updateR_gps(Type* pR_gps_data, int size)
	{
		if (size == M2)
		{
			// 所有元素置零
			R_gps.zero();
			// 设定对角线元素
			for (int i = 0; i <= M2 - 1; i++)
			{
				R_gps(i, i) = pR_gps_data[i];
			}
		}
		else
		{
			while (1) {};
		}
	}

	// 一步预测&更新预测状态估计协方差
	void step_prediction(Type dt)
	{
		xhat = xhat + xdot * dt;
		P = PHI * P * PHI.transpose() + Qk;
	}
	// 航向角测量融合
	void fuse_yaw_measurement(void)
	{
		if (isMagYawVaild == 1)
		{
			//
			Type H_yaw_data[] = { 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };
			for (int i = 0; i <= N - 1; i++)
			{
				H_yaw(0, i) = H_yaw_data[i];
			}
			//
			K_yaw = P * H_yaw.transpose() * inv(static_cast<SquareMatrix<Type, M1>>(H_yaw*P*H_yaw.transpose()+ R_yaw));
			// 解缠绕
			while (Z_yaw(0, 0) > (xhat(2, 0) + PI)) { Z_yaw(0, 0) = Z_yaw(0, 0) - 2 * PI; }
			while (Z_yaw(0, 0) < (xhat(2, 0) - PI)) { Z_yaw(0, 0) = Z_yaw(0, 0) + 2 * PI; }
			xhat = xhat + K_yaw * (Z_yaw - H_yaw * xhat);
			//
			P = (I_nxn - K_yaw * H_yaw) * P;
		}
	}
	// GPS测量融合
	void fuse_gps_measurement(void)
	{
		if (isGPSvalid == 1)
		{
			Type H_gps_data[] = {
				 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,
				 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
				 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0,
				 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0
			};
			for (int i = 0; i <= M2 - 1; i++)
			{
				for (int j = 0; j <= N - 1; j++)
				{
					H_gps(i, j) = H_gps_data[i*N+j];
				}
			}
			K_gps = P * H_gps.transpose() * inv(static_cast<SquareMatrix<Type, M2>>(H_gps * P * H_gps.transpose() + R_gps));//static_cast<SquareMatrix<Type, M2>>
			xhat = xhat + K_gps * (Z_gps - H_gps * xhat);
			//
			P = (I_nxn - K_gps * H_gps) * P;
		}

	}
	// 滤波结果输出
	void output_result(Type * outcome, int outcome_dim)
	{
		if (outcome_dim == N)
		{
			for (int i = 0; i <= N - 1; i++)
			{
				outcome[i] = xhat(i, 0);
			}
		}
		else
		{
			while (1) {};
		}
	}
	
	void Run(
		Type * outcome,	int outcome_dim,			// 滤波结果
		Type dt,									// 时间间隔
		Type gx_rps, Type gy_rps, Type gz_rps,		// 陀螺仪
		Type fx_mps2, Type fy_mps2, Type fz_mps2,	// 加速度计 
		int isYawUpdate, Type yaw_rad,				// 航向角测量值
		int isGPSUpdate, Type Pn_m, Type Pe_m, Type h_msl_m, Type Vn_mps, Type Ve_mps, Type Vd_mps, // GPS测量值
		Type * pQt_data, int Qt_dim,				// Q 模型噪声描述
		Type * pR_yaw_data, int Ryaw_dim,			// R_yaw 航向角测量方差描述
		Type * pR_gps_data, int Rgps_dim			// R_gps 位置速度测量方差描述
	)
	{
		static int runOnce = 0;
		if (runOnce == 0)
		{
			init();
			runOnce = 1;
		}
		// Q R 参数调谐
		updateQt(pQt_data, Qt_dim);

		updateR_yaw(pR_yaw_data, Ryaw_dim);

		updateR_gps(pR_gps_data, Rgps_dim);
		// IMU 数据更新
		imu_measure(gx_rps, gy_rps, gz_rps, fx_mps2, fy_mps2, fz_mps2);
		// 航向角测量更新
		yaw_measure(isYawUpdate, yaw_rad);
		// GPS(三轴位置、速度更新)
		gps_measure(isGPSUpdate, Pn_m, Pe_m, h_msl_m, Vn_mps, Ve_mps, Vd_mps);
		// 更新状态向量的导数
		updateXdot();
		// 计算一步预测矩阵PHI & 离散化的过程噪声矩阵Qk
		updatePHIQk(dt);
		// 一步预测 & 更新协方差
		step_prediction(dt);
		// 磁强计测量融合
		fuse_yaw_measurement();
		// GPS 测量更新
		fuse_gps_measurement();
		// 输出
		output_result(outcome, outcome_dim);
	}
};
