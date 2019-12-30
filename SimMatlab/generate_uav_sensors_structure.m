%% 导入真实值数据
load('UavData_SquareSpiral.mat');

%% 定义传感器测量特性
% 重力常量
gravity_mps2 = 9.81;

% 初始化uavSensors结构体，并赋值时间
uavSensors=[];
uavSensors.sigmas=[];
uavSensors.biases=[];
uavSensors.params=[];
uavSensors.time_s = uavTruth.time_s;

% 设置传感器测量噪声特性
uavSensors.sigmas.gyro_noise_rps            = 0.01;     % 陀螺仪测量高斯噪声特性，rad/s
uavSensors.sigmas.gyro_bias_rps             = 0.01;     % 陀螺仪测量常值偏差特性，rad/s
uavSensors.sigmas.accel_noise_mps2          = 0.1;      % 加速度计测量高斯噪声特性，m/s^2
uavSensors.sigmas.accel_bias_mps2           = 0.1;      % 加速度计测量常值偏差特性，m/s^2
uavSensors.sigmas.GPSpos_noise_m            = 2;        % GPS位置测量高斯噪声特性，m
uavSensors.sigmas.GPSpos_bias_m             = 0;        % GPS位置测量常值偏差，m
uavSensors.sigmas.GPSvel_noise_mps          = 1;        % GPS速度测量高斯噪声特性, m/s
uavSensors.sigmas.GPSvel_bias_mps           = 0;        % GPS速度测量常值偏差，m/s
uavSensors.sigmas.mag3D_unitVector_noise    = 0.02;     % 磁强计3D测量噪声特性（单位向量）
uavSensors.sigmas.mag3D_unitVector_bias     = 0;     % 磁强计3D测量常值偏差特性（单位向量）

% 设置 GPS 更新频率
uavSensors.params.dt_GPS_s = 1.0;  % GPS更新时间间隔，s
% 此处我们将磁偏角置零
uavSensors.params.mag_declination_deg = 0*randn;    

%% 传感器测量值 = 真实值 + 常值偏差项*randn + 高斯噪声项*randn
%% 生成GPS测量数据
uavSensors.biases.GPS_east_m    = uavSensors.sigmas.GPSpos_bias_m*randn;    % Matlab中，rand是0-1的均匀随机分布，而randn是均值为0方差为1的正态分布
uavSensors.biases.GPS_north_m   = uavSensors.sigmas.GPSpos_bias_m*randn;
uavSensors.biases.GPS_h_m       = uavSensors.sigmas.GPSpos_bias_m*randn;
uavSensors.biases.GPS_Vned_mps= uavSensors.sigmas.GPSvel_bias_mps*randn(1,3);

uavSensors.GPS_east_m = uavTruth.east_m ... 
                            + uavSensors.biases.GPS_east_m ...
                            + uavSensors.sigmas.GPSpos_noise_m*randn(size(uavTruth.time_s));
uavSensors.GPS_north_m = uavTruth.north_m ...
                            + uavSensors.biases.GPS_north_m ...
                            + uavSensors.sigmas.GPSpos_noise_m*randn(size(uavTruth.time_s));
uavSensors.GPS_h_msl_m = uavTruth.h_msl_m ... 
                            + uavSensors.biases.GPS_h_m ...
                            + uavSensors.sigmas.GPSpos_noise_m*randn(size(uavTruth.time_s));
for m=[1 2 3] % 1: North, 2: East, 3: Down                        
    uavSensors.GPS_v_ned_mps(:,m) = uavTruth.v_ned_mps(:,m) ...
                                + uavSensors.biases.GPS_Vned_mps(m) ...
                                + uavSensors.sigmas.GPSvel_noise_mps*randn(size(uavTruth.time_s));
end
% - GPS_valid:1-有效；0-无效  
% - GPS_valid=0时，GPS测量值为NaN
kGPS = interp1(uavTruth.time_s,1:length(uavTruth.time_s),uavTruth.time_s(1):uavSensors.params.dt_GPS_s:uavTruth.time_s(end),'nearest'); % interp1(x,v,xq,method) 使用线性插值返回一维函数在特定查询点的插入值。向量 x 包含样本点，v 包含对应值 v(x)。向量 xq 包含查询点的坐标
uavSensors.GPS_valid = false(size(uavSensors.time_s)); % \  Initialze GPS_valid to a vector of "false" values.
uavSensors.GPS_valid(kGPS)=true;                       % /  Set GPS_valid to "true" for GPS update times.
uavSensors.GPS_east_m(~uavSensors.GPS_valid)      = NaN; % \ 
uavSensors.GPS_north_m(~uavSensors.GPS_valid)     = NaN; %  | GPS values are set to NaN
uavSensors.GPS_h_msl_m(~uavSensors.GPS_valid)     = NaN; %  |   if not valid
uavSensors.GPS_v_ned_mps(~uavSensors.GPS_valid,:) = NaN; % /

% IMU 陀螺仪角速度测量值（机体坐标系下）
uavSensors.biases.gyro_wb_rps = uavSensors.sigmas.gyro_bias_rps*randn(1,3);
for m=[1 2 3] % 1: body-x, 2: body-y, 3: body-z                        
    uavSensors.gyro_wb_rps(:,m) = uavTruth.wb_rps(:,m) ...
                            + uavSensors.biases.gyro_wb_rps(m) ...
                            + uavSensors.sigmas.gyro_noise_rps*randn(size(uavTruth.time_s));
end

% IMU 加速度计测量值（机体坐标系下）
uavSensors.biases.accel_fb_mps2 = uavSensors.sigmas.accel_bias_mps2*randn(1,3);
% 梯度法求导航系（NED）下的加速度
vdot_n_mps2 = gradient(uavTruth.v_ned_mps(:,1),uavTruth.time_s);
vdot_e_mps2 = gradient(uavTruth.v_ned_mps(:,2),uavTruth.time_s);
vdot_d_mps2 = gradient(uavTruth.v_ned_mps(:,3),uavTruth.time_s);
vdot_ned_mps2 = [vdot_n_mps2 vdot_e_mps2 vdot_d_mps2]; % [nx3 double]
for kTime=1:length(uavTruth.time_s)
    % Compute NED-to-body Direction Cosine Matrix for time index kTime.
    % See rotation_examples.m for more information about converting between
    % different rotation representations.  Here we use the method in
    % Section B of rotation_examples.m to convert from Euler angles to DCM.
    yaw_rad   = pi/180*uavTruth.yaw_deg(kTime);
    pitch_rad = pi/180*uavTruth.pitch_deg(kTime);
    roll_rad  = pi/180*uavTruth.roll_deg(kTime);
    % 旋转矩阵：导航系（NED）到机体坐标系
    C_ned2b  = [cos(pitch_rad)*cos(yaw_rad)                                             cos(pitch_rad)*sin(yaw_rad)                                           -sin(pitch_rad); ...
                sin(roll_rad)*sin(pitch_rad)*cos(yaw_rad)-cos(roll_rad)*sin(yaw_rad)    sin(roll_rad)*sin(pitch_rad)*sin(yaw_rad)+cos(roll_rad)*cos(yaw_rad)   sin(roll_rad)*cos(pitch_rad); ...
                cos(roll_rad)*sin(pitch_rad)*cos(yaw_rad)+sin(roll_rad)*sin(yaw_rad)    cos(roll_rad)*sin(pitch_rad)*sin(yaw_rad)-sin(roll_rad)*cos(yaw_rad)   cos(roll_rad)*cos(pitch_rad)];

    % 计算载体坐标系下的加速度真实值
    fb_mps2_perfect = C_ned2b*vdot_ned_mps2(kTime,:)' - C_ned2b*[0;0;gravity_mps2];
    
    % 计算载体坐标系下的三轴加速度测量值（含有噪声）
    uavSensors.accel_fb_mps2(kTime,:) = fb_mps2_perfect' ...
                            + uavSensors.biases.accel_fb_mps2 ...
                            + uavSensors.sigmas.accel_noise_mps2*randn(1,3);    
end

% 三轴磁强计测量得到的偏航角
uavSensors.biases.mag3D_unitVector = uavSensors.sigmas.mag3D_unitVector_bias*randn(1,3);
for kTime=1:length(uavTruth.time_s)
    yaw_rad   = pi/180*uavTruth.yaw_deg(kTime);
    pitch_rad = pi/180*uavTruth.pitch_deg(kTime);
    roll_rad  = pi/180*uavTruth.roll_deg(kTime);
    % 与加速度计测量值的处理方式类似，先计算旋转矩阵
    C_ned2b  = [cos(pitch_rad)*cos(yaw_rad)                                             cos(pitch_rad)*sin(yaw_rad)                                           -sin(pitch_rad); ...
                sin(roll_rad)*sin(pitch_rad)*cos(yaw_rad)-cos(roll_rad)*sin(yaw_rad)    sin(roll_rad)*sin(pitch_rad)*sin(yaw_rad)+cos(roll_rad)*cos(yaw_rad)   sin(roll_rad)*cos(pitch_rad); ...
                cos(roll_rad)*sin(pitch_rad)*cos(yaw_rad)+sin(roll_rad)*sin(yaw_rad)    cos(roll_rad)*sin(pitch_rad)*sin(yaw_rad)-sin(roll_rad)*cos(yaw_rad)   cos(roll_rad)*cos(pitch_rad)];
    % 考虑磁偏角的影响
    C_mag2ned = [ cos(-pi/180*uavSensors.params.mag_declination_deg)  sin(-pi/180*uavSensors.params.mag_declination_deg)   0; ...
                 -sin(-pi/180*uavSensors.params.mag_declination_deg)  cos(-pi/180*uavSensors.params.mag_declination_deg)   0; ...
                 0                                                    0                                                    1];
    % [1;0;0] 表示
    % 北-东-地，载体坐标系与北-东-地坐标系对齐时（俯仰角、横滚角、偏航角为零），此时三轴磁强计的的X轴测量最大，Y轴和Z轴测量值为0，归一化之后为[1;0;0]（此处忽略磁倾角的影响，但实际工程中磁倾角影响很大）
    M = C_ned2b*C_mag2ned*[1;0;0] ...
              + uavSensors.biases.mag3D_unitVector' ...
              + uavSensors.sigmas.mag3D_unitVector_noise*randn(1,3)';
    M = M/norm(M);
    
    % 计算得到三轴磁强计的测量值
    uavSensors.mag3D_unitVector_in_body(kTime,:) = M';
    
    % 也可以更进一步的假设：三轴磁强计直接输出航向角
    uavSensors.mag2D_yaw_deg(kTime,1) = 180/pi*atan2(-M(2),M(1)) + uavSensors.params.mag_declination_deg;
    uavSensors.mag2D_yaw_deg(kTime,1) = mod(uavSensors.mag2D_yaw_deg(kTime,1)+180,360)-180; % -180 <= yaw <=180 
    
    % 若直接得到航向角，将三轴磁强计的测量噪声直接等效为航向角的测量误差
    if kTime==1
        uavSensors.sigmas.mag2D_yaw_noise_rad = 1*uavSensors.sigmas.mag3D_unitVector_noise;
        uavSensors.sigmas.mag2D_yaw_bias_rad  = 1*uavSensors.sigmas.mag3D_unitVector_bias;
    end
end

