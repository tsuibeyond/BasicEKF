%% ������ʵֵ����
load('UavData_SquareSpiral.mat');

%% ���崫������������
% ��������
gravity_mps2 = 9.81;

% ��ʼ��uavSensors�ṹ�壬����ֵʱ��
uavSensors=[];
uavSensors.sigmas=[];
uavSensors.biases=[];
uavSensors.params=[];
uavSensors.time_s = uavTruth.time_s;

% ���ô�����������������
uavSensors.sigmas.gyro_noise_rps            = 0.01;     % �����ǲ�����˹�������ԣ�rad/s
uavSensors.sigmas.gyro_bias_rps             = 0.01;     % �����ǲ�����ֵƫ�����ԣ�rad/s
uavSensors.sigmas.accel_noise_mps2          = 0.1;      % ���ٶȼƲ�����˹�������ԣ�m/s^2
uavSensors.sigmas.accel_bias_mps2           = 0.1;      % ���ٶȼƲ�����ֵƫ�����ԣ�m/s^2
uavSensors.sigmas.GPSpos_noise_m            = 2;        % GPSλ�ò�����˹�������ԣ�m
uavSensors.sigmas.GPSpos_bias_m             = 0;        % GPSλ�ò�����ֵƫ�m
uavSensors.sigmas.GPSvel_noise_mps          = 1;        % GPS�ٶȲ�����˹��������, m/s
uavSensors.sigmas.GPSvel_bias_mps           = 0;        % GPS�ٶȲ�����ֵƫ�m/s
uavSensors.sigmas.mag3D_unitVector_noise    = 0.02;     % ��ǿ��3D�����������ԣ���λ������
uavSensors.sigmas.mag3D_unitVector_bias     = 0;     % ��ǿ��3D������ֵƫ�����ԣ���λ������

% ���� GPS ����Ƶ��
uavSensors.params.dt_GPS_s = 1.0;  % GPS����ʱ������s
% �˴����ǽ���ƫ������
uavSensors.params.mag_declination_deg = 0*randn;    

%% ����������ֵ = ��ʵֵ + ��ֵƫ����*randn + ��˹������*randn
%% ����GPS��������
uavSensors.biases.GPS_east_m    = uavSensors.sigmas.GPSpos_bias_m*randn;    % Matlab�У�rand��0-1�ľ�������ֲ�����randn�Ǿ�ֵΪ0����Ϊ1����̬�ֲ�
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
% - GPS_valid:1-��Ч��0-��Ч  
% - GPS_valid=0ʱ��GPS����ֵΪNaN
kGPS = interp1(uavTruth.time_s,1:length(uavTruth.time_s),uavTruth.time_s(1):uavSensors.params.dt_GPS_s:uavTruth.time_s(end),'nearest'); % interp1(x,v,xq,method) ʹ�����Բ�ֵ����һά�������ض���ѯ��Ĳ���ֵ������ x ���������㣬v ������Ӧֵ v(x)������ xq ������ѯ�������
uavSensors.GPS_valid = false(size(uavSensors.time_s)); % \  Initialze GPS_valid to a vector of "false" values.
uavSensors.GPS_valid(kGPS)=true;                       % /  Set GPS_valid to "true" for GPS update times.
uavSensors.GPS_east_m(~uavSensors.GPS_valid)      = NaN; % \ 
uavSensors.GPS_north_m(~uavSensors.GPS_valid)     = NaN; %  | GPS values are set to NaN
uavSensors.GPS_h_msl_m(~uavSensors.GPS_valid)     = NaN; %  |   if not valid
uavSensors.GPS_v_ned_mps(~uavSensors.GPS_valid,:) = NaN; % /

% IMU �����ǽ��ٶȲ���ֵ����������ϵ�£�
uavSensors.biases.gyro_wb_rps = uavSensors.sigmas.gyro_bias_rps*randn(1,3);
for m=[1 2 3] % 1: body-x, 2: body-y, 3: body-z                        
    uavSensors.gyro_wb_rps(:,m) = uavTruth.wb_rps(:,m) ...
                            + uavSensors.biases.gyro_wb_rps(m) ...
                            + uavSensors.sigmas.gyro_noise_rps*randn(size(uavTruth.time_s));
end

% IMU ���ٶȼƲ���ֵ����������ϵ�£�
uavSensors.biases.accel_fb_mps2 = uavSensors.sigmas.accel_bias_mps2*randn(1,3);
% �ݶȷ��󵼺�ϵ��NED���µļ��ٶ�
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
    % ��ת���󣺵���ϵ��NED������������ϵ
    C_ned2b  = [cos(pitch_rad)*cos(yaw_rad)                                             cos(pitch_rad)*sin(yaw_rad)                                           -sin(pitch_rad); ...
                sin(roll_rad)*sin(pitch_rad)*cos(yaw_rad)-cos(roll_rad)*sin(yaw_rad)    sin(roll_rad)*sin(pitch_rad)*sin(yaw_rad)+cos(roll_rad)*cos(yaw_rad)   sin(roll_rad)*cos(pitch_rad); ...
                cos(roll_rad)*sin(pitch_rad)*cos(yaw_rad)+sin(roll_rad)*sin(yaw_rad)    cos(roll_rad)*sin(pitch_rad)*sin(yaw_rad)-sin(roll_rad)*cos(yaw_rad)   cos(roll_rad)*cos(pitch_rad)];

    % ������������ϵ�µļ��ٶ���ʵֵ
    fb_mps2_perfect = C_ned2b*vdot_ned_mps2(kTime,:)' - C_ned2b*[0;0;gravity_mps2];
    
    % ������������ϵ�µ�������ٶȲ���ֵ������������
    uavSensors.accel_fb_mps2(kTime,:) = fb_mps2_perfect' ...
                            + uavSensors.biases.accel_fb_mps2 ...
                            + uavSensors.sigmas.accel_noise_mps2*randn(1,3);    
end

% �����ǿ�Ʋ����õ���ƫ����
uavSensors.biases.mag3D_unitVector = uavSensors.sigmas.mag3D_unitVector_bias*randn(1,3);
for kTime=1:length(uavTruth.time_s)
    yaw_rad   = pi/180*uavTruth.yaw_deg(kTime);
    pitch_rad = pi/180*uavTruth.pitch_deg(kTime);
    roll_rad  = pi/180*uavTruth.roll_deg(kTime);
    % ����ٶȼƲ���ֵ�Ĵ���ʽ���ƣ��ȼ�����ת����
    C_ned2b  = [cos(pitch_rad)*cos(yaw_rad)                                             cos(pitch_rad)*sin(yaw_rad)                                           -sin(pitch_rad); ...
                sin(roll_rad)*sin(pitch_rad)*cos(yaw_rad)-cos(roll_rad)*sin(yaw_rad)    sin(roll_rad)*sin(pitch_rad)*sin(yaw_rad)+cos(roll_rad)*cos(yaw_rad)   sin(roll_rad)*cos(pitch_rad); ...
                cos(roll_rad)*sin(pitch_rad)*cos(yaw_rad)+sin(roll_rad)*sin(yaw_rad)    cos(roll_rad)*sin(pitch_rad)*sin(yaw_rad)-sin(roll_rad)*cos(yaw_rad)   cos(roll_rad)*cos(pitch_rad)];
    % ���Ǵ�ƫ�ǵ�Ӱ��
    C_mag2ned = [ cos(-pi/180*uavSensors.params.mag_declination_deg)  sin(-pi/180*uavSensors.params.mag_declination_deg)   0; ...
                 -sin(-pi/180*uavSensors.params.mag_declination_deg)  cos(-pi/180*uavSensors.params.mag_declination_deg)   0; ...
                 0                                                    0                                                    1];
    % [1;0;0] ��ʾ
    % ��-��-�أ���������ϵ�뱱-��-������ϵ����ʱ�������ǡ�����ǡ�ƫ����Ϊ�㣩����ʱ�����ǿ�Ƶĵ�X��������Y���Z�����ֵΪ0����һ��֮��Ϊ[1;0;0]���˴����Դ���ǵ�Ӱ�죬��ʵ�ʹ����д����Ӱ��ܴ�
    M = C_ned2b*C_mag2ned*[1;0;0] ...
              + uavSensors.biases.mag3D_unitVector' ...
              + uavSensors.sigmas.mag3D_unitVector_noise*randn(1,3)';
    M = M/norm(M);
    
    % ����õ������ǿ�ƵĲ���ֵ
    uavSensors.mag3D_unitVector_in_body(kTime,:) = M';
    
    % Ҳ���Ը���һ���ļ��裺�����ǿ��ֱ����������
    uavSensors.mag2D_yaw_deg(kTime,1) = 180/pi*atan2(-M(2),M(1)) + uavSensors.params.mag_declination_deg;
    uavSensors.mag2D_yaw_deg(kTime,1) = mod(uavSensors.mag2D_yaw_deg(kTime,1)+180,360)-180; % -180 <= yaw <=180 
    
    % ��ֱ�ӵõ�����ǣ��������ǿ�ƵĲ�������ֱ�ӵ�ЧΪ����ǵĲ������
    if kTime==1
        uavSensors.sigmas.mag2D_yaw_noise_rad = 1*uavSensors.sigmas.mag3D_unitVector_noise;
        uavSensors.sigmas.mag2D_yaw_bias_rad  = 1*uavSensors.sigmas.mag3D_unitVector_bias;
    end
end

