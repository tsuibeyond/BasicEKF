clc
clear 

% ���ɷ�������
generate_uav_sensors_structure;

% �������ٶȳ�����m/s^2
g_mps2 = 9.81;  

% ����directEKF�㷨��ر���
ekf = [];
% ��ǿ�ƺ�GPS�����������ƣ����ڹ��� R ����
ekf.sigmas.mag2D_meas_rad        = sqrt(uavSensors.sigmas.mag2D_yaw_noise_rad^2);
ekf.sigmas.pos_meas_m            = sqrt(uavSensors.sigmas.GPSpos_noise_m^2);
ekf.sigmas.vel_meas_mps          = sqrt(uavSensors.sigmas.GPSvel_noise_mps^2);
ekf.R_gps = diag([...
        [1 1 1]*ekf.sigmas.pos_meas_m ...            % North-East-Alt position measurement noise, m
        [1 1 1]*ekf.sigmas.vel_meas_mps ...          % NED velocity measurement noise, m/s
    ].^2);
ekf.R_gps = max(ekf.R_gps,(1e-3)^2*eye(size(ekf.R_gps)));    
ekf.R_mag2d = ekf.sigmas.mag2D_meas_rad;
ekf.R_mag2d = max(ekf.R_mag2d,(1e-3)^2*eye(size(ekf.R_mag2d)));    

% �����Ǻͼ��ٶȼƳ�ʼ���ƵĲ�ȷ���ȣ����ڹ���P����Э�������
ekf.sigmas.gyro_bias_rps         = uavSensors.sigmas.gyro_bias_rps;
ekf.sigmas.accel_bias_mps2       = uavSensors.sigmas.accel_bias_mps2;

% EKF�����������ƣ����ڹ��� Q ����
% Q ������Ҫ���ݾ�����ڡ�Q ֵ���������Ƕ�״̬ģ�͵����ų̶�
ekf.sigmas.attitude_process_noise_rad = 0.002;%0.002;% Euler angle process noise, rad
ekf.sigmas.pos_process_noise_m = 0.005;%0;           % Position process noise, m
ekf.sigmas.vel_process_noise_mps = 0.001;%2;         % Velocity process noise, m/s
ekf.sigmas.gyroBias_process_noise_rps = 1e-5;%1e-6; % Gyro bias process noise, rad/s
ekf.sigmas.accelBias_process_noise_mps2=1e-5;%1e-6; % Accel bias process noise, m/s^2
ekf.Q = diag([...
            [1 1 1]*ekf.sigmas.attitude_process_noise_rad ...   % Euler angle process noise, rad
            [1 1 1]*ekf.sigmas.pos_process_noise_m ...          % North-East-Alt position process noise, m
            [1 1 1]*ekf.sigmas.vel_process_noise_mps ...        % NED velocity process noise, m/s
            [1 1 1]*ekf.sigmas.gyroBias_process_noise_rps ...   % XYZ gyro bias process noise, rad/s
            [1 1 1]*ekf.sigmas.accelBias_process_noise_mps2 ... % XYZ accel bias process noise, m/s^2
        ].^2);
ekf.Q = max(ekf.Q,(1e-3)^2*eye(size(ekf.Q)));    

% �趨EKF״̬������ʼֵ
phi = 0;                                   % ����ǳ�ʼֵ��rad
theta=0;                                   % �����ǳ�ʼֵ��rad
psi = uavSensors.mag2D_yaw_deg(1)*pi/180;  % ƫ���ǳ�ʼֵ, rad
euler_init = [ phi theta psi ];   
pos_init = [ uavSensors.GPS_north_m(1) uavSensors.GPS_east_m(1) uavSensors.GPS_h_msl_m(1) ];
vel_init = uavSensors.GPS_v_ned_mps(1,:);
ekf.xhat = [ ...
        euler_init ...                              %   1-3: �� Euler angle (Roll, Pitch, Yaw) ״̬���г�ʼ����rad
        pos_init ...                                %   4-6: �� North-East-Alt λ��״̬���г�ʼ��, m
        vel_init ...                                %   7-9: �� NED �ٶ�״̬��ʼ��, m/s
        [0 0 0] ...                                 % 10-12: �� XYZ ���������ǵĳ�ֵƫ����г�ʼ��, rad/s
        [0 0 0] ...                                 % 13-15: �� XYZ ������ٶȼƵĳ���ƫ����г�ʼ��, m/s^2
    ]';
clear psi theta phi euler_init pos_init vel_init

% �Գ�ʼ��Э�������P�����趨�ϴ�ֵ���Ա�EKF�㷨�ڳ�ʼ�ļ����˲������и�������Ų�������
large_angle_uncertainty_rad = 30*pi/180;
large_pos_uncertainty_m = 100;
large_vel_uncertainty_mps = 10;
ekf.P = diag([...
        [1 1 1]*large_angle_uncertainty_rad ...     % init Euler angle (NED-to-body) uncertainties, rad
        [1 1 1]*large_pos_uncertainty_m ...         % init North-East-Alt position uncertainties, m
        [1 1 1]*large_vel_uncertainty_mps ...       % init NED velocity uncertainties, m/s
        [1 1 1]*ekf.sigmas.gyro_bias_rps ...        % init XYZ gyro bias uncertainties, rad/s
        [1 1 1]*ekf.sigmas.accel_bias_mps2 ...      % init XYZ accel bias uncertainties, m/s^2
    ].^2);
ekf.P = max(ekf.P,(1e-3)^2*eye(size(ekf.P)));    % Ϊ�˱���gyro_bias_rps��accel_bias_mps2�趨Ϊ0
clear large_angle_uncertainty_rad large_pos_uncertainty_m large_vel_uncertainty_mps

% ����uavEst�ṹ�岢�����ڴ�ռ䣬���ڴ洢���Ƶ��м���
uavEst.xhat = zeros(length(uavTruth.time_s),length(ekf.xhat));
uavEst.P    = zeros(length(uavTruth.time_s),length(ekf.xhat));

for kTime = 1 : length(uavSensors.time_s)
    % ���������˲���ʱ����
    time_s                  = uavSensors.time_s(kTime);
    if kTime == 1
        dt_s = uavSensors.time_s(2) - uavSensors.time_s(1);
    else
        dt_s = time_s - uavSensors.time_s(kTime-1);
    end
    % ��ȡ�˲��������Ƶ�״̬��
    phi=ekf.xhat(1);    theta=ekf.xhat(2);  psi=ekf.xhat(3);	% Roll, Pitch, Yaw Euler angles, rad
    Pn=ekf.xhat(4); Pe=ekf.xhat(5); Alt=ekf.xhat(6);            % Position, North/East/Altitude, meters
    Vn=ekf.xhat(7); Ve=ekf.xhat(8); Vd=ekf.xhat(9);             % Velocity, North/East/Down, m/s
    bwx=ekf.xhat(10);   bwy=ekf.xhat(11);   bwz=ekf.xhat(12);	% Gyro biases, rad/s
    bax=ekf.xhat(13);   bay=ekf.xhat(14);   baz=ekf.xhat(15);	% Accelerometer biases, m/s^2
    % IMU ����
    gyro_wb_rps             = uavSensors.gyro_wb_rps(kTime,:)';                     % 3x1 vector
    accel_fb_mps2           = uavSensors.accel_fb_mps2(kTime,:)';                   % 3x1 vector
    wx = gyro_wb_rps(1);    wy = gyro_wb_rps(2);    wz = gyro_wb_rps(3);            % rad/s
    fx = accel_fb_mps2(1);    fy = accel_fb_mps2(2);    fz = accel_fb_mps2(3);      % m/s^2
    % ���� xdot
    C_bodyrate2eulerdot  = [1      sin(phi)*tan(theta)    cos(phi)*tan(theta); ...
                            0           cos(phi)                -sin(phi)    ; ...
                            0      sin(phi)*sec(theta)    cos(phi)*sec(theta)];
    C_ned2b  = [cos(theta)*cos(psi)                               cos(theta)*sin(psi)                             -sin(theta); ...
                sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi)    sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi)   sin(phi)*cos(theta); ...
                cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi)    cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi)   cos(phi)*cos(theta)];
    C_b2ned=transpose(C_ned2b);
    xdot = [ ...
            C_bodyrate2eulerdot*([wx;wy;wz]-[bwx;bwy;bwz]); ...  % Derivative of [roll; pitch; yaw]
            [Vn; Ve; -Vd]; ...                                   % Derivative of [Pn; Pe; Alt]
            C_b2ned*([fx;fy;fz]-[bax;bay;baz])+[0;0;g_mps2]; ... % Derivative of [Vn; Ve; Vd]
            [0;0;0]; ...                                         % Derivative of [bwx; bwy; bwz]
            [0;0;0]; ...                                         % Derivative of [bax; bay; baz]
           ];
    % ���� F ����
    F = [ ...
        [                                                                 sin(phi)*tan(theta)*(bwz - wz) - cos(phi)*tan(theta)*(bwy - wy),                                  - cos(phi)*(bwz - wz)*(tan(theta)^2 + 1) - sin(phi)*(bwy - wy)*(tan(theta)^2 + 1),                                                                                                                                                              0, 0, 0, 0, 0, 0,  0, -1, -sin(phi)*tan(theta), -cos(phi)*tan(theta),                    0,                                                  0,                                                  0]
        [                                                                                       cos(phi)*(bwz - wz) + sin(phi)*(bwy - wy),                                                                                                                  0,                                                                                                                                                              0, 0, 0, 0, 0, 0,  0,  0,            -cos(phi),             sin(phi),                    0,                                                  0,                                                  0]
        [                                                             (sin(phi)*(bwz - wz))/cos(theta) - (cos(phi)*(bwy - wy))/cos(theta),                    - (cos(phi)*sin(theta)*(bwz - wz))/cos(theta)^2 - (sin(phi)*sin(theta)*(bwy - wy))/cos(theta)^2,                                                                                                                                                              0, 0, 0, 0, 0, 0,  0,  0, -sin(phi)/cos(theta), -cos(phi)/cos(theta),                    0,                                                  0,                                                  0]
        [                                                                                                                               0,                                                                                                                  0,                                                                                                                                                              0, 0, 0, 0, 1, 0,  0,  0,                    0,                    0,                    0,                                                  0,                                                  0]
        [                                                                                                                               0,                                                                                                                  0,                                                                                                                                                              0, 0, 0, 0, 0, 1,  0,  0,                    0,                    0,                    0,                                                  0,                                                  0]
        [                                                                                                                               0,                                                                                                                  0,                                                                                                                                                              0, 0, 0, 0, 0, 0, -1,  0,                    0,                    0,                    0,                                                  0,                                                  0]
        [ - (sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta))*(bay - fy) - (cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta))*(baz - fz), cos(psi)*sin(theta)*(bax - fx) - cos(phi)*cos(psi)*cos(theta)*(baz - fz) - cos(psi)*cos(theta)*sin(phi)*(bay - fy), (cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))*(bay - fy) - (cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta))*(baz - fz) + cos(theta)*sin(psi)*(bax - fx), 0, 0, 0, 0, 0,  0,  0,                    0,                    0, -cos(psi)*cos(theta),   cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta), - sin(phi)*sin(psi) - cos(phi)*cos(psi)*sin(theta)]
        [   (cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta))*(bay - fy) + (cos(phi)*cos(psi) + sin(phi)*sin(psi)*sin(theta))*(baz - fz), sin(psi)*sin(theta)*(bax - fx) - cos(phi)*cos(theta)*sin(psi)*(baz - fz) - cos(theta)*sin(phi)*sin(psi)*(bay - fy), (cos(phi)*sin(psi) - cos(psi)*sin(phi)*sin(theta))*(bay - fy) - (sin(phi)*sin(psi) + cos(phi)*cos(psi)*sin(theta))*(baz - fz) - cos(psi)*cos(theta)*(bax - fx), 0, 0, 0, 0, 0,  0,  0,                    0,                    0, -cos(theta)*sin(psi), - cos(phi)*cos(psi) - sin(phi)*sin(psi)*sin(theta),   cos(psi)*sin(phi) - cos(phi)*sin(psi)*sin(theta)]
        [                                                                 cos(theta)*sin(phi)*(baz - fz) - cos(phi)*cos(theta)*(bay - fy),                            cos(theta)*(bax - fx) + cos(phi)*sin(theta)*(baz - fz) + sin(phi)*sin(theta)*(bay - fy),                                                                                                                                                              0, 0, 0, 0, 0, 0,  0,  0,                    0,                    0,           sin(theta),                               -cos(theta)*sin(phi),                               -cos(phi)*cos(theta)]
        [                                                                                                                               0,                                                                                                                  0,                                                                                                                                                              0, 0, 0, 0, 0, 0,  0,  0,                    0,                    0,                    0,                                                  0,                                                  0]
        [                                                                                                                               0,                                                                                                                  0,                                                                                                                                                              0, 0, 0, 0, 0, 0,  0,  0,                    0,                    0,                    0,                                                  0,                                                  0]
        [                                                                                                                               0,                                                                                                                  0,                                                                                                                                                              0, 0, 0, 0, 0, 0,  0,  0,                    0,                    0,                    0,                                                  0,                                                  0]
        [                                                                                                                               0,                                                                                                                  0,                                                                                                                                                              0, 0, 0, 0, 0, 0,  0,  0,                    0,                    0,                    0,                                                  0,                                                  0]
        [                                                                                                                               0,                                                                                                                  0,                                                                                                                                                              0, 0, 0, 0, 0, 0,  0,  0,                    0,                    0,                    0,                                                  0,                                                  0]
        [                                                                                                                               0,                                                                                                                  0,                                                                                                                                                              0, 0, 0, 0, 0, 0,  0,  0,                    0,                    0,                    0,                                                  0,                                                  0]
    ];    
    % ���� PHI_k & Q_k ����һ��Ԥ����� �� ��ɢ��ʽ�Ĺ�����������
    nStates = length(ekf.xhat);
    AA = [-F ekf.Q; zeros(nStates,nStates) F']*dt_s;
    BB = expm(AA); % <- Matrix exponential!
    PHI_k = BB(nStates+1:2*nStates,nStates+1:2*nStates)';
    Q_k = PHI_k*BB(1:nStates,nStates+1:2*nStates);
    % һ��Ԥ�� & ����Э����
    ekf.xhat = ekf.xhat + xdot*dt_s;
    ekf.P = PHI_k*ekf.P*PHI_k' + Q_k;
    % ��ǿ�Ʋ����ںϣ�����Ƶ����IMUһ�£�
    mag2D_yaw_rad           = uavSensors.mag2D_yaw_deg(kTime)*pi/180;
    ekf.H_mag2d = [ 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
    ekf.K_mag2d = (ekf.P*ekf.H_mag2d')/(ekf.H_mag2d*ekf.P*ekf.H_mag2d'+ekf.R_mag2d);
    % ����yaw�ǵķ�ΧԼ��
    while mag2D_yaw_rad>ekf.xhat(3)+pi, mag2D_yaw_rad=mag2D_yaw_rad-2*pi; end
    while mag2D_yaw_rad<ekf.xhat(3)-pi, mag2D_yaw_rad=mag2D_yaw_rad+2*pi; end
    ekf.Z_mag2d = mag2D_yaw_rad;
    ekf.xhat = ekf.xhat + ekf.K_mag2d*(ekf.Z_mag2d - ekf.H_mag2d*ekf.xhat);
    ekf.P = (eye(length(ekf.xhat))-ekf.K_mag2d*ekf.H_mag2d)*ekf.P;
    % GPS �������� ������Ƶ�ʽϵͣ�
    GPS_east_m              = uavSensors.GPS_east_m(kTime);
    GPS_north_m             = uavSensors.GPS_north_m(kTime);
    GPS_h_msl_m             = uavSensors.GPS_h_msl_m(kTime);
    GPS_v_ned_mps           = uavSensors.GPS_v_ned_mps(kTime,:)';            % 3x1 vector
    GPS_valid               = uavSensors.GPS_valid(kTime);
    if GPS_valid
        ekf.H_gps = [ ...
            [ 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            [ 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            [ 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]
            [ 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0]
            [ 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0]
            [ 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0]
        ];
        ekf.K_gps = (ekf.P*ekf.H_gps')/(ekf.H_gps*ekf.P*ekf.H_gps'+ekf.R_gps);
        ekf.Z_gps = [ ...
            [GPS_north_m; GPS_east_m; GPS_h_msl_m];
            GPS_v_ned_mps
        ];
        ekf.xhat = ekf.xhat + ekf.K_gps*(ekf.Z_gps - ekf.H_gps*ekf.xhat);
        ekf.P = (eye(length(ekf.xhat))-ekf.K_gps*ekf.H_gps)*ekf.P;
    end
    % 
    uavEst.xhat(kTime,:) = ekf.xhat';
end

uavEst.time_s = uavSensors.time_s;
uavEst.north_m = uavEst.xhat(:,4);
uavEst.east_m = uavEst.xhat(:,5);
uavEst.h_msl_m = uavEst.xhat(:,6);
uavEst.v_ned_mps = uavEst.xhat(:,7:9);
uavEst.roll_deg = rad2deg(uavEst.xhat(:,1));
uavEst.pitch_deg = rad2deg(uavEst.xhat(:,2));
uavEst.yaw_deg = rad2deg(uavEst.xhat(:,3));

% ���ƽ�����ӻ�
figure(1)
clf;
% ��ͼ�����ͺ͵�������
truthLineType    = 'b-';
sensedLineType   = 'gd';
estimateLineType = 'r.';
markerSize       = 6;
linewidth        = 2;

% ��-�������λ�ù��ƽ��
subplot(4,2,[1 3])
plot(uavSensors.GPS_east_m, uavSensors.GPS_north_m, sensedLineType, ...
     uavEst.east_m,         uavEst.north_m,         estimateLineType, ...
     uavTruth.east_m,       uavTruth.north_m,       truthLineType, ...
     'markersize',markerSize, 'linewidth', linewidth);
hold on
plot(uavTruth.east_m(1),uavTruth.north_m(1),'co', ...
     'markersize',1.5*markerSize, 'linewidth', 2*linewidth);
hold off
axis equal
grid on
xlabel('East, m'); ylabel('North, m')
title('UAV Position')
legend('GPS Meas.','Est. Pos.','True Pos.','Start','location','best')

% ���θ߶ȹ��ƽ��
subplot(4,2,5)
plot(uavSensors.time_s, uavSensors.GPS_h_msl_m, sensedLineType, ...
     uavEst.time_s,     uavEst.h_msl_m,         estimateLineType, ...
     uavTruth.time_s,   uavTruth.h_msl_m,       truthLineType, ...
     'markersize',markerSize, 'linewidth', linewidth);
grid on
xlabel('Time, s'); ylabel('Altitude, above MeanSeaLevel, m')
legend('GPS Meas.','Est. Alt','True Alt','location','best')

% ��-��-�� �ٶȹ��ƽ��
subplot(4,2,7)
mag = @(v)(sqrt(sum(v.^2,2))); % magnitude of each row
plot(uavSensors.time_s, mag(uavSensors.GPS_v_ned_mps), sensedLineType, ...
     uavEst.time_s,     mag(uavEst.v_ned_mps),         estimateLineType, ...
     uavTruth.time_s,   mag(uavTruth.v_ned_mps),       truthLineType, ...
     'markersize',markerSize, 'linewidth', linewidth);
grid on
xlabel('Time, s'); ylabel('Inertial Speed, m/s')
legend('GPS Meas.','Est. Speed','True Speed','location','best')

% ����ǹ��ƽ��
subplot(3,2,2); 
plot(uavEst.time_s,     uavEst.roll_deg,    estimateLineType, ...
     uavTruth.time_s,   uavTruth.roll_deg,  truthLineType, ...    
     'markersize',markerSize, 'linewidth', linewidth);
grid on
xlabel('Time, s'); ylabel('Roll Angle, deg')
title('UAV Attitude: Roll')
legend('Est. Roll','True Roll','location','best')

% �����ǹ��ƽ��
subplot(3,2,4); 
plot(uavEst.time_s,     uavEst.pitch_deg,    estimateLineType, ...
     uavTruth.time_s,   uavTruth.pitch_deg,  truthLineType, ...     
     'markersize',markerSize, 'linewidth', linewidth);
grid on
xlabel('Time, s'); ylabel('Pitch Angle, deg')
title('UAV Attitude: Pitch')
legend('Est. Pitch','True Pitch','location','best')

% ƫ���ǹ��ƽ����-180,180��
subplot(3,2,6);
wrap = @(angleDeg)(mod(angleDeg+180,360)-180); % Wrap from -180deg to 180deg.
plot(uavSensors.time_s, wrap(uavSensors.mag2D_yaw_deg), sensedLineType, ...
     uavEst.time_s,     wrap(uavEst.yaw_deg),           estimateLineType, ...
     uavTruth.time_s,   wrap(uavTruth.yaw_deg),         truthLineType, ...     
     'markersize',markerSize, 'linewidth', linewidth);
grid on
xlabel('Time, s'); ylabel('Yaw Angle, deg')
title('UAV Attitude: Yaw')
legend('Magnetometer','Est. Yaw','True Yaw','location','best')
