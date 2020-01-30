clc
clear 

generate_uav_sensors_structure;

EstData = load('EstInsGpsEulerDirectEKF.txt');
uavEst.time_s = EstData(:,1);
uavEst.roll_deg = rad2deg(EstData(:,2));
uavEst.pitch_deg = rad2deg(EstData(:,3));
uavEst.yaw_deg = rad2deg(EstData(:,4));
uavEst.north_m = EstData(:,5);
uavEst.east_m = EstData(:,6);
uavEst.h_msl_m = EstData(:,7);
uavEst.v_ned_mps = EstData(:,8:10);


% 估计结果可视化
figure(1)
clf;
% 绘图的线型和点型设置
truthLineType    = 'b-';
sensedLineType   = 'gd';
estimateLineType = 'r.';
markerSize       = 6;
linewidth        = 2;

% 北-东方向的位置估计结果
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

% 海拔高度估计结果
subplot(4,2,5)
plot(uavSensors.time_s, uavSensors.GPS_h_msl_m, sensedLineType, ...
     uavEst.time_s,     uavEst.h_msl_m,         estimateLineType, ...
     uavTruth.time_s,   uavTruth.h_msl_m,       truthLineType, ...
     'markersize',markerSize, 'linewidth', linewidth);
grid on
xlabel('Time, s'); ylabel('Altitude, above MeanSeaLevel, m')
legend('GPS Meas.','Est. Alt','True Alt','location','best')

% 北-东-地 速度估计结果
subplot(4,2,7)
mag = @(v)(sqrt(sum(v.^2,2))); % magnitude of each row
plot(uavSensors.time_s, mag(uavSensors.GPS_v_ned_mps), sensedLineType, ...
     uavEst.time_s,     mag(uavEst.v_ned_mps),         estimateLineType, ...
     uavTruth.time_s,   mag(uavTruth.v_ned_mps),       truthLineType, ...
     'markersize',markerSize, 'linewidth', linewidth);
grid on
xlabel('Time, s'); ylabel('Inertial Speed, m/s')
legend('GPS Meas.','Est. Speed','True Speed','location','best')

% 横滚角估计结果
subplot(3,2,2); 
plot(uavEst.time_s,     uavEst.roll_deg,    estimateLineType, ...
     uavTruth.time_s,   uavTruth.roll_deg,  truthLineType, ...    
     'markersize',markerSize, 'linewidth', linewidth);
grid on
xlabel('Time, s'); ylabel('Roll Angle, deg')
title('UAV Attitude: Roll')
legend('Est. Roll','True Roll','location','best')

% 俯仰角估计结果
subplot(3,2,4); 
plot(uavEst.time_s,     uavEst.pitch_deg,    estimateLineType, ...
     uavTruth.time_s,   uavTruth.pitch_deg,  truthLineType, ...     
     'markersize',markerSize, 'linewidth', linewidth);
grid on
xlabel('Time, s'); ylabel('Pitch Angle, deg')
title('UAV Attitude: Pitch')
legend('Est. Pitch','True Pitch','location','best')

% 偏航角估计结果（-180,180）
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
