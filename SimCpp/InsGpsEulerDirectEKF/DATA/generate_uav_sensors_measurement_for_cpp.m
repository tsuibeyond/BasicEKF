clc
clear

generate_uav_sensors_structure;
% 打开文件，准备写入数据
fid = fopen('UavData_SquareSpiral.txt','wt');
if(fid<0)   % 若打开文件失败，直接退出程序
    disp('Cannot open UavData_SquareSpiral.txt ！！！');
    disp('Exit :( !!!');
    return;             
else 
    disp('successful ！！Open UavData_SquareSpiral.txt ！！！');
    disp('Ready for writing :) > > > > >');
    disp('Wait until echo Finish .....');
end
% 写入数据
for i=1:1:length(uavSensors.time_s)
    time_s = uavSensors.time_s(i);
    gx_rps = uavSensors.gyro_wb_rps(i,1);
    gy_rps = uavSensors.gyro_wb_rps(i,2);
    gz_rps = uavSensors.gyro_wb_rps(i,3);
    fx_mps2 = uavSensors.accel_fb_mps2(i,1);
    fy_mps2 = uavSensors.accel_fb_mps2(i,2);
    fz_mps2 = uavSensors.accel_fb_mps2(i,3);
    isYawUpdate = 1;
    yaw_rad = uavSensors.mag2D_yaw_deg(i);
    isGPSUpdate = uavSensors.GPS_valid(i);
    if isGPSUpdate==0
        Pn_m = 0;
        Pe_m = 0;
        h_msl_m = 0;
        Vn_mps = 0;
        Ve_mps = 0;
        Vd_mps = 0;
    else
        Pn_m = uavSensors.GPS_north_m(i);
        Pe_m = uavSensors.GPS_east_m(i);
        h_msl_m = uavSensors.GPS_h_msl_m(i);
        Vn_mps = uavSensors.GPS_v_ned_mps(i,1);
        Ve_mps = uavSensors.GPS_v_ned_mps(i,2);
        Vd_mps = uavSensors.GPS_v_ned_mps(i,3);
    end
    fprintf(fid,'%f,%f,%f,%f,%f,%f,%f,%d,%f,%d,%f,%f,%f,%f,%f,%f\n',time_s,gx_rps,gy_rps,gz_rps,fx_mps2,fy_mps2,fz_mps2,isYawUpdate,yaw_rad,isGPSUpdate,Pn_m,Pe_m,h_msl_m,Vn_mps,Ve_mps,Vd_mps);
    % 不能用%lf，否则写入为空
end

fclose(fid);
disp('Finish :)')



%% 辅助
% 1. fid用于存储文件句柄值，如果fid>0，这说明文件打开成功
% 2. 打开方式如下：
%     -‘r’：只读方式打开文件（默认的方式），该文件必须已存在。?
%     -‘r+’：读写方式打开文件，打开后先读后写。该文件必须已存在。?
%     -‘w’：打开后写入数据。该文件已存在则更新；不存在则创建。?
%     -‘w+’：读写方式打开文件。先读后写。该文件已存在则更新；不存在则创建。?
%     -‘a’：在打开的文件末端添加数据。文件不存在则创建。?
%     -‘a+’：打开文件后，先读入数据再添加数据。文件不存在则创建。?
%     - 在这些字符串后添加一个“t”，如‘rt’或‘wt+’，则将该文件以文本方式打开；如果添加的是“b”，则以二进制格式打开，这也是fopen函数默认的打开方式。
% 注：当用‘a’时，如果文本中已经存在数据，不会清空数据，而是在数据之后写入，而‘w’会清空原本的数据，重新写入