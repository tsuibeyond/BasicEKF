clc
clear

generate_uav_sensors_structure;
% ���ļ���׼��д������
fid = fopen('UavData_SquareSpiral.txt','wt');
if(fid<0)   % �����ļ�ʧ�ܣ�ֱ���˳�����
    disp('Cannot open UavData_SquareSpiral.txt ������');
    disp('Exit :( !!!');
    return;             
else 
    disp('successful ����Open UavData_SquareSpiral.txt ������');
    disp('Ready for writing :) > > > > >');
    disp('Wait until echo Finish .....');
end
% д������
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
    % ������%lf������д��Ϊ��
end

fclose(fid);
disp('Finish :)')



%% ����
% 1. fid���ڴ洢�ļ����ֵ�����fid>0����˵���ļ��򿪳ɹ�
% 2. �򿪷�ʽ���£�
%     -��r����ֻ����ʽ���ļ���Ĭ�ϵķ�ʽ�������ļ������Ѵ��ڡ�?
%     -��r+������д��ʽ���ļ����򿪺��ȶ���д�����ļ������Ѵ��ڡ�?
%     -��w�����򿪺�д�����ݡ����ļ��Ѵ�������£��������򴴽���?
%     -��w+������д��ʽ���ļ����ȶ���д�����ļ��Ѵ�������£��������򴴽���?
%     -��a�����ڴ򿪵��ļ�ĩ��������ݡ��ļ��������򴴽���?
%     -��a+�������ļ����ȶ���������������ݡ��ļ��������򴴽���?
%     - ����Щ�ַ��������һ����t�����确rt����wt+�����򽫸��ļ����ı���ʽ�򿪣������ӵ��ǡ�b�������Զ����Ƹ�ʽ�򿪣���Ҳ��fopen����Ĭ�ϵĴ򿪷�ʽ��
% ע�����á�a��ʱ������ı����Ѿ��������ݣ�����������ݣ�����������֮��д�룬����w�������ԭ�������ݣ�����д��