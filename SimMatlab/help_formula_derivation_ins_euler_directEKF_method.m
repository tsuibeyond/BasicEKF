clc
clear

% 定义符号变量
syms phi theta psi                  % 横滚角，俯仰角，偏航角
syms Pn Pe Alt                      % 北向位置，东向位置，高度位置
syms Vn Ve Vd                       % 北向速度，东向速度，高度方向速度
syms wx wy wz                       % 三轴角速度
syms fx fy fz                       % 三轴加速度
syms bwx bwy bwz                    % 在线估计的角速度常值偏差
syms bax bay baz                    % 在线估计的加速度常值偏差
syms g_mps2                         % 重力加速度常量

% 状态向量
X = [phi; theta; psi; Pn; Pe; Alt; Vn; Ve; Vd; bwx; bwy; bwz; bax; bay; baz];
%% 状态方程
% 四元数微分方程
C_bodyrate2euldot  = [1      sin(phi)*tan(theta)    cos(phi)*tan(theta); ...
                      0           cos(phi)                -sin(phi)    ; ...
                      0      sin(phi)*sec(theta)    cos(phi)*sec(theta)];
% 旋转矩阵
C_ned2b  = [cos(theta)*cos(psi)                               cos(theta)*sin(psi)                             -sin(theta); ...
            sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi)    sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi)   sin(phi)*cos(theta); ...
            cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi)    cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi)   cos(phi)*cos(theta)];
C_b2ned=transpose(C_ned2b);
% 状态微分方程
Xdot = [ ...
        C_bodyrate2euldot*([wx;wy;wz]-[bwx;bwy;bwz]); ...      % 四元数微分方程
        [Vn; Ve; -Vd]; ...
        C_b2ned*([fx;fy;fz]-[bax;bay;baz])+[0;0;g_mps2]; ... % 速度微分方程
        [0;0;0]; ...                                         % 陀螺仪常值偏差微分方程
        [0;0;0]; ...                                         % 加速度计常值偏差微分方程
       ];
for n = 1:length(X)
    F(:,n) = diff(Xdot,X(n));
end
disp('F矩阵')
disp(F)
%% 量测方程
Zhat1 = psi;
Zhat2 = [ ...
    [Pn; Pe; Alt]; ...
    [Vn; Ve; Vd]; ...             % GPS速度观测 （tip:$PUBX,……）
    ];
Zhat = [Zhat1;Zhat2];
for n = 1:length(X)
   H1(:,n) = diff(Zhat1,X(n));
   H2(:,n) = diff(Zhat2,X(n));
end
disp('H1矩阵(MAG2D)')
disp(H1)
disp('H2矩阵(GPS)')
disp(H2)

