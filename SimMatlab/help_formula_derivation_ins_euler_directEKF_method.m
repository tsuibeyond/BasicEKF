clc
clear

% ������ű���
syms phi theta psi                  % ����ǣ������ǣ�ƫ����
syms Pn Pe Alt                      % ����λ�ã�����λ�ã��߶�λ��
syms Vn Ve Vd                       % �����ٶȣ������ٶȣ��߶ȷ����ٶ�
syms wx wy wz                       % ������ٶ�
syms fx fy fz                       % ������ٶ�
syms bwx bwy bwz                    % ���߹��ƵĽ��ٶȳ�ֵƫ��
syms bax bay baz                    % ���߹��Ƶļ��ٶȳ�ֵƫ��
syms g_mps2                         % �������ٶȳ���

% ״̬����
X = [phi; theta; psi; Pn; Pe; Alt; Vn; Ve; Vd; bwx; bwy; bwz; bax; bay; baz];
%% ״̬����
% ��Ԫ��΢�ַ���
C_bodyrate2euldot  = [1      sin(phi)*tan(theta)    cos(phi)*tan(theta); ...
                      0           cos(phi)                -sin(phi)    ; ...
                      0      sin(phi)*sec(theta)    cos(phi)*sec(theta)];
% ��ת����
C_ned2b  = [cos(theta)*cos(psi)                               cos(theta)*sin(psi)                             -sin(theta); ...
            sin(phi)*sin(theta)*cos(psi)-cos(phi)*sin(psi)    sin(phi)*sin(theta)*sin(psi)+cos(phi)*cos(psi)   sin(phi)*cos(theta); ...
            cos(phi)*sin(theta)*cos(psi)+sin(phi)*sin(psi)    cos(phi)*sin(theta)*sin(psi)-sin(phi)*cos(psi)   cos(phi)*cos(theta)];
C_b2ned=transpose(C_ned2b);
% ״̬΢�ַ���
Xdot = [ ...
        C_bodyrate2euldot*([wx;wy;wz]-[bwx;bwy;bwz]); ...      % ��Ԫ��΢�ַ���
        [Vn; Ve; -Vd]; ...
        C_b2ned*([fx;fy;fz]-[bax;bay;baz])+[0;0;g_mps2]; ... % �ٶ�΢�ַ���
        [0;0;0]; ...                                         % �����ǳ�ֵƫ��΢�ַ���
        [0;0;0]; ...                                         % ���ٶȼƳ�ֵƫ��΢�ַ���
       ];
for n = 1:length(X)
    F(:,n) = diff(Xdot,X(n));
end
disp('F����')
disp(F)
%% ���ⷽ��
Zhat1 = psi;
Zhat2 = [ ...
    [Pn; Pe; Alt]; ...
    [Vn; Ve; Vd]; ...             % GPS�ٶȹ۲� ��tip:$PUBX,������
    ];
Zhat = [Zhat1;Zhat2];
for n = 1:length(X)
   H1(:,n) = diff(Zhat1,X(n));
   H2(:,n) = diff(Zhat2,X(n));
end
disp('H1����(MAG2D)')
disp(H1)
disp('H2����(GPS)')
disp(H2)

