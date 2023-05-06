function [qis, utc0] = viosimu(avp, ns, mubs, utc0)
%   vio 姿态四元数 'qis' 模拟器
%
% Prototype: [qis, utc0] = cnssimu(avp, datt, mubs, utc0)
% Inputs: avp - avp info, always from trajectory simulation- avp 信息，总是来自轨迹模拟
%         ns - VIO attitude observation noise              - VIO 姿态观测噪声
%                                                           ns=rk2 = [[10;10;30]*glv.sec];
%                                                           glv.sec = glv.min/60; % arcsec ——`单位:弧秒`
%         mubs - install angles to form DCM 'Cbs'          - 安装角度以形成 DCM 'Cbs' 
%                                                          -
%         utc0 - = [year,month,day,second,dUT1,dTAI]
% Outputs: qis - = [q1, q2, q3, t] array
%         utc0 - UTC info
%          
% Example:
%   [qis, utc0] = cnssimu(avp, [10;10;30]*glv.sec, [10;10;30]*glv.min, [2021;11;22;12*3600; -0.1;37]);
%
% See also  gpssimu, trjsimu, odsimu, bhsimu.

% Copyright(c) 2009-2021, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 22/11/2021
global glv
    if nargin<4,  utc0=[2021;11;22;0.0;-0.1;37];   end
    if nargin<3,  mubs=[0;0;0];   end
    if nargin<2,  ns=[0;0;0];   end

    Cie = cnsCie(utc0(1:3), utc0(4), utc0(5), utc0(6));
    Cbs = a2mat(mubs);           %将欧拉角转换为方向余弦矩阵（DCM）
    qis = avp(:,[1:3,end]);
    len = length(avp);   %len=13820;

%     timebar(1, len, 'VIO quaternion(qis) Simulation.');
    for k=1:len

%          视觉导航加入20倍误差
%         if k>=25000 && k<=40000
%         Cniose = a2mat(2000*ns.*randn(3,1));
%         else
%         Cniose = a2mat(ns.*randn(3,1)); 
%         end

% %         %不加入误差
      Cniose = a2mat(ns.*randn(3,1));
        qq = m2qua(Cie*rxyz(glv.wie*avp(k,end))*pos2cen(avp(k,7:9)')*a2mat(avp(k,1:3)')*Cbs*Cniose);  
        %m2qua 将方向余弦 转化为四元数    %Cie为i-e的旋转矩阵   %rxyz为矩阵 glv.wie：地球自转角速率  avp(k,end)：对应k时刻的时间                                 
                            %pos2cen：为Ce-n从e-n的旋转矩阵   %a2mat：输出Cn-b
                            %%Cb-s                           %Cniose
        if qq(1)>=0, 
            qis(k,1:3) = qq(2:4)';
        else, 
            qis(k,1:3) = -qq(2:4)';
        end  % make sure q(1)>=0

%         timebar;
    end
% 
%     myfig, plot(qis(:,end), qis(:,1:3)); xygo('q ^i_s (2,3,4)');
