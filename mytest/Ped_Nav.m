clc;
close all;
%% 轨迹仿真
glvs
ts = 0.01;       % sampling interval 采样间隔
glv.pos0 = [40.143045; 116.239318; 380];
avp0 = [[0;0;0]; [0;0;0]; glv.pos0]; % init avp 经度:116.239318纬度:40.143045
 
% trajectory segment setting
xxx = [];
seg = trjsegment(xxx, 'init',         0);
% seg = trjsegment(seg, 'uniform',      20);
seg = trjsegment(seg, 'accelerate',   2.5, xxx, 1);
seg = trjsegment(seg, 'uniform',      100);
seg = trjsegment(seg, 'headdown',        1, 10, xxx, 1);
seg = trjsegment(seg, 'uniform',      5);
seg = trjsegment(seg, 'headup',        1, 10, xxx, 1);
seg = trjsegment(seg, 'uniform',      100);


seg = trjsegment(seg, 'turnleft',   9, 10, xxx, 1);
seg = trjsegment(seg, 'uniform',      200);
seg = trjsegment(seg, 'turnleft',  9, 10, xxx, 1);
seg = trjsegment(seg, 'uniform',      200);
seg = trjsegment(seg, 'turnleft',    9, 10, xxx, 1);
seg = trjsegment(seg, 'uniform',      80);
seg = trjsegment(seg, 'turnleft',    9, 10, xxx, 1);
% seg = trjsegment(seg, 'deaccelerate', 1.1,  xxx, 1.1);   %如果加上就是完整的闭环
seg = trjsegment(seg, 'headup',        1, 10, xxx, 1);
seg = trjsegment(seg, 'uniform',      10);
seg = trjsegment(seg, 'headdown',        1, 10, xxx, 1);
seg = trjsegment(seg, 'uniform',      150);

% seg = trjsegment(seg, 'turnleft',    0.5, 180, xxx, 1);
% seg = trjsegment(seg, 'turnleft',    0.5, 180, xxx, 1);
% seg = trjsegment(seg, 'climb',        1, 20, xxx, 1);
% seg = trjsegment(seg, 'uniform',      1);

% generate, save & plot
trj = trjsimu(avp0, seg.wat, ts, 1);
trjfile('trj10ms.mat', trj);
insplot(trj.avp);
figure();
plot3(trj.avp(:,7),trj.avp(:,8),trj.avp(:,9),'k','linewidth',1);
xlabel('Latitude/°');ylabel('Longitude/°');zlabel('Altitude/m');

hold on
plot3(trj.avp(1,7),trj.avp(1,8),trj.avp(1,9),'kO');

hold on
plot3(trj.avp(end,7),trj.avp(end,8),trj.avp(end,9),'k*');

grid on
% imuplot(trj.imu);

%% add noise and simulator
[nn, ts, nts] = nnts(2, trj.ts);                %nn=2; ts=trj.ts; nts=nn*trj.ts;

% GPS simulator
lever = [1; 2; 3]*0; 
rk1 = vperrset(0.1, 1);

davp0 = avperrset([0.5;-0.5;20], 0.01, [1;1;3]);    %输入对应 平台未对准角度(单位弧分 1°的六十分之一) 速度误差（m/s） 位置误差(m)
                  %davp0 输出
                  %位置误差[1 1 3] 经纬度 变成弧度？？  高度单位还是m
gps = gpssimu(trj.avp, davp0(4:6), davp0(7:9), 1, lever, 0.0);    

%gpsplot(gps)


% VIO simulator
mu = 0*[1;2;5]*glv.min; 
rk2 = [[10;10;30]*glv.sec];                                   %加入的误差
[qis, utc0] = viosimu(trj.avp, rk2, mu, [2021;11;22;12*3600; -0.1;37]);

%qis前三列为 方向余弦矩阵转化为四元数的后三位数  %%最后一列为时间100hz 

% SINS init
imuerr = imuerrset(0.03, [100;100;100], 0.001, 5);
imu = imuadderr(trj.imu, imuerr);                                          % imuplot(imu)
ins = insinit(trj.avp0(1:9), ts, davp0);  
len = length(imu);
Cie0 = cnsCie(utc0(1:3), utc0(4),  utc0(5), utc0(6));                       %cnsCie坐标转换矩阵？？

%%  centralized KF setting
psinstypedef(156);
ckf = [];
ckf.Phikk_1 = eye(15);                                 %矩阵F
ckf.Qt = diag([imuerr.web; imuerr.wdb; zeros(9,1)])^2; % 15-state  web角度随机游走 wdb速度随机游走 
ckf.Rk = diag([rk1; rk2])^2;                           %量测噪声
ckf.Pxk = diag([davp0; imuerr.eb; imuerr.db]*1.0)^2;   %eb为陀螺常值零偏 db为加计常值零偏  误差协方差矩阵P
ckf.Hk = [ zeros(6,3), eye(6), zeros(6,6); ...       % SINS/GPS Hk
          eye(3), zeros(3,12) ];  %ckf.Hk(7,7)=1;     % SINS/VIO Hk

ckf = kfinit0(ckf, nts); 

% federated KF setting
fins = ins;
fkf = fkfinit(ckf, {1:15,1:15}, {1:6,7:9}, [1/24; 22/24]);

% alloc memory分配内存
[avpc, xkpkc] = prealloc(fix(len/nn), 10, 2*ckf.n+1);
[avp1, xkpk1] = prealloc(fix(len/nn), 10, 2*fkf{1}.n+1);
[avp2, xkpk2] = prealloc(fix(len/nn), 10, 2*fkf{2}.n+1);
[avpf, xkpkf] = prealloc(fix(len/nn), 10, 2*fkf{3}.n+1);
ki = timebar(nn, len, '15-state SINS/GPS/VIO centralized v.s. federated KF simulation.');
imugpssyn(imu(:,7), gps(:,end)); 

%% %% federated KF with 15-state

for k=1:nn:len-nn+1             %nn=2; len = length(imu) 
    k1 = k+nn-1;                %k1是从2 到 len

    wvm = imu(k:k1,1:6);  t = imu(k1,end);   % '角度增量和速度增量，采样时刻'

    ins = insupdate(ins, wvm);       % '惯导更新'
    fins = insupdate(fins, wvm);

    ckf.Phikk_1 = kffk(ins);   %'计算状态转移矩阵' ckf.phikk_1= F为x点=Fx+u 中的F 
    
    [kgps, dt] = imugpssyn(k, k1, 'F');

    if kgps>0 && mod(t,1)<1.5*ts && norm(ins.wnb)<1*glv.dps  % SINS/GPS/VIO &&与的逻辑符：只要有一个条件不满足就判断错误，不在判断后面是否正确
                                                             % glv.dps=pi/180;
        Cns = cnsCns(qis(k1,1:3)', ins.pos, Cie0, t);    %in CNS应用中，坐标转换矩阵 n to 载体系
        ckf.Hk(8:9,8) = -[ins.eth.cl; ins.eth.sl];

        zk = [ins.vn-gps(kgps,1:3)'; ins.pos-gps(kgps,4:6)'; qq2phi(ins.qnb,m2qua(Cns))]; %测量值 
        ckf = kfupdate(ckf, zk);  % centralized KF  量测更新
        
        Cns = cnsCns(qis(k1,1:3)', fins.pos, Cie0, t);
        ckf.Hk(8:9,8) = -[fins.eth.cl; fins.eth.sl];
        zk = [fins.vn-gps(kgps,1:3)'; fins.pos-gps(kgps,4:6)'; qq2phi(fins.qnb,m2qua(Cns))];
        fkf = fkfupdate(ckf, fkf, zk);  % federated KF
    else
        ckf = kfupdate(ckf);
        fkf = fkfupdate(ckf, fkf);
    end

    [ckf, ins] = kffeedback(ckf, ins, 1, 'avp');
    [fkf{3}, fins] = kffeedback(fkf{3}, fins, 1, 'avp');

    % save results
    avpc(ki,:)  = [ins.avp', t];
    avpf(ki,:)  = [fins.avp', t];
    xkpkc(ki,:) = [ckf.xk; diag(ckf.Pxk); t]';
	xkpk1(ki,:) = [fkf{1}.xk; diag(fkf{1}.Pxk); t]';
	xkpk2(ki,:) = [fkf{2}.xk; diag(fkf{2}.Pxk); t]';
    xkpkf(ki,:) = [fkf{3}.xk; diag(fkf{3}.Pxk); t]';
    ki = timebar;
end

[avpc,avpf,xkpkc,xkpk1,xkpk2,xkpkf] = no0s(avpc,avpf,xkpkc,xkpk1,xkpk2,xkpkf);
insplot(avpc);
insplot(avpf);
% avpcmpplot(avpc, avpf);
% kfplot(xkpkc);
% kfplot(xkpkf);
kfplot(xkpk1);
kfplot(xkpk2);











