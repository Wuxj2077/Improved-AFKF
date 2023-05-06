%这个程序为fkf.beta为1/2  不是自适应调节
%
%

clc;
close all;
clear;
glvs
trj = trjfile('trj10ms.mat');
[nn, ts, nts] = nnts(2, trj.ts);

% GPS simulator
lever = [1; 2; 3]*0; rk1 = vperrset(0.5, 3);
davp0 = avperrset([0.5;-0.5;20], 0.5, [3;3;3]);
gps = gpssimu(trj.avp, davp0(4:6), davp0(7:9), 1, lever, 0.0);  % gpsplot(gps)
% gpsplot(gps);

% vio simulator
mu = 0*[1;2;5]*glv.min;  rk2 = [[0.3;0.3;0.3]*glv.sec];
[qis, utc0] = viosimu(trj.avp, rk2, mu, [2021;11;22;12*3600; -0.1;37]);

% SINS init
imuerr = imuerrset(0.5, [1000;1000;1000], 0.5, 100);
imu = imuadderr(trj.imu, imuerr);  % imuplot(imu)
ins = insinit(trj.avp0(1:9), ts, davp0);  
len = length(imu);
Cie0 = cnsCie(utc0(1:3), utc0(4),  utc0(5), utc0(6));
% imuplot(imu)

% centralized KF setting
psinstypedef(156);
ckf = [];
ckf.Phikk_1 = eye(15); 
ckf.Qt = diag([imuerr.web; imuerr.wdb; zeros(9,1)])^2;  % 15-state
ckf.Rk = diag([rk1; rk2])^2;
ckf.Pxk = diag([davp0; imuerr.eb; imuerr.db]*1.0)^2;
ckf.Hk = [ zeros(6,3), eye(6), zeros(6,6); ...    % SINS/GPS Hk
          eye(3), zeros(3,12) ];                  % SINS/VIO Hk
ckf.Hk(7,7)=1;  
ckf = kfinit0(ckf, nts);

% federated KF setting
fins = ins;
fkf = fkfinit(ckf, {1:15,1:15}, {1:6,7:9}, [1/2; 1/2]);

% alloc memory 分配内存
[avpc, xkpkc] = prealloc(fix(len/nn), 10, 2*ckf.n+1);
[avp1, xkpk1] = prealloc(fix(len/nn), 10, 2*fkf{1}.n+1);
[avp2, xkpk2] = prealloc(fix(len/nn), 10, 2*fkf{2}.n+1);
[avpf, xkpkf] = prealloc(fix(len/nn), 10, 2*fkf{3}.n+1);
ki = timebar(nn, len, '15-state SINS/GPS/VIO federated KF simulation.');

imugpssyn(imu(:,7), gps(:,end)); 
tic
for k=1:nn:len-nn+1
    k1 = k+nn-1;

    wvm = imu(k:k1,1:6);  %wvm为imu的数据 加计和陀螺    
    t = imu(k1,end);

    ins = insupdate(ins, wvm);
    fins = insupdate(fins, wvm);
    ckf.Phikk_1 = kffk(fins);

    [kgps, dt] = imugpssyn(k, k1, 'F');

    if kgps>0 && mod(t,1)<1.5*ts && norm(ins.wnb)<1*glv.dps  % SINS/GPS/VIO

        Cns = cnsCns(qis(k1,1:3)', fins.pos, Cie0, t);
        ckf.Hk(8:9,8) = -[fins.eth.cl; fins.eth.sl];
        zk = [fins.vn-gps(kgps,1:3)'; fins.pos-gps(kgps,4:6)'; qq2phi(fins.qnb,m2qua(Cns))];
        fkf = fkfupdate(ckf, fkf, zk);  % federated KF
    else
        
        fkf = fkfupdate(ckf, fkf);
    end
    
    [fkf{3}, fins] = kffeedback(fkf{3}, fins, 1, 'avp');

    % save results
    avpc(ki,:)  = [ins.avp', t];
    avpf(ki,:)  = [fins.avp', t];

	xkpk1(ki,:) = [fkf{1}.xk; diag(fkf{1}.Pxk); t]';
	xkpk2(ki,:) = [fkf{2}.xk; diag(fkf{2}.Pxk); t]';
    xkpkf(ki,:) = [fkf{3}.xk; diag(fkf{3}.Pxk); t]';
    ki = timebar;
    
end
toc
%% 误差分析及画图

[avpc,avpf,xkpk1,xkpk2,xkpkf] = no0s(avpc,avpf,xkpk1,xkpk2,xkpkf);

insplot(avpf);
% avpcmpplot(avpc, avpf);   %这个是惯导积分得到的轨迹 与fkf得到的轨迹比较 没有意义
% kfplot(xkpkc);
% kfplot(xkpkf);
% kfplot(xkpk1);
% kfplot(xkpk2);

fkfpose3=pos2dxyz(avpf(:,7:9));
t=avpf(:,end);
truepose3=pos2dxyz(trj.avp(1:2:end,7:9));
es=fkfpose3(:,1:3)-truepose3(:,1:3);
esv=avpf(:,4:6)-trj.avp(1:2:end,4:6);
figure()
subplot(2,1,[1 2]), hold on, plot(t,es); legend('dE','dU','dN');
figure()
subplot(2,1,[1 2]), hold on, plot(t,esv); legend('dEv','dUv','dNv');



