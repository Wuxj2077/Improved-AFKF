clc;
clear all;
close all;

[a]=xlsread('0305fangzhenjieguo.xlsx');
t=length(a);
times=[];

for i=1:1:t
    if i<2
    times(i,1)=0;
    else
    times(i,1)=(i-1)*0.02;
    end
end

d=[];
d=[a times];

figure()

subplot(3,1, 1 ), hold on, plot(d(:,13),d(:,1),'k:','linewidth',1.5); hold on, plot(d(:,13),d(:,2),'k-','linewidth',1);xlabel('Times（s）');ylabel('dPe（m）'); legend('Classic FKF','Improved AFKF');

subplot(3,1, 2), hold on, plot(d(:,13),d(:,3),'k:','linewidth',1.5); hold on, plot(d(:,13),d(:,4),'k-','linewidth',1); xlabel('Times（s）');ylabel('dPn（m）');legend('Classic FKF','Improved AFKF');

subplot(3,1, 3), hold on, plot(d(:,13),d(:,5),'k:','linewidth',1.5); hold on, plot(d(:,13),d(:,6),'k-','linewidth',1);xlabel('Times（s）');ylabel('dPu（m）'); legend('Classic FKF','Improved AFKF');

figure()
subplot(3,1, 1), hold on, plot(d(:,13),d(:,7),'k:','linewidth',1.5); hold on, plot(d(:,13),d(:,8),'k-','linewidth',1);xlabel('Times（s）');ylabel('dVe（m/s）');legend('Classic FKF','Improved AFKF');

subplot(3,1, 2), hold on, plot(d(:,13),d(:,9),'k:','linewidth',1.5); hold on, plot(d(:,13),d(:,10),'k-','linewidth',1);xlabel('Times（s）');ylabel('dVn（m/s）');legend('Classic FKF','Improved AFKF');

subplot(3,1, 3), hold on, plot(d(:,13),d(:,11),'k:','linewidth',1.5); hold on, plot(d(:,13),d(:,12),'k-','linewidth',1);xlabel('Times（s）');ylabel('dVu（m/s）');legend('Classic FKF','Improved AFKF');

d=abs(d);
%% 求位置误差的均值
%250-400
Pe_mean1=mean(mean(d(12500:20000,1))); Pe_adpmean1=mean(mean(d(12500:20000,2))); Pn_mean1=mean(mean(d(12500:20000,3))); Pn_adpmean1=mean(mean(d(12500:20000,4)));Pu_mean1=mean(mean(d(12500:20000,5))); Pu_adpmean1=mean(mean(d(12500:20000,6)));

%600-750s
Pe_mean2=mean(mean(d(30000:32500,1))); Pe_adpmean2=mean(mean(d(30000:32500,2))); Pn_mean2=mean(mean(d(30000:32500,3))); Pn_adpmean2=mean(mean(d(30000:32500,4))); Pu_mean2=mean(mean(d(30000:32500,5))); Pu_adpmean2=mean(mean(d(30000:32500,6)));

%无干扰时间段
% Pe_mean=mean(mean(d(25000:40000,1))); Pe_adpmean=mean(mean(d(25000:40000,2)));
wd1=d(1:12499,:); wd2=d(20001:29999,:); wd3=d(35001:end,:);  wd=[wd1;wd2;wd3];
Pe_mean3=mean(mean(wd(:,1)));     Pe_adpmean3=mean(mean(wd(:,2))); Pn_mean3=mean(mean(wd(:,3)));     Pn_adpmean3=mean(mean(wd(:,4))); Pu_mean3=mean(mean(wd(:,5)));     Pu_adpmean3=mean(mean(wd(:,6)));



%全部时间段
Pe_mean4=mean(mean(d(:,1)));     Pe_adpmean4=mean(mean(d(:,2))); Pn_mean4=mean(mean(d(:,3)));     Pn_adpmean4=mean(mean(d(:,4))); Pu_mean4=mean(mean(d(:,5)));     Pu_adpmean4=mean(mean(d(:,6)));

%% 求速度误差的均值
%250-400
Ve_mean1=mean(mean(d(12500:20000,7))); Ve_adpmean1=mean(mean(d(12500:20000,8))); Vn_mean1=mean(mean(d(12500:20000,9))); Vn_adpmean1=mean(mean(d(12500:20000,10)));Vu_mean1=mean(mean(d(12500:20000,11))); Vu_adpmean1=mean(mean(d(12500:20000,12)));

%600-700s
Ve_mean2=mean(mean(d(30000:32500,7))); Ve_adpmean2=mean(mean(d(30000:32500,8))); Vn_mean2=mean(mean(d(30000:32500,9))); Vn_adpmean2=mean(mean(d(30000:32500,10))); Vu_mean2=mean(mean(d(30000:32500,11))); Vu_adpmean2=mean(mean(d(30000:32500,12)));

%无干扰时间段
Ve_mean3=mean(mean(wd(:,7)));           Ve_adpmean3=mean(mean(wd(:,8)));           Vn_mean3=mean(mean(wd(:,9)));            Vn_adpmean3=mean(mean(wd(:,10)));        Vu_mean3=mean(mean(wd(:,11)));          Vu_adpmean3=mean(mean(wd(:,12)));

%全部时间段
Ve_mean4=mean(mean(d(:,7)));           Ve_adpmean4=mean(mean(d(:,8)));           Vn_mean4=mean(mean(d(:,9)));            Vn_adpmean4=mean(mean(d(:,10)));        Vu_mean4=mean(mean(d(:,11)));          Vu_adpmean4=mean(mean(d(:,12)));
