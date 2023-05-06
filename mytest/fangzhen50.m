clc;
clear all;
close all;

[a]=xlsread('50fangzhen.xlsx');
t=length(a);

figure()

subplot(3,1, 1 ), hold on, plot(a(:,1),'k.-','linewidth',0.5); hold on, plot(a(:,4),'k.:','linewidth',1);ylabel('Rmse-Pe/m'); legend('Classic FKF','Improved AFKF');

subplot(3,1, 2), hold on, plot(a(:,2),'k.-','linewidth',0.5); hold on, plot(a(:,5),'k.:','linewidth',1);ylabel('Rmse-Pn/m');legend('Classic FKF','Improved AFKF');

subplot(3,1, 3), hold on, plot(a(:,3),'k.-','linewidth',0.5); hold on, plot(a(:,6),'k.:','linewidth',1);xlabel('Number of simulations');ylabel('Rmse-Pu/m'); legend('Classic FKF','Improved AFKF');

figure()

subplot(3,1, 1 ), hold on, plot(a(:,7),'k.-','linewidth',0.5); hold on, plot(a(:,10),'k.:','linewidth',1);ylabel('Rmse-Ve（m/s）'); legend('Classic FKF','Improved AFKF');

subplot(3,1, 2), hold on, plot(a(:,8),'k.-','linewidth',0.5); hold on, plot(a(:,11),'k.:','linewidth',1);ylabel('Rmse-Vn（m/s）');legend('Classic FKF','Improved AFKF');

subplot(3,1, 3), hold on, plot(a(:,9),'k.-','linewidth',0.5); hold on, plot(a(:,12),'k.:','linewidth',1);xlabel('Number of simulations');ylabel('Rmse-Vu（m/s）'); legend('Classic FKF','Improved AFKF');

p1=mean(a(:,1)); p4=mean(a(:,4));
p2=mean(a(:,2)); p5=mean(a(:,5));
p3=mean(a(:,3)); p6=mean(a(:,6));

v1=mean(a(:,7)); v4=mean(a(:,10));
v2=mean(a(:,8)); v5=mean(a(:,11));
v3=mean(a(:,9)); v6=mean(a(:,12));