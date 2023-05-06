function [fkf, testHk] = errofkfupdate(ckf, fkf, yk)
% federated Kalman filter updating.
%
% Prototype: fkf = fkfupdate(ckf, fkf, zk)
% Inputs: ckf - centralized KF structure, just for FKF Phikk_1/Hk setting;  ckf - 集中式 KF 结构，仅用于 FKF Phikk_1/Hk 设置； 
%         fkf - federated Kalman filter structure cell                      fkf- 联邦卡尔曼滤波器结构单元
%         yk - measurement array same as input to ckf                       yk - 测量数组与 ckf 的输入相同
%
% Outputs: fkf - federated Kalman filter structure cell                     fkf - 联合卡尔曼滤波器结构单元
%          testHk - must be 0, indicating CKF measurement consistency with  FKF表明 CKF 测量与 FKF 的一致性
%
% See also  fkfinit, kfupdate, kfinit.

% Copyright(c) 2009-2022, by Gongmin Yan, All rights reserved.
% Northwestern Polytechnical University, Xi An, P.R.China
% 29/01/2022
    len = length(fkf);
    zk = cell(len,1);

    if nargin<3, 
        for k=1:len, 
            zk{k}=[];
        end;
    else, 
        for k=1:len, 
            zk{k}=yk(fkf{k}.subm);       %子滤波器对应的测量值   子滤波器1：[1,2,3,4,5,6] 2:[7,8,9]
        end;   
    end

    iPg = 0;   iPkxk = 0; fHk = ckf.Hk*0;

     for k=1:len   %len=3;
        if fkf{k}.fkfbeta<1e-10, continue; end  % fkfbeta<=0 表示无信息融合 fkf{k}  
                                                
        sm = fkf{k}.subm;  %subn={1:15,1:15,1:15};subm={1:6,7:9,[]}
        sn = fkf{k}.subn; 
        ncom = fkf{k}.ncom;   % 在每次循环中ncom=15
        fkf{k}.Phikk_1 = ckf.Phikk_1(sn,sn);  %x=Fx+w  中的 转移矩阵F
        fkf{k}.Qk = ckf.Qk(sn,sn);   % NOTE: ckf.Gammak=1
        fkf{k}.Hk = ckf.Hk(sm,sn);       
        fkf{k}.Rk = ckf.Rk(sm,sm);
        fHk(sm,sn) = ckf.Hk(sm,sn);

       %fkf{k}.fkfreset==1 一直存在   %根据主滤波器来重置子滤波器X, P, Q  
        if fkf{k}.fkfreset==1 
            if k<len  %子滤波器重置X, P, Q       %221003 感觉这个if end没啥用   最终目的是让fkf{k}.Pxk=fkf{end}.Pxk                        
                fkf{k}.xk(1:ncom) = fkf{end}.xk;   %注意：公共系统状态必须是前 1：nfkf elements
                sP = chol(fkf{k}.Pxk(1:ncom,1:ncom));sPfkf = chol(fkf{end}.Pxk); A = sP^-1 * sPfkf;    
                fkf{k}.Pxk = [fkf{end}.Pxk, A'*fkf{k}.Pxk(1:ncom,ncom+1:end); fkf{k}.Pxk(ncom+1:end,1:ncom)*A, fkf{k}.Pxk(ncom+1:end,ncom+1:end)];
            end
            fkf{k}.Pxk = PQbeta(fkf{k}.Pxk, fkf{k}.fkfbeta, ncom, fkf{k}.pmethod);  % enlarge P,Q   %fkf{k}.pmethod一直等于3
            fkf{k}.Qk  = PQbeta(fkf{k}.Qk,  fkf{k}.fkfbeta, ncom, fkf{k}.pmethod);
        end      
      
        %% 时间更新 量测更新  
        % 如果zk{k}为空的，只进行时间更新，不进行量测更新
        if isempty(zk{k}), 
            fkf{k} = kfupdate(fkf{k});  % KF update
        else  
            fkf{k} = kfupdate(fkf{k}, zk{k}); 
        end
        
%    if k<len
        iPk = invbc(fkf{k}.Pxk(1:ncom,1:ncom));
        iPg = iPg + iPk;                        %子滤波器与主滤波器的P逆的和  % FKF info fusion 初值iPg=0
        iPkxk = iPkxk + iPk*fkf{k}.xk(1:ncom);  %子滤波器与主滤波器的P逆*Xi的和 初值iPkxk=0
%    end

    end    

 %% 自适应分配因子
% 
%         Pz = 0;
%         for i=1:2, Pz = Pz+1/((trace((fkf{1,i}.Pxk)*((fkf{1,i}.Pxk))') )^(1/2));end
%         for i=1:2, fkf{1,i}.fkfbeta=(1/((trace((fkf{1,i}.Pxk)*((fkf{1,i}.Pxk))') )^(1/2)))/Pz;end

 %%  主滤波器的协方差矩阵P 状态向量X  下一时刻的

    fkf{end}.Pxk = invbc(iPg); 
    fkf{end}.xk = fkf{end}.Pxk * iPkxk;  
    
    testHk = norm(fHk-ckf.Hk,inf);   
    %% PQbete函数 就是为了求子滤波器的P Q

function Pxk = PQbeta(Pxk, beta, ncom, pmethod)   % enlarge P,Q, Eq.(6.11.26b)
	if pmethod==1,
        Pxk = Pxk * (1/beta);
	elseif pmethod==2,
        Pxk(1:ncom,1:ncom) = Pxk(1:ncom,1:ncom) * (1/beta);
	elseif pmethod==3,
        sb = sqrt(1/beta);
        Pxk(1:ncom,:) = Pxk(1:ncom,:) * sb; 
        Pxk(:,1:ncom) = Pxk(:,1:ncom) * sb;
    end

