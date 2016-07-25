%Sa1
clr
Dr = 0.001:0.001:0.04;
Sa = 0.01:0.01:10;
dIM = 0.01;
for jj = 1:max(size(Dr))
reg = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Regression details_Medina_Krawinkler\Bilinear_Regression\RD\SaT1.5\beta.txt');
sa1_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Hazard curves\Hazard curves with eps\SaT1.5.txt'));%OpenSHA_SaT1.5_hazard
%AFE_sa1 = exp(interp1(log(sa1_temp(:,1)),log(sa1_temp(:,2)),log(Sa),'spline'));
% AFE_sa1(1) = 0.999999;
% AFE_sa1 = -log(1-AFE_sa1)/50;
dsa1 = abs((Differentiation((dIM),sa1_temp(:,2))));
%std_Sa = 0.2892*ones(1,max(size(Sa)));%sqrt(exp(-2.55-2.126*Sa+0.7047*Sa.^2));
prob_ex_Sa = 1-normcdf(log(Dr(jj)),reg(2)*log(Sa)+reg(1),reg(3));
bayes_Sa = prob_ex_Sa.*dsa1./sum(prob_ex_Sa.*dsa1*dIM);
bayes_Sa = bayes_Sa/trapz(Sa,bayes_Sa);
for ii = 1:max(size(Sa))
    temp = zeros(ii,1);
    temp1 = ones(max(size(Sa))-ii,1);
    psuedo_Sa = vertcat(temp,temp1);
        temp2 = sum(abs(psuedo_Sa-prob_ex_Sa').^2*dIM).^(1/2);
        Div_Sa(ii) = (temp2);
%     temp2 = abs(psuedo_Sa-prob_ex_Sa');
%     Div_Sa(ii) = max(temp2);
end
   A_Sa(jj) = sum(bayes_Sa./Div_Sa*dIM);
end
%A_Sa = 1./A_Sa;
%plot(Sa,prob_ex_Sa)
plot(Dr,A_Sa)
%% Test for epsilon
% Sa1
clr
Dr = 0.003:0.001:0.05;
Sa = 0.01:0.01:10;
dIM = 0.01;
%for jj = 1:max(size(Sa))
    
%end
for jj = 1:max(size(Dr))
sa1_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\SaT0.22_hazard.txt'));%OpenSHA_SaT1.5_hazard
AFE_sa1 = exp(interp1(log(sa1_temp(:,1)),log(sa1_temp(:,2)),log(Sa),'spline'));
% AFE_sa1(1) = 0.999999;
% AFE_sa1 = -log(1-AFE_sa1)/50;
dsa1 = abs((Differentiation((dIM),AFE_sa1)));
std_Sa = 0.2892*ones(1,max(size(Sa)));%sqrt(exp(-2.55-2.126*Sa+0.7047*Sa.^2));
prob_ex_Sa = 1-normcdf((log(Dr(jj))*ones(max(size(Sa)),1)'-(0.4978*log(Sa)-3.5869))./std_Sa,0,1);
sum_temp = 0;
for ii = -4:0.01:4
temp = 1-normcdf(log(Dr(jj)),-3.2672+0.6336*log(Sa)-0.1258*ii,0.2823);
sum_temp = sum_temp+temp*normpdf(ii,0,1)*0.01;
end
prob_ex_eps = sum_temp;
bayes_Sa = prob_ex_Sa.*dsa1./sum(prob_ex_Sa.*dsa1*dIM);
bayes_Sa = bayes_Sa/trapz(Sa,bayes_Sa);
for ii = 1:max(size(Sa))
    temp = zeros(ii,1);
    temp1 = ones(max(size(Sa))-ii,1);
    psuedo_Sa = vertcat(temp,temp1);
        temp2 = sum(abs(prob_ex_eps'-prob_ex_Sa').^2*dIM).^(1/2);
        Div_Sa(ii) = (temp2);
%     temp2 = abs(psuedo_Sa-prob_ex_Sa');
%     Div_Sa(ii) = max(temp2);
end
div(jj) = 1./sum(abs(prob_ex_eps'-prob_ex_Sa').^2*dIM).^(1/2);
   A_Sa(jj) = sum(bayes_Sa./Div_Sa*dIM);
end
%% PGV test cm/sec
clr
Dr = 0.01;
Sa = 1:0.1:1000;
dIM = 0.1;
for jj = 1:max(size(Dr))
sa1_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\PGV_hazard.txt'));%OpenSHA_SaT1.5_hazard
AFE_sa1 = exp(interp1(log(sa1_temp(:,1)),log(sa1_temp(:,2)),log(Sa),'spline'));
% AFE_sa1(1) = 0.999999;
% AFE_sa1 = -log(1-AFE_sa1)/50;
dsa1 = abs((Differentiation((dIM),AFE_sa1)));
std_Sa = 0.2914*ones(1,max(size(Sa)));%sqrt(exp(-2.55-2.126*Sa+0.7047*Sa.^2));
prob_ex_Sa = 1-normcdf((log(Dr(jj))*ones(max(size(Sa)),1)'-(0.8174*log(Sa)-7.4765))./std_Sa,0,1);
bayes_Sa = prob_ex_Sa.*dsa1./sum(prob_ex_Sa.*dsa1*dIM);
prob = normpdf((log(Dr(jj))*ones(max(size(Sa)),1)'-(0.8174*log(Sa)-7.4765))./std_Sa,0,1);
bayes_Sa = bayes_Sa/trapz(Sa,bayes_Sa);
for ii = 1:max(size(Sa))
    temp = zeros(ii,1);
    temp1 = ones(max(size(Sa))-ii,1);
    psuedo_Sa = vertcat(temp,temp1);
    temp2 = sum(abs(psuedo_Sa-prob_ex_Sa').^1*dIM).^(1/1);
    Div_Sa(ii) = (temp2);%temp2 = sum(abs(psuedo_Sa-prob_ex_Sa').^2*dIM).^(1/2);
%     temp2 = abs(psuedo_Sa-prob_ex_Sa');
%     Div_Sa(ii) = mean(temp2);
end
%Div_Sa = -log(sqrt(prob));
   A_Sa(jj) = sum(bayes_Sa./Div_Sa*dIM);
end
%% PGV test m/s
clr
Dr = 0.01;
Sa = 0.01:0.001:10;
dIM = 0.001;
for jj = 1:max(size(Dr))
sa1_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\PGV_hazard.txt'));%OpenSHA_SaT1.5_hazard
AFE_sa1 = exp(interp1(log(sa1_temp(:,1)*0.01),log(sa1_temp(:,2)),log(Sa),'spline'));
% AFE_sa1(1) = 0.999999;
% AFE_sa1 = -log(1-AFE_sa1)/50;
dsa1 = abs((Differentiation((dIM),AFE_sa1)));
std_Sa = 0.2914*ones(1,max(size(Sa)));%sqrt(exp(-2.55-2.126*Sa+0.7047*Sa.^2));
prob_ex_Sa = 1-normcdf((log(Dr(jj))*ones(max(size(Sa)),1)'-(0.8174*log(Sa*100)-7.4765))./std_Sa,0,1);
bayes_Sa = prob_ex_Sa.*dsa1./sum(prob_ex_Sa.*dsa1*dIM);
prob = normpdf((log(Dr(jj))*ones(max(size(Sa)),1)'-(0.8174*log(Sa*100)-7.4765))./std_Sa,0,1);
bayes_Sa = bayes_Sa/trapz(Sa,bayes_Sa);
for ii = 1:max(size(Sa))
    temp = zeros(ii,1);
    temp1 = ones(max(size(Sa))-ii,1);
    psuedo_Sa = vertcat(temp,temp1);
    temp2 = sum(abs(psuedo_Sa-prob_ex_Sa').^1*dIM).^(1/1);
    Div_Sa(ii) = (temp2);
%     temp2 = abs(psuedo_Sa-prob_ex_Sa');
%     Div_Sa(ii) = mean(temp2);
end
   A_Sa(jj) = sum(bayes_Sa./Div_Sa*dIM);
end
%% Plots
drs = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\smaller dataset\results_PGA_SaT1T2_drifts.txt');
pfas = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\smaller dataset\results_PGA_SaT1T2T3_PFAs.txt');
DR = 0.001:0.001:0.05;
PFA = (1:10:772)*2.54e-3;
% Red-PGA
% Blue-SaT1
% Default MATLAB blue - SaT2
% Cyan - SaT3
% Default MATLAB orange - SaGMT1T2
subplot(3,2,1)
plot(DR,drs(1,:),'red',DR,drs(2,:),'blue',DR,drs(3,:),DR,drs(4,:),'cyan',DR,drs(5,:),'linewidth',2)
title('Roof Drift')
grid on
legend('PGA','SaT1','SaT2','SaT3','SaGMT1,T2')
subplot(3,2,2)
plot(DR,drs(16,:),'red',DR,drs(17,:),'blue',DR,drs(18,:),DR,drs(19,:),'cyan',DR,drs(20,:),'linewidth',2)
title('Joint Rotation')
grid on
subplot(3,2,3)
plot(DR,drs(6,:),'red',DR,drs(7,:),'blue',DR,drs(8,:),DR,drs(9,:),'cyan',DR,drs(10,:),'linewidth',2)
title('IDR1')
grid on
subplot(3,2,4)
plot(DR,drs(11,:),'red',DR,drs(12,:),'blue',DR,drs(13,:),DR,drs(14,:),'cyan',DR,drs(15,:),'linewidth',2)
title('IDR4')
grid on
subplot(3,2,5)
plot(PFA,pfas(1,:),'red',PFA,pfas(2,:),'blue',PFA,pfas(3,:),PFA,pfas(4,:),'cyan',PFA,pfas(5,:),'linewidth',2)
title('PFA1')
grid on
subplot(3,2,6)
plot(PFA,pfas(6,:),'red',PFA,pfas(7,:),'blue',PFA,pfas(8,:),PFA,pfas(9,:),'cyan',PFA,pfas(10,:),'linewidth',2)
title('PFA4')
grid on