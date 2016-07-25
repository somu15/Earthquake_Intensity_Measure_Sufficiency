%%
clear all
clc
file_path  = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\RD\';
Sa = importdata(strcat(file_path,'x_sa1.5sec.txt'));
eps = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\eps_sa1.txt');

Dr = importdata(strcat(file_path,'y.txt'));
M = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\','Mw.txt'));
R = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\','Rjb.txt'));
Vs = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\','Vs30.txt'));
% R(35) = 0.1;
% R(36) = 0.1;
X = [ones(max(size(Sa)),1) log(Sa)];
X_M = [ones(max(size(Sa)),1) log(Sa) M];
X_R = [ones(max(size(Sa)),1) log(Sa) R];
X_eps = [ones(max(size(Sa)),1) log(Sa) eps];
X_Vs = [ones(max(size(Sa)),1) log(Sa) Vs];
beta = regress(log(Dr),X);
beta_M = regress(log(Dr),X_M);
beta_R = regress(log(Dr),X_R);
beta_eps = regress(log(Dr),X_eps);
beta_Vs = regress(log(Dr),X_Vs);
res = log(Dr)-(beta(1)+beta(2)*log(Sa));
res_M = log(Dr)-(beta_M(1)+beta_M(2)*log(Sa)+beta_M(3)*M);
res_R = log(Dr)-(beta_R(1)+beta_R(2)*log(Sa)+beta_R(3)*R);
res_eps = log(Dr)-(beta_eps(1)+beta_eps(2)*log(Sa)+beta_eps(3)*eps);
res_Vs = log(Dr)-(beta_Vs(1)+beta_Vs(2)*log(Sa)+beta_Vs(3)*Vs);
stdd = std(res);
std_M = std(res_M);
std_R = std(res_R);
std_eps = std(res_eps);
std_Vs = std(res_Vs);
beta(3) = stdd;
beta_M(4) = std_M;
beta_R(4) = std_R;
beta_eps(4) = std_eps;
beta_Vs(4) = std_Vs;

str_M = regstats(res,M,'linear',{'tstat'});
str_R = regstats(res,R,'linear',{'tstat'});
str_eps = regstats(res,eps,'linear',{'tstat'});
str_Vs = regstats(res,Vs,'linear',{'tstat'});

colldiag(X)
colldiag(X_M)
colldiag(X_R)
colldiag(X_eps)
colldiag(X_Vs)
cov_matrix = stdd^2*inv(X'*X)
cov_matrix_M = std_M^2*inv(X_M'*X_M)
cov_matrix_R = std_R^2*inv(X_R'*X_R)
cov_matrix_eps = std_eps^2*inv(X_eps'*X_eps)
cov_matrix_Vs = std_Vs^2*inv(X_Vs'*X_Vs)
KS = kstest(res/stdd)
AD = adtest(res/stdd)
KS_M = kstest(res_M/std_M)
AD_M = adtest(res_M/std_M)
KS_R = kstest(res_R/std_R)
AD_R = adtest(res_R/std_R)
KS_eps = kstest(res_eps/std_eps)
AD_eps = adtest(res_eps/std_eps)
KS_Vs = kstest(res_Vs/std_Vs)
AD_Vs = adtest(res_Vs/std_Vs)
str_M.tstat.pval
str_R.tstat.pval
str_eps.tstat.pval
str_Vs.tstat.pval