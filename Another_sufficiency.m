clr
format long
% Data imports
str1 = 'PGA';
str2 = 'pga';
beta = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Regression details\Bilinear_Regression\RD\',str1,'\beta.txt'));
beta_M = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Regression details\Bilinear_Regression\RD\',str1,'\beta_M.txt'));
beta_R = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Regression details\Bilinear_Regression\RD\',str1,'\beta_R.txt'));
beta_eps = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Regression details\Bilinear_Regression\RD\',str1,'\beta_eps.txt'));
eps = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\eps_',str2,'.txt'));
file_path  = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\RD\';

M = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\','Mw.txt'));
R = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\','Rjb.txt'));
Sa = importdata(strcat(file_path,'x_',str2,'.txt'));
Dr = importdata(strcat(file_path,'y.txt'));
D_M = abs(sum(log(normpdf(log(Dr),beta(1)+beta(2)*log(Sa),beta(3))))-sum(log(normpdf(log(Dr),beta_M(1)+beta_M(2)*log(Sa)+beta_M(3)*M,beta_M(4)))));
D_R = abs(sum(log(normpdf(log(Dr),beta(1)+beta(2)*log(Sa),beta(3))))-sum(log(normpdf(log(Dr),beta_R(1)+beta_R(2)*log(Sa)+beta_R(3)*R,beta_R(4)))));
D_eps = abs(sum(log(normpdf(log(Dr),beta(1)+beta(2)*log(Sa),beta(3))))-sum(log(normpdf(log(Dr),beta_eps(1)+beta_eps(2)*log(Sa)+beta_eps(3)*eps,beta_eps(4)))));
D = (D_M+D_R+D_eps)/max(size(Dr))

