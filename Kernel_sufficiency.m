clear all
clc
file_path  = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\RD\';
Sa = importdata(strcat(file_path,'x_sa1.5sec.txt'));
eps = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\eps_sa1.txt');

Dr = importdata(strcat(file_path,'y.txt'));
M = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\','Mw.txt'));
R = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\','Rjb.txt'));
Vs = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\','Vs30.txt'));
X = [ones(max(size(Sa)),1) log(Sa)];
beta = regress(log(Dr),X);
res = log(Dr)-(beta(1)+beta(2)*log(Sa));
stdd = std(res);
[yres,xres] = ksdensity(res,'npoints',256);
[ym,xm] = ksdensity(M,'npoints',256);
[xx,yy] = meshgrid(xm,xres);
[px,py] = meshgrid(ym,yres);
pxy = px.*py;
figure()
surf(xx,yy,pxy,'LineStyle','none')