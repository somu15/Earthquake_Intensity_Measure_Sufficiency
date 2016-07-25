clear all
clc
file_path  = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\RD\';
name = 'sa1';
Sa = importdata(strcat(file_path,'x_',name,'.txt'));
eps = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\eps_',name,'.txt'));

Dr = importdata(strcat(file_path,'y.txt'));
M = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\','Mw.txt'));
R = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\','Rjb.txt'));
Vs = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\','Vs30.txt'));

X = [ones(max(size(Sa)),1) log(Sa)];
beta = regress(log(Dr),X);
res = log(Dr)-(beta(1)+beta(2)*log(Sa));
% R
[bandwidth,density2d_R,X_res1,X_R] = kde2d([res,R]);
x_R = X_R(:,1);
x_res1 = X_res1(1,:)';
dR = abs(x_R(2)-x_R(1));
dres1 = abs(x_res1(2)-x_res1(1));
density_R = ksdensity(R,x_R);
density_res1 = ksdensity(res,x_res1);
[density1,density2] = meshgrid(density_res1,density_R);
density = density1.*density2;
suff_R = density2d_R.*log2(density2d_R./density)*dR*dres1;
suff_R(isnan(suff_R)) = 0;
    suff_R(isinf(suff_R)) = 0;
    suff_R = sum(sum(suff_R));
    
    % M
[bandwidth,density2d_R,X_res1,X_R] = kde2d([res,M]);
x_R = X_R(:,1);
x_res1 = X_res1(1,:)';
dR = abs(x_R(2)-x_R(1));
dres1 = abs(x_res1(2)-x_res1(1));
density_R = ksdensity(M,x_R);
density_res1 = ksdensity(res,x_res1);
[density1,density2] = meshgrid(density_res1,density_R);
density = density1.*density2;
suff_M = density12.*log2(density12./density)*dR*dres1;
suff_M(isnan(suff_M)) = 0;
    suff_M(isinf(suff_M)) = 0;
    suff_M = sum(sum(suff_M));

   % eps
[bandwidth,density2d_R,X_res1,X_R] = kde2d([res,eps]);
x_R = X_R(:,1);
x_res1 = X_res1(1,:)';
dR = abs(x_R(2)-x_R(1));
dres1 = abs(x_res1(2)-x_res1(1));
density_R = ksdensity(eps,x_R);
density_res1 = ksdensity(res,x_res1);
[density1,density2] = meshgrid(density2,density1);
density11 = density1.*density2;
suff_eps = density2d_R.*log2(density2d_R./density)*dR*dres1;
suff_eps(isnan(suff_eps)) = 0;
    suff_eps(isinf(suff_eps)) = 0;
    suff_eps = sum(sum(suff_eps));
suff = suff_M+suff_R+suff_eps;
suff = [suff_M suff_R suff_eps suff];



 %surf(X_res,Y_res,density,'LineStyle','none')
res1 = linspace(-0.8788,0.8297,256);
M1 = linspace(5.94,8.18,256);
S = cov(M,res);
S(1,2) = 0;
S(2,1) = 0;
d = 2;
n = max(size(M));
h  = (4/(d+2))^(1/(d+4))*n^(-1/(d+4));
for ii = 1:256
    for jj = 1:256
    y = [M1(ii);res1(jj)];
    sum1 = 0;
    for kk = 1:max(size(M))
    yi = [M(kk);res(kk)];
        u = (y-yi)'*inv(S)*(y-yi)/h^2;
        K = 1/((2*pi)^(d/2)*h^d*det(S)^0.5)*exp(-u/2);
        sum1 = sum1+K;
    end
    density12(ii,jj) = sum1/n;
    end
end
% [xx1,yy1] = meshgrid(res1,M1);
% figure()
% surf(xx1,yy1,density12,'LineStyle','no')
% d = 1;
% n = max(size(M));
% h  = (4/(d+2))^(1/(d+4))*n^(-1/(d+4));
S = cov(M);
S2 = cov(res);
for ii = 1:256
    y = M1(ii);
    y2 = res1(ii);
    sum1 = 0;
    sum2 = 0;
for jj = 1:max(size(M))
yi = M(jj);
yi2 = res(jj);
u = (y-yi)'*inv(S)*(y-yi)/h^2;
u2 = (y2-yi2)'*inv(S2)*(y2-yi2)/h^2;
K = 1/((2*pi)^(d/2)*h^d*det(S)^0.5)*exp(-u/2);
        sum1 = sum1+K;
        K2 = 1/((2*pi)^(d/2)*h^d*det(S2)^0.5)*exp(-u2/2);
        sum2 = sum2+K2;
end
density1(ii) = sum1/n;
density2(ii) = sum2/n;
end
suff_eps = density12.*log2(density12./density11)*dR*dres1;
suff_eps(isnan(suff_eps)) = 0;
    suff_eps(isinf(suff_eps)) = 0;
    suff_eps = sum(sum(suff_eps))