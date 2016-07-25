clear all
clc
file_path  = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\FEMA records unscaled\PFA1\';
Sa = importdata(strcat(file_path,'x_pga.txt'));
Dr = importdata(strcat(file_path,'y.txt'));
x = log(Sa);
y = log(Dr);
beta = regress(y,[ones(44,1) x]);
beta
res = y-(beta(1)+beta(2)*x);
stdd = std(res);
slope1 = beta(2)-tinv(1-0.025,42)*sqrt(stdd^2/sum(res.^2));
inter1 = beta(1)-tinv(1-0.025,42)*sqrt(stdd^2/44);
slope2 = beta(2)+tinv(1-0.025,42)*sqrt(stdd^2/sum(res.^2));
inter2 = beta(1)+tinv(1-0.025,42)*sqrt(stdd^2/44);
[slope1 slope2]
[inter1 inter2]
% for jj = 1:44
% for ii = 1:1000
% index = randsample(44,jj);
% x1 = Sa(index);
% y1 = Dr(index);
% X1 = [ones(max(size(x1)),1) log(x1)];
% beta1(ii,:) = regress(log(y1),X1);
% res1 = log(y1)-(beta1(ii,1)+beta1(ii,2)*log(x1));
% std1(ii) = std(res1);
% end
% %mean(beta1)
% std_beta(:,jj) = std(beta1);
% %mean(std1)
% std_std(jj) = std(std1);
% end
% subplot(1,2,1)
% plot(5:44,std_sa1(2,5:44),'red',5:44,std_sa2(2,5:44),'b',5:44,std_pga(2,5:44),5:44,std_sa3(2,5:44),5:44,std_sa2s(2,5:44),'linewidth',1.5)
% legend('SaT1.33','SaT0.43','PGA','SaT0.22','SaT2')
% subplot(1,2,2)
% plot(5:44,std_sa1(1,5:44),'red',5:44,std_sa2(1,5:44),'b',5:44,std_pga(1,5:44),5:44,std_sa3(1,5:44),5:44,std_sa2s(1,5:44),'linewidth',1.5)


