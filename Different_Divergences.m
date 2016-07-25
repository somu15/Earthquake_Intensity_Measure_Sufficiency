clr
% PGA
pdf1 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Approx_Bayes_PGA.txt');
pdf2 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Exact_Bayes_PGA.txt');
cdf1 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Approx_CDF_PGA.txt');
cdf2 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Exact_CDF_PGA.txt');
H_PGA = sqrt(0.5*sum((sqrt(pdf2)-sqrt(pdf1)).^2)*0.05);
L2_PGA = sum(abs(cdf1-cdf2).^2*0.05).^(1/2);
temp1 = log2(pdf1./pdf2);
temp2 = log2(pdf2./pdf1);
temp1(isnan(temp1)) = 0;
temp2(isnan(temp2)) = 0;
temp1(isinf(temp1)) = 0;
temp2(isinf(temp2)) = 0;
KLD_PGA = sum(pdf1.*temp1)*0.05 + sum(pdf2.*temp2)*0.05;
% SaT1.33
pdf1 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Approx_Bayes_SaT1.33.txt');
pdf2 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Exact_Bayes_SaT1.33.txt');
cdf1 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Approx_CDF_SaT1.33.txt');
cdf2 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Exact_CDF_SaT1.33.txt');
H_SaT133 = sqrt(0.5*sum((sqrt(pdf2)-sqrt(pdf1)).^2)*0.05);
L2_SaT133 = sum(abs(cdf1-cdf2).^2*0.05).^(1/2);
temp1 = log2(pdf1./pdf2);
temp2 = log2(pdf2./pdf1);
temp1(isnan(temp1)) = 0 ;
temp2(isnan(temp2)) = 0;
temp1(isinf(temp1)) = 0;
temp2(isinf(temp2)) = 0;
KLD_SaT133 = sum(pdf1.*temp1)*0.05 + sum(pdf2.*temp2)*0.05;
%SaT0.43
pdf1 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Approx_Bayes_SaT0.43.txt');
pdf2 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Exact_Bayes_SaT0.43.txt');
cdf1 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Approx_CDF_SaT0.43.txt');
cdf2 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Exact_CDF_SaT0.43.txt');
H_SaT043 = sqrt(0.5*sum((sqrt(pdf2)-sqrt(pdf1)).^2)*0.05);
L2_SaT043 = sum(abs(cdf1-cdf2).^2*0.05).^(1/2);
temp1 = log2(pdf1./pdf2);
temp2 = log2(pdf2./pdf1);
temp1(isnan(temp1)) = 0 ;
temp2(isnan(temp2)) = 0;
temp1(isinf(temp1)) = 0;
temp2(isinf(temp2)) = 0;
KLD_SaT043 = sum(pdf1.*temp1)*0.05 + sum(pdf2.*temp2)*0.05;
%SaT0.22
pdf1 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Approx_Bayes_SaT0.22.txt');
pdf2 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Exact_Bayes_SaT0.22.txt');
cdf1 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Approx_CDF_SaT0.22.txt');
cdf2 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Exact_CDF_SaT0.22.txt');
H_SaT022 = sqrt(0.5*sum((sqrt(pdf2)-sqrt(pdf1)).^2)*0.05);
L2_SaT022 = sum(abs(cdf1-cdf2).^2*0.05).^(1/2);
temp1 = log2(pdf1./pdf2);
temp2 = log2(pdf2./pdf1);
temp1(isnan(temp1)) = 0 ;
temp2(isnan(temp2)) = 0;
temp1(isinf(temp1)) = 0;
temp2(isinf(temp2)) = 0;
KLD_SaT022 = sum(pdf1.*temp1)*0.05 + sum(pdf2.*temp2)*0.05;
%SaT2
pdf1 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Approx_Bayes_SaT2.txt');
pdf2 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Exact_Bayes_SaT2.txt');
cdf1 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Approx_CDF_SaT2.txt');
cdf2 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Exact_CDF_SaT2.txt');
H_SaT2 = sqrt(0.5*sum((sqrt(pdf2)-sqrt(pdf1)).^2)*0.05);
L2_SaT2 = sum(abs(cdf1-cdf2).^2*0.05).^(1/2);
temp1 = log2(pdf1./pdf2);
temp2 = log2(pdf2./pdf1);
temp1(isnan(temp1)) = 0 ;
temp2(isnan(temp2)) = 0;
temp1(isinf(temp1)) = 0;
temp2(isinf(temp2)) = 0;
KLD_SaT2 = sum(pdf1.*temp1)*0.05 + sum(pdf2.*temp2)*0.05;
% PGV
pdf1 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Approx_Bayes_PGV.txt');
pdf2 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Exact_Bayes_PGV.txt');
cdf1 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Approx_CDF_PGV.txt');
cdf2 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Exact_CDF_PGV.txt');
H_PGV = sqrt(0.5*sum((sqrt(pdf2)-sqrt(pdf1)).^2)*5);
L2_PGV = sum(abs(cdf1-cdf2).^2*5).^(1/2);
temp1 = log2(pdf1./pdf2);
temp2 = log2(pdf2./pdf1);
temp1(isnan(temp1)) = 0 ;
temp2(isnan(temp2)) = 0;
temp1(isinf(temp1)) = 0;
temp2(isinf(temp2)) = 0;
KLD_PGV = sum(pdf1.*temp1)*5 + sum(pdf2.*temp2)*5;



Dr = 0.005:0.005:0.04;
subplot(1,3,1)
plot(Dr,L2_SaT133,Dr,L2_SaT043,Dr,L2_SaT022,Dr,L2_SaT2,Dr,L2_PGA,'linewidth',1.5)
title('Eucledian distance on CDF')
ylabel('Distance root g')
subplot(1,3,2)
plot(Dr,H_SaT133,Dr,H_SaT043,Dr,H_SaT022,Dr,H_SaT2,Dr,H_PGA,Dr,H_PGV,'linewidth',1.5)
title('Hellinger divergence on IM distribution')
xlabel('Roof Drift')
ylabel('Divergence')
legend('SaT1.33','SaT0.43','SaT0.22','SaT2','PGA','PGV')
subplot(1,3,3)
plot(Dr,KLD_SaT133,Dr,KLD_SaT043,Dr,KLD_SaT022,Dr,KLD_SaT2,Dr,KLD_PGA,Dr,KLD_PGV,'linewidth',1.5)
title('KL divergence on IM distribution')
ylabel('KLD')
%SaT1.33
pdf1 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Approx_Bayes_SaT1.33.txt');
pdf2 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Exact_Bayes_SaT1.33.txt');
cdf1 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Approx_CDF_SaT1.33.txt');
cdf2 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Exact_CDF_SaT1.33.txt');
pga_val = 0.01:0.05:10;
figure
subplot(1,2,1)
plot(pga_val,cdf1(:,8),pga_val,cdf2(:,8),'linewidth',1.5)
xlabel('SaT1.33 (g)')
ylabel('P(Dr>x|SaT1.33)')
legend('Assuming conditional independence','Exact')
subplot(1,2,2)
plot(pga_val,pdf1(:,8),pga_val,pdf2(:,8),'linewidth',1.5)
xlabel('SaT1.33 (g)')
ylabel('f(SaT1.33|Dr>x)')
%PGA
pdf1 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Approx_Bayes_PGA.txt');
pdf2 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Exact_Bayes_PGA.txt');
cdf1 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Approx_CDF_PGA.txt');
cdf2 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Bayes and CDFs\0.05g\Exact_CDF_PGA.txt');
figure
subplot(1,2,1)
plot(pga_val,cdf1(:,8),pga_val,cdf2(:,8),'linewidth',1.5)
xlabel('PGA (g)')
ylabel('P(Dr>x|PGA)')
legend('Assuming conditional independence','Exact')
subplot(1,2,2)
plot(pga_val,pdf1(:,8),pga_val,pdf2(:,8),'linewidth',1.5)
xlabel('PGA (g)')
ylabel('f(PGA|Dr>x)')