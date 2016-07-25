clr
haz_sa1 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Hazard curves\Hazard curves with eps\SaT1.33.txt');
haz_sa2 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Hazard curves\Hazard curves with eps\SaT0.43.txt');
IM = 0.01:0.01:10;
dsa1 = abs((Differentiation(0.01,haz_sa1(:,2))));
dsa2 = abs((Differentiation(0.01,haz_sa2(:,2))));
Dr = 0.001:0.001:0.04;
for ii = 1:max(size(Dr))
prob_ex_sa1 = 1-normcdf(log(Dr(ii)),-3.57069+0.72330*log(IM),0.16953);
prob_ex_sa2 = 1-normcdf(log(Dr(ii)),-4.31566+0.88909*log(IM),0.48237);
afe_sa1(ii) = sum(prob_ex_sa1.*dsa1)*0.01;
afe_sa2(ii) = sum(prob_ex_sa2.*dsa2)*0.01;
end
loglog(Dr,afe_sa1,Dr,afe_sa2,'linewidth',1.5)
diff = sum(abs(afe_sa1-afe_sa2))