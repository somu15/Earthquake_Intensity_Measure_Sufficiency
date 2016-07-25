clr
haz_sa2 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Hazard curves\Hazard curves with eps\SaT0.43.txt');
haz_pga = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Hazard curves\Hazard curves with eps\PGA.txt');
haz_sa1 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Hazard curves\Hazard curves with eps\SaT1.33.txt');
haz_sa3 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Hazard curves\Hazard curves with eps\SaT0.22.txt');
haz_sa2s = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Hazard curves\Hazard curves with eps\SaT2.txt');
IM = 0.01:0.01:10;
dsa2 = abs((Differentiation(0.01,haz_sa2(:,2))));
dpga = abs((Differentiation(0.01,haz_pga(:,2))));
dsa1 = abs((Differentiation(0.01,haz_sa1(:,2))));
dsa3 = abs((Differentiation(0.01,haz_sa3(:,2))));
dsa2s = abs((Differentiation(0.01,haz_sa2s(:,2))));
Dr = 0.001:0.001:0.04;
% FEMA
for ii = 1:max(size(Dr))
prob_ex_sa1 = 1-normcdf(log(Dr(ii)),-3.57069+0.72330*log(IM),0.16953);
afe_sa1(ii) = sum(prob_ex_sa1.*dsa1)*0.01;

prob_ex_sa2 = 1-normcdf(log(Dr(ii)),-4.40940+0.27561*log(IM),0.38218);
afe_sa2(ii) = sum(prob_ex_sa2.*dsa2)*0.01;

prob_ex_sa3 = 1-normcdf(log(Dr(ii)),-4.43620+0.24382*log(IM),0.39605);
afe_sa3(ii) = sum(prob_ex_sa3.*dsa3)*0.01;

prob_ex_sa2s = 1-normcdf(log(Dr(ii)),-3.58685+0.49777*log(IM),0.28918);
afe_sa2s(ii) = sum(prob_ex_sa2s.*dsa2s)*0.01;

prob_ex_pga = 1-normcdf(log(Dr(ii)),-4.0770+0.39839*log(IM),0.38007);
afe_pga(ii) = sum(prob_ex_pga.*dpga)*0.01;
end
ES_FEMA = 1./[1;7.4074;7.8691;1.6534;3.8284];
ES_FEMA = ES_FEMA/sum(ES_FEMA);
afe_FEMA = afe_sa1*ES_FEMA(1)+afe_sa2*ES_FEMA(2)+afe_sa3*ES_FEMA(3)+afe_sa2s*ES_FEMA(4)+afe_pga*ES_FEMA(5);
% M-K
for ii = 1:max(size(Dr))
prob_ex_sa1 = 1-normcdf(log(Dr(ii)),-3.31410+0.89976*log(IM),0.13914);
afe_sa1_1(ii) = sum(prob_ex_sa1.*dsa1)*0.01;

prob_ex_sa2 = 1-normcdf(log(Dr(ii)),-4.3156+0.88909*log(IM),0.48237);
afe_sa2_1(ii) = sum(prob_ex_sa2.*dsa2)*0.01;

prob_ex_sa3 = 1-normcdf(log(Dr(ii)),-4.5359+0.7540*log(IM),0.5679);
afe_sa3_1(ii) = sum(prob_ex_sa3.*dsa3)*0.01;

prob_ex_sa2s = 1-normcdf(log(Dr(ii)),-3.1891+0.76179*log(IM),0.41973);
afe_sa2s_1(ii) = sum(prob_ex_sa2s.*dsa2s)*0.01;

prob_ex_pga = 1-normcdf(log(Dr(ii)),-3.64084+0.89563*log(IM),0.50970);
afe_pga_1(ii) = sum(prob_ex_pga.*dpga)*0.01;
end
ES_MK = 1./[2.003;1.8619;4.4477;2.1768;2.0847];
ES_MK = ES_MK/sum(ES_MK);
afe_MK = afe_sa1_1*ES_MK(1)+afe_sa2_1*ES_MK(2)+afe_sa3_1*ES_MK(3)+afe_sa2s_1*ES_MK(4)+afe_pga_1*ES_MK(5);
loglog(Dr,afe_sa1,Dr,afe_MK,'linewidth',1.5)