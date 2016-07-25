%% Sufficiency
clr
Dr = 0.005:0.005:0.04;
sno = 1;%%%%%%%%
Tp = 2;%%%%%%%
IM_min = 1;
IM_max = 1000;
dIM = 5;
pga_val = IM_min:dIM:IM_max;
reg_full = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Regression details\RD\PGV_MRe.txt');
reg = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Regression details\RD\PGV.txt');
sa1_temp = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Hazard curves\Hazard curves with eps\PGV.txt');
dsa1 = abs((Differentiation((5),sa1_temp(:,2))));
%dsa1 = dsa1(1:5:1000);%%%%%%%%%
for hh = 1:max(size(Dr))
for l = 1:max(size(pga_val))
pga = pga_val(l);
M = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\M.txt');
R = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\R.txt');
eps_max = 4;
intervals = 20;
eps_min = -eps_max;
deps = (eps_max-eps_min)/(intervals-1);
eps = eps_min:deps:eps_max;
dM = M(2)-M(1);
dR = R(2)-R(1);
[fMR] = PSHA_BA08_IMSelection_Paper(pga,sno,Tp,0,760);%importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\SaT1.33_Deagg\f(M,R)_',num2str(pga),'g.txt'));
prob_ex_suff = zeros(20,20);
 for ii = 1:20
     for jj = 1:20    
%         eps_sum = 0;
temp_sum = 0;
        for kk = 1:20
temp = 1-normcdf(log(Dr(hh)),reg_full(1)+reg_full(2)*log(pga_val(l))+reg_full(3)*M(ii)+reg_full(4)*(R(jj))+reg_full(5)*eps(kk),reg_full(6));
       full_sum(ii,jj,kk) = temp;
        end
%temp_sum = temp_sum + temp*fMR(ii,jj);
     end
 end
full_sum = sum(sum(sum(full_sum.*fMR)));
cdf_req(l) = full_sum;
end
prob_ex_approx = 1-normcdf(log(Dr(hh)),reg(2)*log(pga_val)+reg(1),reg(3));
for ii = 1:max(size(pga_val))
    temp = zeros(ii,1);
    temp1 = ones(max(size(pga_val))-ii,1);
    psuedo_Sa = vertcat(temp,temp1);
        temp2 = sum(abs(psuedo_Sa-prob_ex_approx').^2*dIM).^(1/2);
        Div_NC(ii) = (temp2);
        
    temp = zeros(ii,1);
    temp1 = ones(max(size(pga_val))-ii,1);
    psuedo_Sa = vertcat(temp,temp1);
        temp2 = sum(abs(psuedo_Sa-cdf_req').^2*dIM).^(1/2);
        Div_C(ii) = (temp2);
end
bayes_approx = prob_ex_approx.*dsa1./sum(prob_ex_approx.*dsa1*dIM);
bayes_approx = bayes_approx/trapz(pga_val,bayes_approx);
bayes_exact = cdf_req.*dsa1./sum(cdf_req.*dsa1*dIM);
bayes_exact = bayes_exact/trapz(pga_val,bayes_exact);
apr_b(:,hh) = bayes_approx;
ex_b(:,hh) = bayes_exact;
apr_cdf(:,hh) = prob_ex_approx;
ex_cdf(:,hh) = cdf_req;
% A_NC(hh) = sum(bayes_NC.*Div_NC*dIM);
% A_C(hh) = sum(bayes_C.*Div_C*dIM);
% H(hh) = (sqrt(0.5*sum((sqrt(bayes_NC)-sqrt(bayes_C)).^2)*dIM));
% B(hh) = -log(sum(sqrt(bayes_NC.*bayes_C)*dIM));
end
%suff = exp(max([A_NC A_C]))/exp(min([A_NC A_C]))
%plot(Dr,div)
%plot(pga_val,prob_ex_approx,pga_val,cdf_req)
% figure
% plot(Dr,div)