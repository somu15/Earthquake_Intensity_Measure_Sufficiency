%% Sufficiency
clr
Dr = 10:100:772;
pga_val = 0.05:0.5:15;
reg_full = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Regression details\PFA1\PGA_MRe.txt');
reg = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Regression details\PFA1\PGA.txt');
sno = 2;
Tp = 0.22;
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
fMR = PSHA_BA08_IMSelection_Paper(pga,sno,Tp,0,760);%importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\SaT1.33_Deagg\f(M,R)_',num2str(pga),'g.txt'));
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
        temp2 = sum(abs(psuedo_Sa-prob_ex_approx').^2*0.5).^(1/2);
        Div_approx(ii) = (temp2);
end
div(hh) = sum(abs(prob_ex_approx-cdf_req).^2*0.5).^(1/2);
end
%plot(Dr,div)
plot(pga_val,prob_ex_approx,pga_val,cdf_req)
figure
plot(Dr,div)
%% Efficiency
clr
Dr = 10:100:772;
pga_val = 0.01:0.01:150;
reg = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Regression details\PFA1\PGA.txt');
for hh = 1:max(size(Dr))
prob_ex_approx = 1-normcdf(log(Dr(hh)),reg(2)*log(pga_val)+reg(1),reg(3));
cdf_req = 1-normcdf(log(Dr(hh)),reg(2)*log(pga_val)+reg(1),0.000001);
for ii = 1:max(size(pga_val))
    temp = zeros(ii,1);
    temp1 = ones(max(size(pga_val))-ii,1);
    psuedo_Sa = vertcat(temp,temp1);
        temp2 = sum(abs(psuedo_Sa-prob_ex_approx').^2*0.5).^(1/2);
        Div_approx(ii) = (temp2);
end
div(hh) = sum(abs(prob_ex_approx-cdf_req).^2*0.5).^(1/2);
end
plot(Dr,div)