clr
format long
% Data imports
for iii = 1:7
names = {'SaT1.33','SaT0.43','SaT0.22','SaT2','PGA','SaT1','SaT1.5'};
Hazard = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\Hazard curves\Hazard curves with eps\',names{iii},'.txt'));
beta = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\Regression details_FEMA\Bilinear_Regression\IDR1\',names{iii},'\beta.txt'));
beta_all = importdata(strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\Regression details_FEMA\Multilinear_Regression\IDR1\',names{iii},'\beta.txt'));
file = strcat('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\Deaggregation\Hypothetical site\',names{iii},'\deagg_');
Dr = 0.005:0.005:0.04;
if strcmp('PGV',names{iii}) == 1
    IM = 1:1:600;
    dIM = 1;
else
    IM = 0.01:0.01:6;
    dIM = 0.01;
end
upper_index = 600;

% blah blah blah
% Computations
M = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\M.txt');
R = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\R.txt');
eps_max = 4;
intervals = 20;
eps_min = -eps_max;
deps = (eps_max-eps_min)/(intervals-1);
eps = eps_min:deps:eps_max;
dM = M(2)-M(1);
dR = R(2)-R(1);
prob_ex_M = zeros(max(size(IM)),1);prob_ex_R = zeros(max(size(IM)),1);prob_ex_eps = zeros(max(size(IM)),1);
diff_IM = abs((Differentiation(dIM,Hazard(:,2))));
diff_IM = diff_IM(1:upper_index);
for ii = 1:max(size(Dr))
    for jj = 1:max(size(IM))
        pga = IM(jj);
fMRe = importdata(strcat(file,num2str(pga),'.txt'));
fMRe = reshape(fMRe,[20 20 20]);
 for iii = 1:20
     for jjj = 1:20
         for kkk = 1:20
             temp_all(iii,jjj,kkk) = 1-normcdf(log(Dr(ii)),beta_all(1)+beta_all(2)*log(IM(jj))+beta_all(3)*M(iii)+beta_all(4)*R(jjj)+beta_all(5)*eps(kkk),beta_all(6));
         end
     end
 end
 prob_ex_all(jj) = sum(sum(sum(temp_all.*fMRe)));
    end
    prob_ex_approx = 1-normcdf(log(Dr(ii)),beta(1)+beta(2)*log(IM),beta(3));
    bayes(:,ii) = prob_ex_approx.*diff_IM/sum(prob_ex_approx.*diff_IM*dIM);
    bayes(:,ii) = bayes(:,ii)/trapz(IM,bayes(:,ii));
    
    
    
    bayes_all(:,ii) = prob_ex_all.*diff_IM/sum(prob_ex_all.*diff_IM*dIM);
    bayes_all(:,ii) = bayes_all(:,ii)/trapz(IM,bayes_all(:,ii));
    
    temp1 = log2(bayes_all(:,ii)./bayes(:,ii));
    temp1(isnan(temp1)) = 0;
    temp1(isinf(temp1)) = 0;
    KLD_all(ii) = sum(bayes_all(:,ii).*temp1)*dIM;
end
all(:,iii) = KLD_all;
end