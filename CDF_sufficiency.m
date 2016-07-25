clr
format long
% Data imports
Hazard = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Hazard curves\Hazard curves with eps\PGA.txt');
beta = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Regression details_Medina_Krawinkler\Bilinear_Regression\RD\PGA\beta.txt');
beta_M = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Regression details_Medina_Krawinkler\Bilinear_Regression\RD\PGA\beta_M.txt');
beta_R = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Regression details_Medina_Krawinkler\Bilinear_Regression\RD\PGA\beta_R.txt');
beta_eps = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Regression details_Medina_Krawinkler\Bilinear_Regression\RD\PGA\beta_eps.txt');
file = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Deaggregation\PGA\deagg_';
Dr = 0.005:0.005:0.04;
IM = 0.01:0.01:10;
dIM = 0.01;
upper_index = 1000;


% Computations
M = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\M.txt');
R = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\R.txt');
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
fM = zeros(20,1); fR = zeros(20,1); fe = zeros(20,1);
temp_M = zeros(20,1);temp_R = zeros(20,1);temp_eps = zeros(20,1);
 for kk = 1:20
 fM(kk) = sum(sum(fMRe(kk,:,:)));
  fR(kk) = sum(sum(fMRe(:,kk,:)));
   fe(kk) = sum(sum(fMRe(:,:,kk)));
   temp_M(kk) = 1-normcdf(log(Dr(ii)),beta_M(1)+beta_M(2)*log(IM(jj))+beta_M(3)*M(kk),beta_M(4));
 temp_R(kk) = 1-normcdf(log(Dr(ii)),beta_R(1)+beta_R(2)*log(IM(jj))+beta_R(3)*R(kk),beta_R(4));
 temp_eps(kk) = 1-normcdf(log(Dr(ii)),beta_eps(1)+beta_eps(2)*log(IM(jj))+beta_eps(3)*eps(kk),beta_eps(4));
 end
 prob_ex_M(jj) = sum(temp_M.*fM);
 prob_ex_R(jj) = sum(temp_R.*fR);
 prob_ex_eps(jj) = sum(temp_eps.*fe);
    end
    prob_ex_approx = 1-normcdf(log(Dr(ii)),beta(1)+beta(2)*log(IM),beta(3));
    
    bayes = prob_ex_approx.*diff_IM/sum(prob_ex_approx.*diff_IM*dIM);
    bayes = bayes/trapz(IM,bayes);
    
    bayes_M = prob_ex_M.*diff_IM'/sum(prob_ex_M.*diff_IM'*dIM);
    bayes_M = bayes_M/trapz(IM,bayes_M);
    
    bayes_R = prob_ex_R.*diff_IM'/sum(prob_ex_R.*diff_IM'*dIM);
    bayes_R = bayes_R/trapz(IM,bayes_R);
    
    bayes_eps = prob_ex_eps.*diff_IM'/sum(prob_ex_eps.*diff_IM'*dIM);
    bayes_eps = bayes_eps/trapz(IM,bayes_eps);
    
    for jj = 1:max(size(IM))
    temp = zeros(jj,1);
    temp1 = ones(max(size(IM))-jj,1);
    psuedo_Sa = vertcat(temp,temp1);
        Div(jj) = sum(abs(psuedo_Sa-prob_ex_approx').^2*dIM).^(1/2);
        Div_M(jj) = sum(abs(psuedo_Sa-prob_ex_M).^2*dIM).^(1/2);
        Div_R(jj) = sum(abs(psuedo_Sa-prob_ex_R).^2*dIM).^(1/2);
        Div_eps(jj) = sum(abs(psuedo_Sa-prob_ex_eps).^2*dIM).^(1/2);
    end
    A(ii) = sum(Div.*bayes*dIM);
    A_M(ii) = sum(Div_M.*bayes_M'*dIM);
    A_R(ii) = sum(Div_R.*bayes_R'*dIM);
    A_eps(ii) = sum(Div_eps.*bayes_eps'*dIM);
    %sum(bayes_Sa./Div_Sa*dIM);
end
m_M = abs(A-A_M);
m_R = abs(A-A_R);
m_eps = abs(A-A_eps);
%Div = Div_M+Div_R+Div_eps;
m = m_M+m_R+m_eps;