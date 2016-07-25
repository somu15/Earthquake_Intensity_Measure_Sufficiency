clr
format long
% Data imports
Hazard = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Hazard curves\Hazard curves with eps\SaT1.33.txt');
beta = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Regression details_FEMA\Bilinear_Regression\RD\SaT1.33\beta.txt');
beta_M = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Regression details_FEMA\Bilinear_Regression\RD\SaT1.33\beta_M.txt');
beta_R = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Regression details_FEMA\Bilinear_Regression\RD\SaT1.33\beta_R.txt');
beta_eps = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Regression details_FEMA\Bilinear_Regression\RD\SaT1.33\beta_eps.txt');
file = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency\Deaggregation\SaT1.33\deagg_';
Dr = 0.001:0.001:0.04;
IM = 0.01:0.01:6;
dIM = 0.01;
upper_index = 600;

% blah blah blah
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
    prob = normpdf(log(Dr(ii)),beta(1)+beta(2)*log(IM),beta(3));
    bayes(:,ii) = prob_ex_approx.*diff_IM/sum(prob_ex_approx.*diff_IM*dIM);
    bayes(:,ii) = bayes(:,ii)/trapz(IM,bayes(:,ii));
    %joint(:,ii) = bayes;
    bayes_M(:,ii) = prob_ex_M.*diff_IM'/sum(prob_ex_M.*diff_IM'*dIM);
    bayes_M(:,ii) = bayes_M(:,ii)/trapz(IM,bayes_M(:,ii));
    %joint_M(:,ii) = bayes_M;
    bayes_R(:,ii) = prob_ex_R.*diff_IM'/sum(prob_ex_R.*diff_IM'*dIM);
    bayes_R(:,ii) = bayes_R(:,ii)/trapz(IM,bayes_R(:,ii));
    %joint_R(:,ii) = bayes_R;
    bayes_eps(:,ii) = prob_ex_eps.*diff_IM'/sum(prob_ex_eps.*diff_IM'*dIM);
    bayes_eps(:,ii) = bayes_eps(:,ii)/trapz(IM,bayes_eps(:,ii));
    %joint_eps(:,ii) = bayes_eps;
    % Compute the information gain
    temp1 = log2(bayes_M(:,ii)./bayes(:,ii));
    temp1(isnan(temp1)) = 0;
    temp1(isinf(temp1)) = 0;
    KLD_M(ii) = sum(bayes_M(:,ii).*temp1)*dIM;
    
    
    temp1 = log2(bayes_R(:,ii)./bayes(:,ii));
    temp1(isnan(temp1)) = 0;
    temp1(isinf(temp1)) = 0;
    KLD_R(ii) = sum(bayes_R(:,ii).*temp1)*dIM;
    
    
    temp1 = log2(bayes_eps(:,ii)./bayes(:,ii));
    temp1(isnan(temp1)) = 0;
    temp1(isinf(temp1)) = 0;
    KLD_eps(ii) = sum(bayes_eps(:,ii).*temp1)*dIM;
end
KLD_eps(KLD_eps<0) = 0;
KLD_R(KLD_R<0) = 0;
KLD_M(KLD_M<0) = 0;
KLD = KLD_M + KLD_R + KLD_eps;
plot(Dr,KLD)
%plot(IM,bayes,'red',IM,bayes_M,'blue',IM,bayes_R,'green',IM,bayes_eps,'linewidth',1.5)
% subplot(1,2,1)
% plot(IM,bayes(:,1),'red',IM,bayes_M(:,1),'blue',IM,bayes_R(:,1),'green',IM,bayes_eps(:,1),'linewidth',1.5)
% xlabel('IM')
% ylabel('f(IM | Dr>0.005)')
% title('Roof drift = 0.005')
% subplot(1,2,2)
% plot(IM,bayes(:,2),'red',IM,bayes_M(:,2),'blue',IM,bayes_R(:,2),'green',IM,bayes_eps(:,2),'linewidth',1.5)
% xlabel('IM')
% ylabel('f(IM | Dr>0.04)')
% title('Roof drift = 0.04')
% legend('Only IM','IM & M','IM & R','IM & eps')