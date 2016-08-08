clr
format long
file = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\Deaggregation\Hypothetical site\PGV\deagg_';
Hazard = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\Hazard curves\Hazard curves with eps\PGV.txt');
beta = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\Regression details_Medina_Krawinkler\Bilinear_Regression\PFA4\PGV\beta.txt');
beta_M = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\Regression details_Medina_Krawinkler\Bilinear_Regression\PFA4\PGV\beta_M.txt');
beta_R = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\Regression details_Medina_Krawinkler\Bilinear_Regression\PFA4\PGV\beta_R.txt');
beta_eps = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\Regression details_Medina_Krawinkler\Bilinear_Regression\PFA4\PGV\beta_eps.txt');

Dr = 10:100:772;
IM = [1 100 200 300 400 500 600];
for ii = 1:max(size(IM))
pga = IM(ii);
fMRe = importdata(strcat(file,num2str(pga),'.txt'));
fMRe = reshape(fMRe,[20 20 20]);
for kk = 1:20
 fM(kk) = sum(sum(fMRe(kk,:,:)));
  fR(kk) = sum(sum(fMRe(:,kk,:)));
   fe(kk) = sum(sum(fMRe(:,:,kk)));
end
Pm(:,ii) = fM;
Pr(:,ii) = fR;
Pe(:,ii) = fe;
end
for ii = 1:20
coeff_M(ii,:) = coeffvalues(fit(IM',Pm(ii,:)','exp2'));
coeff_R(ii,:) = coeffvalues(fit(IM',Pr(ii,:)','exp2'));
coeff_eps(ii,:) = coeffvalues(fit(IM',Pe(ii,:)','exp2'));
end
IM = 1:1:600;%%%%%%%%%%%%%%
dIM = 1;
upper_index = 600;
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
temp_M = zeros(20,1);temp_R = zeros(20,1);temp_eps = zeros(20,1);
 for kk = 1:20
 fM(kk) = coeff_M(kk,1)*exp(coeff_M(kk,2)*pga)+coeff_M(kk,3)*exp(coeff_M(kk,4)*pga);
   fR(kk) = coeff_R(kk,1)*exp(coeff_R(kk,2)*pga)+coeff_R(kk,3)*exp(coeff_R(kk,4)*pga);
   fe(kk) = coeff_eps(kk,1)*exp(coeff_eps(kk,2)*pga)+coeff_eps(kk,3)*exp(coeff_eps(kk,4)*pga);
   temp_M(kk) = 1-normcdf(log(Dr(ii)),beta_M(1)+beta_M(2)*log(IM(jj))+beta_M(3)*M(kk),beta_M(4));
 temp_R(kk) = 1-normcdf(log(Dr(ii)),beta_R(1)+beta_R(2)*log(IM(jj))+beta_R(3)*R(kk),beta_R(4));
 temp_eps(kk) = 1-normcdf(log(Dr(ii)),beta_eps(1)+beta_eps(2)*log(IM(jj))+beta_eps(3)*eps(kk),beta_eps(4));
 end
 fM = fM/sum(fM);
fR = fR/sum(fR);
fe = fe/sum(fe);
 prob_ex_M(jj) = sum(temp_M.*fM');
 prob_ex_R(jj) = sum(temp_R.*fR');
 prob_ex_eps(jj) = sum(temp_eps.*fe');
    end
    prob_ex_approx = 1-normcdf(log(Dr(ii)),beta(1)+beta(2)*log(IM),beta(3));
    prob = normpdf(log(Dr(ii)),beta(1)+beta(2)*log(IM),beta(3));
    bayes = prob_ex_approx.*diff_IM/sum(prob_ex_approx.*diff_IM*dIM);
    bayes = bayes/trapz(IM,bayes);
    joint(:,ii) = bayes;
    bayes_M = prob_ex_M.*diff_IM'/sum(prob_ex_M.*diff_IM'*dIM);
    bayes_M = bayes_M/trapz(IM,bayes_M);
    joint_M(:,ii) = bayes_M;
    bayes_R = prob_ex_R.*diff_IM'/sum(prob_ex_R.*diff_IM'*dIM);
    bayes_R = bayes_R/trapz(IM,bayes_R);
    joint_R(:,ii) = bayes_R;
    bayes_eps = prob_ex_eps.*diff_IM'/sum(prob_ex_eps.*diff_IM'*dIM);
    bayes_eps = bayes_eps/trapz(IM,bayes_eps);
    joint_eps(:,ii) = bayes_eps;
    % Compute the information gain
    temp1 = log2(bayes_M./bayes');
    temp1(isnan(temp1)) = 0;
    temp1(isinf(temp1)) = 0;
    KLD_M(ii) = sum(bayes_M.*temp1)*dIM;
    
    
    temp1 = log2(bayes_R./bayes');
    temp1(isnan(temp1)) = 0;
    temp1(isinf(temp1)) = 0;
    KLD_R(ii) = sum(bayes_R.*temp1)*dIM;
    
    
    temp1 = log2(bayes_eps./bayes');
    temp1(isnan(temp1)) = 0;
    temp1(isinf(temp1)) = 0;
    KLD_eps(ii) = sum(bayes_eps.*temp1)*dIM;
end
KLD_eps(KLD_eps<0) = 0;
KLD_R(KLD_R<0) = 0;
KLD_M(KLD_M<0) = 0;
KLD = KLD_M + KLD_R + KLD_eps;
KLD = real(KLD);
plot(Dr,KLD)