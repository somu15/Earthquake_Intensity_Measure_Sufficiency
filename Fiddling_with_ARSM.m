clr
beta_pga = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\Regression details_Medina_Krawinkler\Bilinear_Regression\PFA4\PGA\beta.txt');
beta_sa1 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\Regression details_Medina_Krawinkler\Bilinear_Regression\PFA4\SaT1.33\beta.txt');
beta_sa2 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\Regression details_Medina_Krawinkler\Bilinear_Regression\PFA4\SaT0.43\beta.txt');
beta_sa3 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\Regression details_Medina_Krawinkler\Bilinear_Regression\PFA4\SaT0.22\beta.txt');
beta_sa2s = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\Regression details_Medina_Krawinkler\Bilinear_Regression\PFA4\SaT2\beta.txt');
beta_sa1s = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\Regression details_Medina_Krawinkler\Bilinear_Regression\PFA4\SaT1\beta.txt');
beta_sa4 = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\Regression details_Medina_Krawinkler\Bilinear_Regression\PFA4\SaT1.5\beta.txt');
beta_pgv = importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\Regression details_Medina_Krawinkler\Bilinear_Regression\PFA4\PGV\beta.txt');

file_path  = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Medina Krawinkler records unscaled\PFA4\';

pga = importdata(strcat(file_path,'x_pga.txt'));
sa1 = importdata(strcat(file_path,'x_sa1.txt'));
sa2 = importdata(strcat(file_path,'x_sa2.txt'));
sa3 = importdata(strcat(file_path,'x_sa3.txt'));
sa2s = importdata(strcat(file_path,'x_sa2sec.txt'));
sa1s = importdata(strcat(file_path,'x_sa1sec.txt'));
sa4 = importdata(strcat(file_path,'x_sa1.5sec.txt'));
pgv = importdata(strcat(file_path,'x_pgv.txt'));
Dr = importdata(strcat(file_path,'y.txt'));
for ii = 1:max(size(Dr))
p_pga(ii) = normpdf(log(Dr(ii)),beta_pga(1)+beta_pga(2)*log(pga(ii)),beta_pga(3));
p_sa1(ii) = normpdf(log(Dr(ii)),beta_sa1(1)+beta_sa1(2)*log(sa1(ii)),beta_sa1(3));
p_sa2(ii) = normpdf(log(Dr(ii)),beta_sa2(1)+beta_sa2(2)*log(sa2(ii)),beta_sa2(3));
p_sa3(ii) = normpdf(log(Dr(ii)),beta_sa3(1)+beta_sa3(2)*log(sa3(ii)),beta_sa3(3));
p_sa2s(ii) = normpdf(log(Dr(ii)),beta_sa2s(1)+beta_sa2s(2)*log(sa2s(ii)),beta_sa2s(3));
p_sa1s(ii) = normpdf(log(Dr(ii)),beta_sa1s(1)+beta_sa1s(2)*log(sa1s(ii)),beta_sa1s(3));
p_sa4(ii) = normpdf(log(Dr(ii)),beta_sa4(1)+beta_sa4(2)*log(sa4(ii)),beta_sa4(3));
p_pgv(ii) = normpdf(log(Dr(ii)),beta_pgv(1)+beta_pgv(2)*log(pgv(ii)),beta_pgv(3));
end
res_sa1 = ((log(Dr)-beta_sa1(1)-beta_sa1(2)*log(sa1)));
res_sa4 = ((log(Dr)-beta_sa4(1)-beta_sa4(2)*log(sa4)));
ARSM_sa1 = sum(log2(p_sa4./p_sa1))/max(size(Dr));
ARSM_pga = sum(log2(p_sa4./p_pga))/max(size(Dr));
ARSM_sa2 = sum(log2(p_sa4./p_sa2))/max(size(Dr));
ARSM_sa3 = sum(log2(p_sa4./p_sa3))/max(size(Dr));
ARSM_sa2s = sum(log2(p_sa4./p_sa2s))/max(size(Dr));
ARSM_sa1s = sum(log2(p_sa4./p_sa1s))/max(size(Dr));
ARSM_sa4 = sum(log2(p_sa4./p_sa4))/max(size(Dr));
stds = [beta_sa1(3) beta_sa2(3) beta_sa3(3) beta_sa2s(3) beta_pga(3) beta_sa1s(3) beta_sa4(3) beta_pgv(3)];
std_ratio = stds/min(stds);
std_ratio = std_ratio';
ARSM = [ARSM_sa1,ARSM_sa2,ARSM_sa3,ARSM_sa2s,ARSM_pga,ARSM_sa1s,ARSM_sa4];