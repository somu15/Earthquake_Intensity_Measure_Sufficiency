%% ROOF DRIFT
clr
Dr = 350;

%    subplot(1,2,1)
%    loglog(Sa,AFE_sa1,'blue','linewidth',4)
%    set(gca, 'FontName', 'Times')
%    set(gca, 'FontSize', 30)
%    xlabel('IM*')
%    ylabel('F(IM*)')
%    grid on
%    subplot(1,2,2)
%    yyaxis left
%    plot(Sa,prob_ex_Sa,['blue'],'linewidth',4)
%    set(gca, 'FontName', 'Times')
%    set(gca, 'FontSize', 30)
%    xlabel('IM*')
%    ylabel('Pr.(thetamax > x | IM*)')
%    grid on
%    yyaxis right
%    plot(Sa,bayes_Sa,['red','--'],'linewidth',4)
%    ylabel('w*')
%    legend('CDF','weight function')
%    hold on
%    scatter(0.22,0,'g','linewidth',10)
%   subplot(3,2,1)
%   loglog(Sa,AFE_sa1,'blue','linewidth',4)
% set(gca, 'FontName', 'Times')
% set(gca, 'FontSize', 25)
%   set(gca,'YTick',[])
%   set(gca,'XTick',[])
%   xlabel('IM')
%   ylabel('Pr.(IM > y)')
%   title('(a)')
%   file_path  = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\data\Roof Drift\';
% Sa1 = importdata(strcat(file_path,'xtemp.txt'));
% Dr = importdata(strcat(file_path,'ytemp.txt'));
%   subplot(3,2,2)
%   plot(log(Sa),(0.6702*log(Sa)-3.6251),'blue','linewidth',4)
%   hold on
%   scatter(log(Sa1),log(Dr),'linewidth',2)
%   set(gca, 'FontName', 'Times')
% set(gca, 'FontSize', 25)
%   set(gca,'YTick',[])
%   set(gca,'XTick',[])
%   xlabel('ln(IM)')
%   ylabel('ln(theta)')
%   title('(b)')
%   subplot(3,2,3)
%   plot(Sa,prob_ex_Sa,'blue','linewidth',4)
%   set(gca,'YTick',[])
%   set(gca,'XTick',[])
%   set(gca, 'FontName', 'Times')
% set(gca, 'FontSize', 25)
%   xlabel('IM')
%   ylabel('Pr.(theta>x|IM)')
%   title('(c)')
%   subplot(3,2,4)
%   plot(Sa,bayes_Sa,'blue','linewidth',4)
%   set(gca,'YTick',[])
%   set(gca,'XTick',[])
%   set(gca, 'FontName', 'Times')
% set(gca, 'FontSize', 25)
%   xlabel('IM')
%   ylabel('w(IM)')
%   title('(d)')
%   subplot(3,2,[5,6])
%   plot(Sa,prob_ex_Sa,'blue','linewidth',4)
%   hold on
%   temp = zeros(100,1);
%     temp1 = ones(max(size(Sa))-100,1);
%     psuedo_Sa = vertcat(temp,temp1);
%   plot(Sa,psuedo_Sa,['red','--'],'linewidth',4)
%   set(gca,'YTick',[])
%   set(gca,'XTick',[])
%   set(gca, 'FontName', 'Times')
% set(gca, 'FontSize', 25)
%   xlabel('IM')
%   ylabel('Pr.(theta>x|IM)')
%   title('(e)')
%   legend('Realistic IM','Psuedo IM
%%
% PGA
pga = 0.01:0.01:10;
pga_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_PGA_hazard.txt'));
AFE_pga = exp(interp1(log(pga_temp(:,1)),log(pga_temp(:,2)),log(pga),'spline'));
dpga = abs((Differentiation((0.01),AFE_pga)));
std_pga = 0.1873*ones(1,max(size(pga)));
prob_ex_pga = 1-normcdf((log(Dr)*ones(max(size(pga)),1)'-(0.6683*log(pga)+6.2827))./std_pga,0,1);
bayes_pga = prob_ex_pga.*dpga./sum(prob_ex_pga.*dpga*0.01);
   for ii = 1:max(size(pga))
    temp = zeros(ii,1);
    temp1 = ones(max(size(pga))-ii,1);
    psuedo_pga = vertcat(temp,temp1);
    temp2 = sum(abs(psuedo_pga-prob_ex_pga').^2*0.01).^0.5;
    Div_pga(ii) = (temp2);
   end
   A_pga = sum(bayes_pga./Div_pga*0.01)
  % PGV
pgv = 1:1:1000;
pgv_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_PGV_hazard.txt'));
AFE_pgv = exp(interp1(log(pgv_temp(:,1)),log(pgv_temp(:,2)),log(pgv),'spline'));
dpgv = abs((Differentiation((0.01),AFE_pgv)));
std_pgv = 0.2836*ones(1,max(size(pgv)));
prob_ex_pgv = 1-normcdf((log(Dr)*ones(max(size(pgv)),1)'-(0.3528*log(pgv)+4.3002))./std_pgv,0,1);
bayes_pgv = prob_ex_pgv.*dpgv./sum(prob_ex_pgv.*dpgv*1);
   for ii = 1:max(size(pgv))
    temp = zeros(ii,1);
    temp1 = ones(max(size(pgv))-ii,1);
    psuedo_pgv = vertcat(temp,temp1);
    temp2 = sum(abs(psuedo_pgv-prob_ex_pgv').^2*1).^0.5;
    Div_pgv(ii) = (temp2);
   end
   A_pgv = sum(bayes_pgv./Div_pgv*1)
   
   % PGD
pgd = 1:1:300;
pgd_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_PGD_hazard.txt'));
AFE_pgd = exp(interp1(log(pgd_temp(:,1)),log(pgd_temp(:,2)),log(pgd),'spline'));
dpgd = abs((Differentiation((0.01),AFE_pgd)));
std_pgd = 0.2934*ones(1,max(size(pgd)));
prob_ex_pgd = 1-normcdf((log(Dr)*ones(max(size(pgd)),1)'-(-0.1696*log(pgd)+6.0913))./std_pgd,0,1);
bayes_pgd = prob_ex_pgd.*dpgd./sum(prob_ex_pgd.*dpgd*1);
   for ii = 1:max(size(pgd))
    temp = zeros(ii,1);
    temp1 = ones(max(size(pgd))-ii,1);
    psuedo_pgd = vertcat(temp,temp1);
    temp2 = sum(abs(psuedo_pgd-prob_ex_pgd').^2*1).^0.5;
    Div_pgd(ii) = (temp2);
   end
   A_pgd = sum(bayes_pgd./Div_pgd*1)
   
   % Sa1
Sa = 0.01:0.01:10;
sa1_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_SaT1.33.txt'));
AFE_sa1 = exp(interp1(log(sa1_temp(:,1)),log(sa1_temp(:,2)),log(Sa),'spline'));
dsa1 = abs((Differentiation((0.01),AFE_sa1)));
std_Sa = 0.296*ones(1,max(size(Sa)));
prob_ex_Sa = 1-normcdf((log(Dr)*ones(max(size(Sa)),1)'-(0.1874*log(Sa)+5.8416))./std_Sa,0,1);
bayes_Sa = prob_ex_Sa.*dsa1./sum(prob_ex_Sa.*dsa1*0.01);
   for ii = 1:max(size(Sa))
    temp = zeros(ii,1);
    temp1 = ones(max(size(Sa))-ii,1);
    psuedo_Sa = vertcat(temp,temp1);
    temp2 = sum(abs(psuedo_Sa-prob_ex_Sa').^2*0.01).^0.5;
    Div_Sa(ii) = (temp2);
   end
   A_Sa1 = sum(bayes_Sa./Div_Sa*0.01)
   
   % Sa2
Sa = 0.01:0.01:10;
sa2_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_Sa0.433_hazard.txt'));
AFE_sa2 = exp(interp1(log(sa2_temp(:,1)),log(sa2_temp(:,2)),log(Sa),'spline'));
dsa2 = abs((Differentiation((0.01),AFE_sa2)));
std_Sa2 = 0.1751*ones(1,max(size(Sa)));
prob_ex_Sa2 = 1-normcdf((log(Dr)*ones(max(size(Sa)),1)'-(0.4932*log(Sa)+5.7396))./std_Sa2,0,1);
bayes_Sa2 = prob_ex_Sa2.*dsa2./sum(prob_ex_Sa2.*dsa2*0.01);
   for ii = 1:max(size(Sa))
    temp = zeros(ii,1);
    temp1 = ones(max(size(Sa))-ii,1);
    psuedo_Sa2 = vertcat(temp,temp1);
    temp2 = sum(abs(psuedo_Sa2-prob_ex_Sa2').^2*0.01).^0.5;
    Div_Sa2(ii) = (temp2);
   end
   A_Sa2 = sum(bayes_Sa2./Div_Sa2*0.01)
   % Sa3
Sa = 0.01:0.01:10;
sa3_temp = (importdata('C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\OpenSHA_Sa0.22_hazard.txt'));
AFE_sa3 = exp(interp1(log(sa3_temp(:,1)),log(sa3_temp(:,2)),log(Sa),'spline'));
dsa3 = abs((Differentiation((0.01),AFE_sa3)));
std_Sa3 = 0.2171*ones(1,max(size(Sa)));
prob_ex_Sa3 = 1-normcdf((log(Dr)*ones(max(size(Sa)),1)'-(0.5087*log(Sa)+5.7111))./std_Sa3,0,1);
bayes_Sa3 = prob_ex_Sa3.*dsa3./sum(prob_ex_Sa3.*dsa3*0.01);
   for ii = 1:max(size(Sa))
    temp = zeros(ii,1);
    temp1 = ones(max(size(Sa))-ii,1);
    psuedo_Sa3 = vertcat(temp,temp1);
    temp2 = sum(abs(psuedo_Sa3-prob_ex_Sa3').^2*0.01).^0.5;
    Div_Sa3(ii) = (temp2);
   end
   A_Sa3 = sum(bayes_Sa3./Div_Sa3*0.01)
  