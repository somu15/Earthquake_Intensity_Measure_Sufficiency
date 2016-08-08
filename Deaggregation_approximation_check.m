clr
file = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\Deaggregation\Hypothetical site\PGV\deagg_';
IM = [1 100 200 300 400 500 600];
index = 1:20;
for ss = 1:max(size(index))
for ii = 1:max(size(IM))
pga = IM(ii);
fMRe = importdata(strcat(file,num2str(pga),'.txt'));
fMRe = reshape(fMRe,[20 20 20]);
for kk = 1:20
 fM(kk) = sum(sum(fMRe(kk,:,:)));
  fR(kk) = sum(sum(fMRe(:,kk,:)));
   fe(kk) = sum(sum(fMRe(:,:,kk)));
end
 Pm(ii) = fM(index(ss));
 Pr(ii) = fR(index(ss));
 Pe(ii) = fe(index(ss));
end
IM1 = 1:600;
% fit_type = fittype( @(a,b,c,d,x) a*exp(b*x)+c*exp(d*x),'independent','x');
coeff_M = coeffvalues(fit(IM',Pm','exp2'));
coeff_R = coeffvalues(fit(IM',Pr','exp2'));
coeff_eps = coeffvalues(fit(IM',Pe','exp2'));
Pm_approx(:,ss) = coeff_M(1)*exp(coeff_M(2)*IM1)+coeff_M(3)*exp(coeff_M(4)*IM1);
Pr_approx(:,ss) = coeff_R(1)*exp(coeff_R(2)*IM1)+coeff_R(3)*exp(coeff_R(4)*IM1);
Pe_approx(:,ss) = coeff_eps(1)*exp(coeff_eps(2)*IM1)+coeff_eps(3)*exp(coeff_eps(4)*IM1);
for ii = 1:max(size(IM1))
pga = IM1(ii);
fMRe = importdata(strcat(file,num2str(pga),'.txt'));
fMRe = reshape(fMRe,[20 20 20]);
for kk = 1:20
 fM(kk) = sum(sum(fMRe(kk,:,:)));
  fR(kk) = sum(sum(fMRe(:,kk,:)));
   fe(kk) = sum(sum(fMRe(:,:,kk)));
end
 Pm1(ii,ss) = fM(index(ss));
 Pr1(ii,ss) = fR(index(ss));
 Pe1(ii,ss) = fe(index(ss));
end
end
[xx,yy] = meshgrid(IM1,index);
figure(1)
h = surf(xx',yy',Pm_approx);
set(h,'LineStyle','none')
hold on
plot3(xx',yy',Pm1,['o','red'])
xlabel('IM (g)')
ylabel('Bin number')
zlabel('Probability')
legend('Approximate','Exact')
title('Magnitude')
figure(2)
h = surf(xx',yy',Pr_approx);
set(h,'LineStyle','none')
hold on
plot3(xx',yy',Pr1,['o','red'])
xlabel('IM (g)')
ylabel('Bin number')
zlabel('Probability')
legend('Approximate','Exact')
title('Distance')
figure(3)
h = surf(xx',yy',Pe_approx);
set(h,'LineStyle','none')
hold on
plot3(xx',yy',Pe1,['o','red'])
xlabel('IM (g)')
ylabel('Bin number')
zlabel('Probability')
legend('Approximate','Exact')
title('Epsilon')