clr
file = 'C:\Users\lakshd5\Dropbox\Bayesian_IM_selection\Accounting for heteroskedasticity\Final analysis 05_28_2016\Exact measure for sufficiency and approximations\Deaggregation\Real site\SaT0.2_finer\';
IM = [0.01 0.25 0.5 0.75 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3];% 3.25 3.5 3.75 4 4.25 4.5 4.75 5];
NR = 50;
NM = 50;

for ii = 1:max(size(IM))
temp = importdata(strcat(file,num2str(IM(ii)),'.txt'));
eval(sprintf('deagg%d = temp',ii));
end
temp = 1:max(size(deagg1))/NR:max(size(deagg1));
R = deagg1(temp,1);
M = deagg1(1:NM,2);
for ii = 1:max(size(IM))
eval(sprintf('temp = deagg%d',ii));
temp = temp(:,3:10);
eval(sprintf('deagg%d = temp',ii));
end
for ii = 1:max(size(IM))
    eval(sprintf('temp = deagg%d',ii));
for jj = 1:8
temp1(:,:,jj) = reshape(temp(:,jj),[NM,NR]);
end
eval(sprintf('deagg%d = temp1',ii));
end
 for jj = 1:max(size(IM))
     eval(sprintf('temp = deagg%d',jj));
 for ii = 1:NM
 temp2(ii) = sum(sum(temp(ii,:,:)));
 end
 eval(sprintf('fM%d = temp2',jj))
 end
 for jj = 1:max(size(IM))
     eval(sprintf('temp = deagg%d',jj));
 for ii = 1:NR
 temp3(ii) = sum(sum(temp(:,ii,:)));
 end
 eval(sprintf('fR%d = temp3',jj))
 end
 for jj = 1:max(size(IM))
     eval(sprintf('temp = deagg%d',jj));
 for ii = 1:8
 temp4(ii) = sum(sum(temp(:,:,ii)));
 end
 eval(sprintf('fe%d = temp4',jj))
 end
 for ii = 1:max(size(IM))
 eval(sprintf('femain(:,ii) = fe%d',ii));
 end
 for ii = 1:max(size(IM))
 eval(sprintf('fMmain(:,ii) = fM%d',ii));
 end
 for ii = 1:max(size(IM))
 eval(sprintf('fRmain(:,ii) = fR%d',ii));
 end
 for ii = 1:NM
 coeff = coeffvalues(fit(IM',fMmain(ii,:)','exp2'));
temp = coeff(1)*exp(coeff(2)*IM)+coeff(3)*exp(coeff(4)*IM);
 temp(temp < 0) = 0;
PM(:,ii) = temp;
 end
 for ii = 1:NR
 coeff = coeffvalues(fit(IM',fRmain(ii,:)','exp2'));
temp = coeff(1)*exp(coeff(2)*IM)+coeff(3)*exp(coeff(4)*IM);
 temp(temp < 0) = 0;
PR(:,ii) = temp;
 end
 for ii = 1:8
 coeff = coeffvalues(fit(IM',femain(ii,:)','exp2'));
temp = coeff(1)*exp(coeff(2)*IM)+coeff(3)*exp(coeff(4)*IM);
 temp(temp < 0) = 0;
Pe(:,ii) = temp;
 end
 PR = PR';
 Pe = Pe';
 PM = PM';
 R_M = 1-sum(sum((fMmain-PM).^2))/sum(sum((mean(mean(fMmain))-fMmain).^2));
 R_R = 1-sum(sum((fRmain-PR).^2))/sum(sum((mean(mean(fRmain))-fRmain).^2));
 R_e = 1-sum(sum((femain-Pe).^2))/sum(sum((mean(mean(femain))-femain).^2));
 [xxm,yym] = meshgrid(M,IM);
 figure(1)
 h = surf(xxm',yym',PM);
%set(h,'LineStyle','none')
xlabel('Magnitude')
ylabel('IM (g)')
zlabel('Probability')
hold on
plot3(xxm,yym,fMmain,['o','red'],'linewidth',2)
figure(2)
[xxr,yyr] = meshgrid(R,IM);
 h = surf(xxr',yyr',PR);
%set(h,'LineStyle','none')
xlabel('Distance')
ylabel('IM (g)')
zlabel('Probability')
hold on
plot3(xxr,yyr,fRmain,['o','red'],'linewidth',2)
figure(3)
[xxe,yye] = meshgrid(1:8,IM);
 h = surf(xxe',yye',Pe);
%set(h,'LineStyle','none')
xlabel('Index')
ylabel('IM (g)')
zlabel('Probability')
hold on
plot3(xxe,yye,femain,['o','red'],'linewidth',2)