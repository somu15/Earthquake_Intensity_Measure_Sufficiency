%% A PSHA tool
% Written by Somayajulu on September 21, 2014
%function [] = PSHA()
% IN AREA SOURCES ONLY RECTANGULAR FAULTS ARE CONSIDERED %
%GLOBAL AXES IS ASSUMED TO BE AT THE SITE %
clr
format long
display('Hi')
global BA08_check
BA08_check = 1;
%%%%%%%%   USER INPUTS %%%%%%%%%%%%%%%%%


prompt = ' How many sources of earthquake to be considered? (NO AREA SOURCES PLEASE)';
Nsources = input(prompt);      % NUMBER OF EARTHQUAKE SOURCES %
prompt = ' of which point sources are? ';
Npoint = input(prompt);    % NUMBER OF POINT SOURCES %
prompt = ' of which line sources are? ';
Nline = input(prompt);      % NUMBER OF LINE SOURCES %
prompt = ' of which area sources are? (PLEASE ENTER ONLY A ZERO HERE AS AREA SOURCES PART OF THE CODE IS NOT FULLY READY) ';
Narea = input(prompt);          % NUMBER OF AREA SOURCES %
if Nsources ~= Npoint+Nline+Narea
    display(' INVALID INPUT, CHECK THE NUMBER OF SOURCES GIVEN ')    % CHECK WHETHER THE NUMBER TOTAL SOURCES ARE EQUAL TO SUM OF INDIVIDUAL SOURCES %
    return
end
if Npoint>0
    for i = 1:Npoint
        prompt = sprintf(' Enter the X co-ordinate of point source %d ',i);    % CO-ORDINATES OF POINT SOURCE %
        eval(sprintf('Xpoint%d = input(prompt)',i));
        prompt = sprintf(' Enter the Y co-ordinate of point source %d ',i);
        eval(sprintf('Ypoint%d = input(prompt)',i)); 
    end
end
if Nline>0
    display(' The geometry of the line source is characterized by the co-ordinates of the points at the ends. ')
    for i = 1:Nline
        prompt = sprintf(' Enter the X co-ordinate of line source %d at one end',i);    % CO-ORDINATES OF LINE SOURCE %
        eval(sprintf('Xlineone%d = input(prompt)',i));
        prompt = sprintf(' Enter the Y co-ordinate of line source %d at one end',i);
        eval(sprintf('Ylineone%d = input(prompt)',i));
        prompt = sprintf(' Enter the X co-ordinate of line source %d at other end',i);
        eval(sprintf('Xlineother%d = input(prompt)',i));
        prompt = sprintf(' Enter the Y co-ordinate of line source %d at other end',i);
        eval(sprintf('Ylineother%d = input(prompt)',i));
    end
end
if Narea>0
    display(' The geometry of the area source is characterized by the co-ordinates of points at the corners of the rectangular fault ')
    for i = 1:Narea
        prompt = sprintf(' Enter the X co-ordinate of area source %d at first corner',i);    % CO-ORDINATES OF AREA SOURCE %
        eval(sprintf('Xareaone%d = input(prompt)',i));
        prompt = sprintf(' Enter the Y co-ordinate of area source %d at first corner',i);
        eval(sprintf('Yareaone%d = input(prompt)',i));
        prompt = sprintf(' Enter the X co-ordinate of area source %d at second corner',i);    
        eval(sprintf('Xareatwo%d = input(prompt)',i));
        prompt = sprintf(' Enter the Y co-ordinate of area source %d at second corner',i);
        eval(sprintf('Yareatwo%d = input(prompt)',i));
        prompt = sprintf(' Enter the X co-ordinate of area source %d at third corner',i);    
        eval(sprintf('Xareathree%d = input(prompt)',i));
        prompt = sprintf(' Enter the Y co-ordinate of area source %d at third corner',i);
        eval(sprintf('Yareathree%d = input(prompt)',i));
        prompt = sprintf(' Enter the X co-ordinate of area source %d at fourth corner',i);    
        eval(sprintf('Xareafour%d = input(prompt)',i));
        prompt = sprintf(' Enter the Y co-ordinate of area source %d at fourth corner',i);
        eval(sprintf('Yareafour%d = input(prompt)',i));
    end
end
for i = 1:Narea
    eval(sprintf('x1 = Xareaone%d',i))
    eval(sprintf('y1 = Yareaone%d',i))
    eval(sprintf('x2 = Xareatwo%d',i))
    eval(sprintf('y2 = Yareatwo%d',i))
    eval(sprintf('x3 = Xareathree%d',i))
    eval(sprintf('y3 = Yareathree%d',i))
    eval(sprintf('x4 = Xareafour%d',i))
    eval(sprintf('y4 = Yareafour%d',i))
    m1 = (y2-y1)/(x2-x1);
    m2 = (y3-y2)/(x3-x2);
    m3 = (y4-y3)/(x4-x3);
    m4 = (y1-y4)/(x1-x4);
    if abs(m1-m2) < 0.1
        m = abs(m3-m4);
        if m < 0.1
        display(' Co-ordinates given for the rectangular fault projection are OK ')
        break
    else
        display(' Check the co-ordinates given for rectangular fault projection ')
        return                                                                                   % FOOLPROOF CHECKING OF THE CO-ORDINATES OF RECTANGULAR FAULT %
        end                                                                        %%%%%% ATTENTION NEEDED HERE %%%%%%%%
    end
    
    if m1*m2 < 0
        m = m3*m4;
        if m < 0
        display(' Co-ordinates given for the rectangular fault projection are OK ')
        break
    else
        display(' Check the co-ordinates given for rectangular fault projection ')
        return
        end
    end
    if isnan(m1*m2)
        m = m3*m4;
        if isnan(m)
        display(' Co-ordinates given for the rectangular fault projection are OK ')
        break
    else
        display(' Check the co-ordinates given for rectangular fault projection ')
        return
        end
    end
    if isnan(m1*m3)
        m = m2*m4;
        if isnan(m)
        display(' Co-ordinates given for the rectangular fault projection are OK ')
        break
    else
        display(' Check the co-ordinates given for rectangular fault projection ')
        return
        end
    end  
    display(' Check the co-ordinates given for rectangular fault projection ')
    return
end
for i = 1:Npoint
    prompt = sprintf('Enter the maximum magnitude of point source %d',i);
    eval(sprintf('Mmaxp%d = input(prompt)',i));
    prompt = sprintf('Enter the minimum magnitude of point source %d',i);
    eval(sprintf('Mminp%d = input(prompt)',i));
end
for i = 1:Nline
    prompt = sprintf('Enter the maximum magnitude of line source %d',i);
    eval(sprintf('Mmaxl%d = input(prompt)',i));
    prompt = sprintf('Enter the minimum magnitude of line source %d',i);                 % ENTER THE MAXIMUM AND MINIMUM MAGNITUDES %
    eval(sprintf('Mminl%d = input(prompt)',i));
end
for i = 1:Narea
    prompt = sprintf('Enter the maximum magnitude of area source %d',i);
    eval(sprintf('Mmaxa%d = input(prompt)',i));
    prompt = sprintf('Enter the minimum magnitude of area source %d',i);
    eval(sprintf('Mmina%d = input(prompt)',i));
end
display('For calculation of uncertainities in size the Gutenberg-Richter recurrence relation is used in this code')
display('The Gutenberg-Richter recurrence relation is of the form log(lambda) = a-b*m  ***(BASE 10)***')
for i = 1:Npoint
    prompt = sprintf('Enter the a value in Gutenberg-Richter recurrence relation of point source %d',i);
    eval(sprintf('aP%d = input(prompt)',i));
    prompt = sprintf('Enter the b value in Gutenberg-Richter recurrence relation of point source %d',i);
    eval(sprintf('bP%d = input(prompt)',i));
end
for i = 1:Nline
    prompt = sprintf('Enter the a value in Gutenberg-Richter recurrence relation of line source %d',i);
    eval(sprintf('aL%d = input(prompt)',i));
    prompt = sprintf('Enter the b value in Gutenberg-Richter recurrence relation of line source %d',i);            % VALUES OF A AND B IN GUTENBERG-RICHTER RECURRENCE RELATION %
    eval(sprintf('bL%d = input(prompt)',i));
end
for i = 1:Narea
    prompt = sprintf('Enter the a value in Gutenberg-Richter recurrence relation of area source %d',i);
    eval(sprintf('aA%d = input(prompt)',i));
    prompt = sprintf('Enter the b value in Gutenberg-Richter recurrence relation of area source %d',i);
    eval(sprintf('bA%d = input(prompt)',i));
end


prompt = 'Whether the same Ground Motion Prediction Equation to be used for all the sources? (1-yes, 0-no)';
check = input(prompt);
if check == 1
    prompt = 'Enter the serial number of the Ground Motion Prediction Equation to be used.';
    serial = input(prompt);
    if serial == 11
        prompt = 'Enter the value of Sa';
        Sa11 = input(prompt);
        prompt = 'Enter the value of Ss';
        Ss11 = input(prompt);
    end
else
    if Npoint>0
    for i = 1:Npoint
        prompt = sprintf('Enter the serial number of the Ground Motion Prediction Equation to be used for point source %d',i);
        serial_point(i) = input(prompt);
    end
    else
        serial_point = 0;
    end
    if Nline>0
    for i = 1:Nline
        prompt = sprintf('Enter the serial number of the Ground Motion Prediction Equation to be used for line source %d',i);
        serial_line(i) = input(prompt);
    end
    else
        serial_line = 0;
    end
    if Narea>0
    for i = 1:Narea
        prompt = sprintf('Enter the serial number of the Ground Motion Prediction Equation to be used for area source %d',i);
        serial_area(i) = input(prompt);
    end
    else
        serial_area = 0;
    end
end
prompt = 'Whether Hazard to be calculated in terms of spectral ordinates?[1-Yes,0-No] (NOTE: ALL SOURCES ARE CONSIDERED AT ONCE, INDIVIDUAL SOURCES NOT CONSIDERED SEPERATELY)';
specord_check = input(prompt);
if specord_check == 1
    prompt = 'Enter the appropriate serial number of the GMPE for predicting spectral ordinates to be used';
    serial_spectral = input(prompt);
else
    serial_spectral = 0;
end
if serial_spectral == 1001
    prompt = 'enter the SA value';
                    SA = input(prompt);
                    prompt = 'enter the SS value';
                    SS = input(prompt);
end
if check<1
serial_unified = horzcat(serial_point,serial_line,serial_area);
serial_unified(serial_unified == 0) = [];
end 
prompt = 'Enter the accuracy level required on a scale of 10';
Nintervals = input(prompt)*10;
prompt = 'For how many years should the poisson probability of atleast on occurance should be calculated';
Time_Period = input(prompt);
prompt = 'Enter the lower limit of the ground acceleration to be used';
Ylow = input(prompt);
prompt = 'Enter the upper limit of the ground acceleration to be used';
Yhigh = input(prompt);
prompt = 'Enter whether the site is soil or rock (0 for rock,1 for soil)';
S = input(prompt);
prompt = 'Whether Deaggregation should be done (1-YES,0-NO)';
Deagg_check = input(prompt);
if Deagg_check == 1
prompt = 'De-Aggregation to be done for how many values of ground acceleration';
NDeagg = input(prompt);
display('All values to be in terms of g only')
for i = 1:NDeagg
    prompt = sprintf('Enter the value of ground acceleration %d',i);
    temp = input(prompt);
    temp = round2(temp,0.01);
    Deagg(i) = temp;
%     temp = input(prompt);
%     eval(sprintf('GDeagg%d = temp',i));
end
deag = size(Deagg);
agg = deag(2);
end
prompt = 'Whether Uniform Hazard Spectrum should be obtained? [1-Yes,0-No] {NOTE: Without calculating the hazard in terms of Spectral Ordinates, this option cannot be selected as YES}';
Spectrum_check = input(prompt);
if Spectrum_check>0
    display('Hazard is defined as the Recurrence of exceeding a particular ground motion in given time period.');
    prompt = 'For how many values of Hazard, spectrum should be obtained?';
    Spectrum_number = input(prompt);
    for i = 1:Spectrum_number
        prompt = sprintf('Enter the value of Hazard %d',i);
        temp = input(prompt);
        Hazard_spectrum(i) = temp;
    end
end
 %%%%%%  END OF USER INPUTS %%%%%
%  clr
%  Nsources = 1;
%  Npoint = 0;
%  Nline = 1;
%  Narea = 0;
%  Xlineone1 = 5;
%  Ylineone1 = 15;
%  Xlineother1 = 7;
%  Ylineother1 = 18;
%  Mmaxl1 = 7.2;
%  Mminl1 = 4;
%  aL1 = 3;
%  bL1 = 0.8;
%  check = 1;
%  serial = 5;
%  Sa11 = 0;
%  Ss11 = 0;
%  specord_check = 0;
%  Nintervals = 10; %%%%%%%%%%
%  Time_Period = 50;
%  Ylow = 0.01;
%  Yhigh = 1;
%  Deagg_check = 0;
%  NDeagg = 1;
%  Spectrum_check = 0;
 %%%%%  DISTANCE UNCERTAINITY %%%%%
 % LINE SOURCE %
 if Nline>0
 for i = 1:Nline
     eval(sprintf('x1 = Xlineone%d',i));
     eval(sprintf('x2 = Xlineother%d',i));
     eval(sprintf('y1 = Ylineone%d',i));
     eval(sprintf('y2 = Ylineother%d',i));
     x = [x1 1;x2 1];
     y = [y1;y2];
     output = inv(x)*y;
     m = output(1);
     c = output(2);
     rr1 = sqrt(x1^2+y1^2);
     rr2 = sqrt(x2^2+y2^2);
     theta1 = asin(y1/rr1);
     theta2 = asin(y2/rr2);
     xx1 = -m*c/(1+m^2);
     yy1 = m*xx1+c;
         xmax = max(x1,x2);
         xmin = min(x1,x2);
         ymax = max(y1,y2);
         ymin = min(y1,y2);
         if inrange(xx1,[xmin,xmax],'includeboth') == 1
             xreq = xx1;
             yreq = yy1;
         
         elseif inrange(yy1,[ymin,ymax],'includeboth') == 1
             xreq = xx1;
             yreq = yy1;
         
         elseif rr1<rr2
             xreq = rr1*cos(theta1);
             yreq = rr1*sin(theta1);
         else
             xreq = rr2*cos(theta2);
             yreq = rr2*sin(theta2);
         end
     eval(sprintf('Rminline%d = sqrt(xreq^2+yreq^2)',i));
     d1 = sqrt(x1^2+y1^2);
     d2 = sqrt(x2^2+y2^2);
     if d1>d2
         eval(sprintf('Rmaxline%d = d1',i));
         Rimp = d2;
         dmax = d1;
         dmin = d2;
     else
         eval(sprintf('Rmaxline%d = d2',i));
         Rimp = d1;
         dmax = d2;
         dmin = d1;
     end
     Nsegmentsline = 1000;                                                   % CHANGE %
     Lsegmentline = sqrt((x1-x2)^2+(y1-y2)^2)/Nsegmentsline;
     eval(sprintf('k1 = Rmaxline%d',i));
     eval(sprintf('k2 = Rminline%d',i));
     dR = (k1-k2)/Nintervals;
     for j = 1:Nintervals
         eval(sprintf('r1 = Rminline%d+(j-1)*dR',i));
         eval(sprintf('r2 = Rminline%d+(j)*dR',i));
         eval(sprintf('Ravg_line%d(j) = (r1+r2)/2',i));
         xtest1a = (-2*m*c+sqrt(4*m^2*c^2-4*(1+m^2)*(c^2-r1^2)))/(2*(1+m^2));
         xt1a = round2(xtest1a,0.01);
         ytest1a = m*xtest1a+c;
         yt1a = round2(ytest1a,0.01);
         xtest1b = (-2*m*c-sqrt(4*m^2*c^2-4*(1+m^2)*(c^2-r1^2)))/(2*(1+m^2));
         xt1b = round2(xtest1b,0.01);
         ytest1b = m*xtest1b+c;
         yt1b = round2(ytest1b,0.01);
         xtest2a = (-2*m*c+sqrt(4*m^2*c^2-4*(1+m^2)*(c^2-r2^2)))/(2*(1+m^2));
         xt2a = round2(xtest2a,0.01);
         ytest2a = m*xtest2a+c;
         yt2a = round2(ytest2a,0.01);
         xtest2b = (-2*m*c-sqrt(4*m^2*c^2-4*(1+m^2)*(c^2-r2^2)))/(2*(1+m^2));
         xt2b = round2(xtest2b,0.01);
         ytest2b = m*xtest2b+c;
         yt2b = round2(ytest2b,0.01);
         z1 = inrange(xt1a,[xmin xmax],'includeboth');
         z2 = inrange(yt1a,[ymin ymax],'includeboth');
         z3 = inrange(xt2a,[xmin xmax],'includeboth');
         z4 = inrange(yt2a,[ymin ymax],'includeboth');
         if z1 == 1
             xtest1 = xtest1a;
         else
             xtest1 = xtest1b;
         end
         if z2 == 1
             ytest1 = ytest1a;
         else
             ytest1 = ytest1b;
         end
         if z3 == 1
             xtest2 = xtest2a;
         else
             xtest2 = xtest2b;
         end
         if z4 == 1
             ytest2 = ytest2a;
         else
             ytest2 = ytest2b;
         end
         eval(sprintf('l = Rminline%d',i));
         if abs(Rimp-l)<0.01
             dreq = sqrt((xtest1-xtest2)^2+(ytest1-ytest2)^2);
         else
             q = inrange(Rimp,[r1 r2],'includeright');
             if q == 1
                 dreq = sqrt(r2^2-l^2)-sqrt(r1^2-l^2)+sqrt(Rimp^2-l^2)-sqrt(r1^2-l^2);
             end
         if Rimp>r1
                 if Rimp>r2
                     dreq = 2*(sqrt(r2^2-l^2)-sqrt(r1^2-l^2));
                 end
         end
                 if Rimp<r1
                     if Rimp<r2
                   dreq = (sqrt(r2^2-l^2)-sqrt(r1^2-l^2));
                     end
                 end
         end
         segmentlength(j) = dreq;      
     end
     Nsegments = segmentlength/Lsegmentline;
     eval(sprintf('Dprob_line%d = Nsegments/Nsegmentsline',i));
     tt = Nsegments/Nsegmentsline';
     figure()
     eval(sprintf('t = Ravg_line%d',i));
     unified_distance2(i,:) = t;
     Deagg_distanceprob2(i,:) = tt;%%%
     eval(sprintf('prob = Dprob_line%d',i));
     bar(prob)
     set(gca,'XTickLabel',t,'XTick',1:numel(t))
     grid on
     grid minor
     title(sprintf('probability distribution in distance of line source %d',i))
     ylabel('probability')
     eval(sprintf('kk = sum(Dprob_line%d)',i));
     if abs(kk-1)>0.1
         i 
         display('CODE ERROR (IN DISTANCE PART OF LINE SOURCE)')
         beep
     else
         i
         display('SUCCESS')
     end
 end
 else
     unified_distance2 = zeros(1,Nintervals);
     Deagg_distanceprob2 = zeros(1,Nintervals);
 end
% Dprob_line1 = [0.25878029	0.136564378	0.110146686	0.089020236	0.06894073	0.068129061	0.06758634	0.06720315	0.06692143	0.066707698];
% Ravg_line1 = [7.39176727059238 16.1019793354353 24.8121914002781 33.5224034651210 42.2326155299639 50.9428275948068 59.6530396596497 68.3632517244925 77.0734637893354 85.7836758541783];
% unified_distance2 = Ravg_line1;

%  Ravg_line1 = Ravg_line1';
%  Dprob_line1 = Dprob_line1';
%  unified_distance2 = unified_distance2';
%  Deagg_distanceprob2 = importdata('C:\Users\SOMAYAJULU\Desktop\PROJECT\finished papers\Objective\results\PSHA_Distance.txt');
%  Deagg_distanceprob2 = Deagg_distanceprob2';
 ff = 0;
 if Nline>0
 for i = 1:Nline
     unified_dummy2(1,1+ff*Nintervals:Nintervals+Nintervals*ff) = unified_distance2(i,:);
     distprob2(1,1+ff*Nintervals:Nintervals+Nintervals*ff) = Deagg_distanceprob2(i,:);
     ff = ff+1;
 end
 unified_distance2 = unified_dummy2;
 Deagg_distprob2 = distprob2;
 else
     Deagg_distprob2 = 0;
 end
    dist2 = unified_distance2;
 % AREA SOURCE %
 % unified_distance3
 % CHANGES TO BE MADE IN THE UNIFIED_DISTANCE PART %
 % CHANGES TO BE MADE IN THE DEAGG_DIST PART AS ONLY UNIFIED DISTANCE 2 IS ONLY CONSIDERED. FOR AREA SOURCE ALSO INCLUDED, UNIFIED DISTANCE2 AND 3 HAVE TO BE CONCATINATED %
 % CHANGES TO BE MADE TO THE DEAGG DIST PROB PART %
 % POINT SOURCE %
 if Npoint>0
 for i = 1:Npoint
     prob = ones(1,Nintervals);
     eval(sprintf('Dprob_point%d = prob',i));
     eval(sprintf('t1 = Xpoint%d',i));
     eval(sprintf('t2 = Ypoint%d',i));
     t = sqrt(t1^2+t2^2);
     eval(sprintf('dist_prob%d = t',i));
     unified_distance1(i,:) = t*prob';
     Deagg_distprob1(i) = 1;
 end
 else
     unified_distance1 = zeros(1,Nintervals);
     Deagg_distprob1 = 0;
 end
 ff = 0;
 if Npoint>0
 for i = 1:Npoint
     unified_dummy1(1,1+ff*Nintervals:Nintervals+Nintervals*ff) = unified_distance1(i,:);
     ff = ff+1;
 end
 unified_distance1 = unified_dummy1;
 end
 unified_distance = horzcat(unified_distance1,unified_distance2);
 unified_distance(unified_distance == 0) = [];
 size1 = size(unified_distance);
 Deagg_distprob = horzcat(Deagg_distprob1,Deagg_distprob2);
 Deagg_distprob(Deagg_distprob == 0) = [];
 
 %%%%%% UNCERTAINITY IN MAGNITUDE %%%%%%%
 % POINT SOURCE %
 if Npoint>0
 for i = 1:Npoint
     eval(sprintf('Mmax = Mmaxp%d',i));
     eval(sprintf('Mmin = Mminp%d',i));
     eval(sprintf('a = aP%d',i));
     eval(sprintf('b = bP%d',i));
     dM = (Mmax-Mmin)/Nintervals;
     for j = 1:Nintervals
         M1 = Mmin+(j-1)*dM;
         M2 = Mmin+j*dM;
         Mavg = (M1+M2)/2;
         m(j) = Mavg;
         fM = 2.303*b*exp(-2.303*b*(Mavg-Mmin))/(1-exp(-2.303*b*(Mmax-Mmin)));
         prob(j) = fM*(M2-M1);
     end
     unified_magnitude1(i,:) = m;
     Deagg_magnprob1(i,:) = prob;
     figure()
     bar(prob,'green')
     set(gca,'XTickLabel',m,'XTick',1:numel(m))
     grid on
     grid minor
     title(sprintf('probability distribution in magnitude of point source %d',i))
     ylabel('probability')
     kk = sum(prob);
     eval(sprintf('Mprob_point%d = prob',i));
     eval(sprintf('Mavg_point%d = m',i));
     if abs(kk-1)>0.1
         i 
         display('CODE ERROR (IN MAGNITUDE PART OF POINT SOURCE)')
         beep
     else
         i
         display('SUCCESS')
     end
 end
 else
     unified_magnitude1 = zeros(1,Nintervals);
     Deagg_magnprob1 = zeros(1,Nintervals);
 end
 ff = 0;
 if Npoint>0
 for i = 1:Npoint
     unified_dummy1(1,1+ff*Nintervals:Nintervals+Nintervals*ff) = unified_magnitude1(i,:);
     Deagg_dummy1(1,1+ff*Nintervals:Nintervals+Nintervals*ff) = Deagg_magnprob1(i,:);
     ff = ff+1;
 end
 unified_magnitude1 = unified_dummy1;
 Deagg_magnprob1 = Deagg_dummy1;
 end
 
 % LINE SOURCE %
 if Nline>0
 for i = 1:Nline
     eval(sprintf('Mmax = Mmaxl%d',i));
     eval(sprintf('Mmin = Mminl%d',i));
     eval(sprintf('a = aL%d',i));
     eval(sprintf('b = bL%d',i));
     dM = (Mmax-Mmin)/Nintervals;
     for j = 1:Nintervals
         M1 = Mmin+(j-1)*dM;
         M2 = Mmin+j*dM;
         Mavg = (M1+M2)/2;
         m(j) = Mavg;
         fM = 2.303*b*exp(-2.303*b*(Mavg-Mmin))/(1-exp(-2.303*b*(Mmax-Mmin)));
         prob(j) = fM*(M2-M1);
     end
     unified_magnitude2(i,:) = m;
     Deagg_magnprob2(i,:) = prob;
     figure()
     bar(prob,'green')
     set(gca,'XTickLabel',m,'XTick',1:numel(m))
     grid on
     grid minor
     title(sprintf('probability distribution in magnitude of line source %d',i))
     ylabel('probability')
     kk = sum(prob);
     eval(sprintf('Mprob_line%d = prob',i));
     eval(sprintf('Mavg_line%d = m',i));
     if abs(kk-1)>0.1
         i 
         display('CODE ERROR (IN MAGNITUDE PART OF LINE SOURCE)')
         beep
     else
         i
         display('SUCCESS')
     end
 end
 else
     unified_magnitude2 = zeros(1,Nintervals);
     Deagg_magnprob2 = zeros(1,Nintervals);
 end
 ff = 0;
 if Nline>0
 for i = 1:Nline
     unified_dummy2(1,1+ff*Nintervals:Nintervals+Nintervals*ff) = unified_magnitude2(i,:);
     Deagg_dummy2(1,1+ff*Nintervals:Nintervals+Nintervals*ff) = Deagg_magnprob2(i,:);
     ff = ff+1;
 end
 unified_magnitude2 = unified_dummy2;
 Deagg_magnprob2 = Deagg_dummy2;
 end
 % AREA SOURCE %
 if Narea>0
 for i = 1:Narea
     eval(sprintf('Mmax = Mmaxa%d',i));
     eval(sprintf('Mmin = Mmina%d',i));
     eval(sprintf('a = aA%d',i));
     eval(sprintf('b = bA%d',i));
     dM = (Mmax-Mmin)/Nintervals;
     for j = 1:Nintervals
         M1 = Mmin+(j-1)*dM;
         M2 = Mmin+j*dM;
         Mavg = (M1+M2)/2;
         m(j) = Mavg;
         fM = 2.303*b*exp(-2.303*b*(Mavg-Mmin))/(1-exp(-2.303*b*(Mmax-Mmin)));
         prob(j) = fM*(M2-M1);
     end
     unified_magnitude3(i,:) = m;
     Deagg_magnprob3(i,:) = prob;
     figure()
     bar(prob,'green')
     set(gca,'XTickLabel',m,'XTick',1:numel(m))
     grid on
     grid minor
     title(sprintf('probability distribution in magnitude of area source %d',i))
     ylabel('probability')
     kk = sum(prob);
     eval(sprintf('Mprob_area%d = prob',i));
     eval(sprintf('Mavg_area%d = m',i));
     if abs(kk-1)>0.1
         i 
         display('CODE ERROR (IN MAGNITUDE PART OF LINE SOURCE)')
         beep
     else
         i
         display('SUCCESS')
     end
 end
 else
     unified_magnitude3 = zeros(1,Nintervals);
     Deagg_magnprob3 = zeros(1,Nintervals);
 end
 ff = 0;
 if Narea>0
 for i = 1:Narea
     unified_dummy3(1,1+ff*Nintervals:Nintervals+Nintervals*ff) = unified_magnitude3(i,:);
     Deagg_dummy3(1,1+ff*Nintervals:Nintervals+Nintervals*ff) = Deagg_magnprob3(i,:);
     ff = ff+1;
 end
 unified_magnitude3 = unified_dummy3;
 Deagg_magnprob3 = Deagg_dummy3;
 end
 unified_magnitude = horzcat(unified_magnitude1,unified_magnitude2,unified_magnitude3);
 unified_magnitude(unified_magnitude == 0) = [];
 Deagg_magnprob = horzcat(Deagg_magnprob1,Deagg_magnprob2,Deagg_magnprob3);
 Deagg_magnprob(Deagg_magnprob == 0) = [];
 magn = horzcat(unified_magnitude1,unified_magnitude2,unified_magnitude3);
 magn(magn == 0) = [];
 size2 = size(unified_magnitude);
 if check<1
 size3 = size(serial_unified);
 end
 if Npoint>0
 for i = 1:Npoint
     dist1(i) = unified_distance(i*Nintervals);
 end
 else
     dist1 = 0;
 end
 Deagg_distance = horzcat(dist1,dist2);
 Deagg_distance(Deagg_distance == 0) = [];
 Deagg_magnitude = magn;
 
 %%%%% UNCERTAINITY IN GROUND MOTION PREDICTION EQUATION %%%%%
%  GROUND MOTION PREDICTION EQUATIONS %
f = 1;
dummy = 1;
Nacc = (Yhigh-Ylow)/0.01+1;
for Y = Ylow:0.01:Yhigh                   
if (check+specord_check) > 0
 for i = 1:size1(2)
     for j = 1:size2(2)
         if check == 1
         if serial == 5
             PGA = (6.74+0.859*unified_magnitude(j)-1.8*log(unified_distance(i)+25));
             sigma = 0.57;
             prob_PGA(j,i) = 1-normcdf(log(Y*981),PGA,sigma);
         end
                 if serial == 2
                    PGA = 0.46+0.63*unified_magnitude(j)-1.1*log10(unified_distance(i));
                    sigma = 0.32;
                    prob_PGA(j,i) = 1-normcdf(log10(Y*981),PGA,sigma);
                 end
                if serial == 4
                    PGA = 3.4+0.89*unified_magnitude(j)-1.17*log(unified_distance(i))-0.2*S;
                    sigma = 0.62;
                    prob_PGA(j,i) = 1-normcdf(log(Y*981),PGA,sigma);
                end
                if serial == 3
                    PGA = 0.0159*exp(0.868*unified_magnitude(j))*(unified_distance(i)+0.0606*exp(0.7*unified_magnitude(j)))^(-1.09);
                    sigma = 0.372;
                    prob_PGA(j,i) = 1-normcdf(log(Y),log(PGA),sigma);
                end
                if serial == 6
                    PGA = 0.0185*exp(1.28*unified_magnitude(j))*(unified_distance(i)+0.147*exp(0.732*unified_magnitude(j)))^(-1.75);
                    sigma = 0.384;
                    prob_PGA(j,i) = 1-normcdf(log(Y),log(PGA),sigma);
                end
                if serial == 7
                    R = sqrt(unified_distance(i)^2+7.3^2);
                    PGA = (-1.02+0.249*unified_magnitude(j)-log10(R)-0.00255*R);
                    sigma = 0.26;
                    prob_PGA(j,i) = 1-normcdf(log10(Y),(PGA),sigma);
                end
                if serial == 8
                    R = sqrt(unified_distance(i)^2+8^2);
                    PGA = 0.43+0.23*(unified_magnitude(j)-6)-log10(R)-0.0027*R;
                    sigma = 0.28;
                    prob_PGA(j,i) = 1-normcdf(log10(Y),(PGA),sigma);
                end
                if serial == 9
                    PGA = 0.41*unified_magnitude(j)-log10(unified_distance(i)+0.032*10^(0.41*unified_magnitude(j)))-0.0034*unified_distance(i)+1.3;
                    sigma = 0.21;
                    prob_PGA(j,i) = 1-normcdf(log10(Y*981),PGA,sigma);
                end
                if serial == 10
                    PGA = -2.501+0.623*unified_magnitude(j)-log(unified_distance(i)+7.28);
                    sigma = 0.506;
                    prob_PGA(j,i) = 1-normcdf(log(Y),PGA,sigma);
                end
                if serial == 11
                    PGA = -0.84+0.219*unified_magnitude(j)-0.954*log10(sqrt(unified_distance(i).^2+4.5^2))+0.117*Sa11+0.124*Ss11;
                    sigma = 0.27;
                    prob_PGA(j,i) = 1-normcdf(log10(Y),PGA,sigma);
                end
                if serial == 1
                    PGA = -3.17+0.508*unified_magnitude(j)-0.885*log10(sqrt(unified_distance(i).^2+4.3^2))+0.117*0+0.124*0;
                    sigma = 0.32;
                    prob_PGA(j,i) = 1-normcdf(log10(Y),PGA,sigma);
                end
                if serial == 15
                    if BA08_check == 1
                    fprintf(   'S.no    Ground motion intensity measure\n')
                    names = {'1'  'PGV';'2'  'PGA';'3'  'Sa'}; 
                    disp(names)
                    prompt = 'Enter the appropriate serial number';
                    sno = input(prompt);
                    if sno == 3
                        prompt = 'Enter the time period for Sa between 0.01 - 10 seconds';
                        Tp = input(prompt);
                    end
                    table3 = [-0.600000000000000,-0.360000000000000,-0.360000000000000,-0.340000000000000,-0.330000000000000,-0.290000000000000,-0.230000000000000,-0.250000000000000,-0.280000000000000,-0.310000000000000,-0.390000000000000,-0.440000000000000,-0.500000000000000,-0.600000000000000,-0.690000000000000,-0.700000000000000,-0.720000000000000,-0.730000000000000,-0.740000000000000,-0.750000000000000,-0.750000000000000,-0.692000000000000,-0.650000000000000;...
                        -0.500000000000000,-0.640000000000000,-0.640000000000000,-0.630000000000000,-0.620000000000000,-0.640000000000000,-0.640000000000000,-0.600000000000000,-0.530000000000000,-0.520000000000000,-0.520000000000000,-0.520000000000000,-0.510000000000000,-0.500000000000000,-0.470000000000000,-0.440000000000000,-0.400000000000000,-0.380000000000000,-0.340000000000000,-0.310000000000000,-0.291000000000000,-0.247000000000000,-0.215000000000000;...
                        -0.0600000000000000,-0.140000000000000,-0.140000000000000,-0.120000000000000,-0.110000000000000,-0.110000000000000,-0.110000000000000,-0.130000000000000,-0.180000000000000,-0.190000000000000,-0.160000000000000,-0.140000000000000,-0.100000000000000,-0.0600000000000000,0,0,0,0,0,0,0,0,0];
                    table6 = [-0.873700000000000,-0.660500000000000,-0.662200000000000,-0.666000000000000,-0.690100000000000,-0.717000000000000,-0.720500000000000,-0.708100000000000,-0.696100000000000,-0.583000000000000,-0.572600000000000,-0.554300000000000,-0.644300000000000,-0.691400000000000,-0.740800000000000,-0.818300000000000,-0.830300000000000,-0.828500000000000,-0.784400000000000,-0.685400000000000,-0.509600000000000,-0.372400000000000,-0.0982400000000000;...
                        0.100600000000000,0.119700000000000,0.120000000000000,0.122800000000000,0.128300000000000,0.131700000000000,0.123700000000000,0.111700000000000,0.0988400000000000,0.0427300000000000,0.0297700000000000,0.0195500000000000,0.0439400000000000,0.0608000000000000,0.0751800000000000,0.102700000000000,0.0979300000000000,0.0943200000000000,0.0728200000000000,0.0375800000000000,-0.0239100000000000,-0.0656800000000000,-0.138000000000000;...
                        -0.00334000000000000,-0.0115100000000000,-0.0115100000000000,-0.0115100000000000,-0.0115100000000000,-0.0115100000000000,-0.0115100000000000,-0.0115100000000000,-0.0111300000000000,-0.00952000000000000,-0.00837000000000000,-0.00750000000000000,-0.00626000000000000,-0.00540000000000000,-0.00409000000000000,-0.00334000000000000,-0.00255000000000000,-0.00217000000000000,-0.00191000000000000,-0.00191000000000000,-0.00191000000000000,-0.00191000000000000,-0.00191000000000000;...
                        2.54000000000000,1.35000000000000,1.35000000000000,1.35000000000000,1.35000000000000,1.35000000000000,1.55000000000000,1.68000000000000,1.86000000000000,1.98000000000000,2.07000000000000,2.14000000000000,2.24000000000000,2.32000000000000,2.46000000000000,2.54000000000000,2.66000000000000,2.73000000000000,2.83000000000000,2.89000000000000,2.93000000000000,3,3.04000000000000];
                    table7 = [5.00121000000000,5.04727000000000,4.63188000000000,5.08210000000000,0.183220000000000,-0.127360000000000,0,8.50000000000000;...
                        -0.538040000000000,-0.503500000000000,-0.754720000000000,-0.509700000000000,0.288050000000000,-0.101640000000000,0,6.75000000000000;...
                        -0.528830000000000,-0.494290000000000,-0.745510000000000,-0.499660000000000,0.288970000000000,-0.100190000000000,0,6.75000000000000;...
                        -0.521920000000000,-0.485080000000000,-0.739060000000000,-0.488950000000000,0.251440000000000,-0.110060000000000,0,6.75000000000000;...
                        -0.452850000000000,-0.418310000000000,-0.667220000000000,-0.422290000000000,0.179760000000000,-0.128580000000000,0,6.75000000000000;...
                        -0.284760000000000,-0.250220000000000,-0.484620000000000,-0.260920000000000,0.0636900000000000,-0.157520000000000,0,6.75000000000000;...
                        0.00767000000000000,0.0491200000000000,-0.205780000000000,0.0270600000000000,0.0117000000000000,-0.170510000000000,0,6.75000000000000;...
                        0.201090000000000,0.231020000000000,0.0305800000000000,0.221930000000000,0.0469700000000000,-0.159480000000000,0,6.75000000000000;...
                        0.461280000000000,0.486610000000000,0.301850000000000,0.493280000000000,0.179900000000000,-0.145390000000000,0,6.75000000000000;...
                        0.571800000000000,0.592530000000000,0.408600000000000,0.614720000000000,0.527290000000000,-0.129640000000000,0.00102000000000000,6.75000000000000;...
                        0.518840000000000,0.534960000000000,0.338800000000000,0.577470000000000,0.608800000000000,-0.138430000000000,0.0860700000000000,6.75000000000000;...
                        0.438250000000000,0.445160000000000,0.253560000000000,0.519900000000000,0.644720000000000,-0.156940000000000,0.106010000000000,6.75000000000000;...
                        0.392200000000000,0.406020000000000,0.213980000000000,0.460800000000000,0.786100000000000,-0.0784300000000000,0.0226200000000000,6.75000000000000;...
                        0.189570000000000,0.198780000000000,0.00967000000000000,0.263370000000000,0.768370000000000,-0.0905400000000000,0,6.75000000000000;...
                        -0.213380000000000,-0.194960000000000,-0.491760000000000,-0.108130000000000,0.751790000000000,-0.140530000000000,0.103020000000000,6.75000000000000;...
                        -0.468960000000000,-0.434430000000000,-0.784650000000000,-0.393300000000000,0.678800000000000,-0.182570000000000,0.0539300000000000,6.75000000000000;...
                        -0.862710000000000,-0.795930000000000,-1.20902000000000,-0.880850000000000,0.706890000000000,-0.259500000000000,0.190820000000000,6.75000000000000;...
                        -1.22652000000000,-1.15514000000000,-1.57697000000000,-1.27669000000000,0.779890000000000,-0.296570000000000,0.298880000000000,6.75000000000000;...
                        -1.82979000000000,-1.74690000000000,-2.22584000000000,-1.91814000000000,0.779660000000000,-0.453840000000000,0.674660000000000,6.75000000000000;...
                        -2.24656000000000,-2.15906000000000,-2.58228000000000,-2.38168000000000,1.24961000000000,-0.358740000000000,0.795080000000000,6.75000000000000;...
                        -1.28408000000000,-1.21270000000000,-1.50904000000000,-1.41093000000000,0.142710000000000,-0.390060000000000,0,8.50000000000000;...
                        -1.43145000000000,-1.31632000000000,-1.81022000000000,-1.59217000000000,0.524070000000000,-0.375780000000000,0,8.50000000000000;...
                        -2.15446000000000,-2.16137000000000,-2.53323000000000,-2.14635000000000,0.403870000000000,-0.484920000000000,0,8.50000000000000];
                    table8 = [0.560000000000000,0.564000000000000,0.566000000000000,0.566000000000000,0.576000000000000,0.589000000000000,0.606000000000000,0.608000000000000,0.594000000000000,0.596000000000000,0.592000000000000,0.608000000000000,0.603000000000000,0.615000000000000,0.645000000000000,0.647000000000000,0.679000000000000,0.700000000000000,0.695000000000000,0.698000000000000,0.744000000000000,0.787000000000000,0.801000000000000];
                        if sno == 1
                            blin = table3(1,1);
                            b1 = table3(2,1);
                            b2 = table3(3,1);
                            c1 = table6(1,1);
                            c2 = table6(2,1);
                            c3 = table6(3,1);
                            h1 = table6(4,1);
                            e1 = table7(1,1);
                            e2 = table7(1,2);
                            e3 = table7(1,3);
                            e4 = table7(1,4);
                            e5 = table7(1,5);
                            e6 = table7(1,6);
                            e7 = table7(1,7);
                            Mh = table7(1,8);
                            sigma = table8(1);
                        elseif sno == 2
                            blin = table3(1,2);
                            b1 = table3(2,2);
                            b2 = table3(3,2);
                            c1 = table6(1,2);
                            c2 = table6(2,2);
                            c3 = table6(3,2);
                            h1 = table6(4,2);
                            e1 = table7(2,1);
                            e2 = table7(2,2);
                            e3 = table7(2,3);
                            e4 = table7(2,4);
                            e5 = table7(2,5);
                            e6 = table7(2,6);
                            e7 = table7(2,7);
                            Mh = table7(2,8);
                            sigma = table8(2);
                        elseif sno == 3
                            Tfunda = [0.01 0.02 0.03 0.05 0.075 0.1 0.15 0.2 0.25 0.3 0.4 0.5 0.75 1 1.5 2 3 4 5 7.5 10];
                            lower_index = (max(find(Tfunda <= Tp)));
                            upper_index = (min(find(Tfunda >= Tp)));
                            lower_Tp = Tfunda(lower_index);
                            upper_Tp = Tfunda(upper_index);
                            blin_lower = table3(1,lower_index+2);
                            b1_lower = table3(2,lower_index+2);
                            b2_lower = table3(3,lower_index+2);
                            c1_lower = table6(1,lower_index+2);
                            c2_lower = table6(2,lower_index+2);
                            c3_lower = table6(3,lower_index+2);
                            h_lower = table6(4,lower_index+2);
                            e1_lower = table7(lower_index+2,1);
                            e2_lower = table7(lower_index+2,2);
                            e3_lower = table7(lower_index+2,3);
                            e4_lower = table7(lower_index+2,4);
                            e5_lower = table7(lower_index+2,5);
                            e6_lower = table7(lower_index+2,6);
                            e7_lower = table7(lower_index+2,7);
                            Mh_lower = table7(lower_index+2,8);
                            sigma_lower = table8(lower_index+2);
                            
                            blin_upper = table3(1,upper_index+2);
                            b1_upper = table3(2,upper_index+2);
                            b2_upper = table3(3,upper_index+2);
                            c1_upper = table6(1,upper_index+2);
                            c2_upper = table6(2,upper_index+2);
                            c3_upper = table6(3,upper_index+2);
                            h_upper = table6(4,upper_index+2);
                            e1_upper = table7(upper_index+2,1);
                            e2_upper = table7(upper_index+2,2);
                            e3_upper = table7(upper_index+2,3);
                            e4_upper = table7(upper_index+2,4);
                            e5_upper = table7(upper_index+2,5);
                            e6_upper = table7(upper_index+2,6);
                            e7_upper = table7(upper_index+2,7);
                            Mh_upper = table7(upper_index+2,8);
                            sigma_upper = table8(upper_index+2);
                            
                            blin = blin_lower + (blin_upper-blin_lower)/(upper_Tp-lower_Tp) * (Tp - lower_Tp);
                            b1 = b1_lower + (b1_upper-b1_lower)/(upper_Tp-lower_Tp) * (Tp - lower_Tp);
                            b2 = b2_lower + (b2_upper-b2_lower)/(upper_Tp-lower_Tp) * (Tp - lower_Tp);
                            c1 = c1_lower + (c1_upper-c1_lower)/(upper_Tp-lower_Tp) * (Tp - lower_Tp);
                            c2 = c2_lower + (c2_upper-c2_lower)/(upper_Tp-lower_Tp) * (Tp - lower_Tp);
                            c3 = c3_lower + (c3_upper-c3_lower)/(upper_Tp-lower_Tp) * (Tp - lower_Tp);
                            h1 = h_lower + (h_upper-h_lower)/(upper_Tp-lower_Tp) * (Tp - lower_Tp);
                            e1 = e1_lower + (e1_upper-e1_lower)/(upper_Tp-lower_Tp) * (Tp - lower_Tp);
                            e2 = e2_lower + (e2_upper-e2_lower)/(upper_Tp-lower_Tp) * (Tp - lower_Tp);
                            e3 = e3_lower + (e3_upper-e3_lower)/(upper_Tp-lower_Tp) * (Tp - lower_Tp);
                            e4 = e4_lower + (e4_upper-e4_lower)/(upper_Tp-lower_Tp) * (Tp - lower_Tp);
                            e5 = e5_lower + (e5_upper-e5_lower)/(upper_Tp-lower_Tp) * (Tp - lower_Tp);
                            e6 = e6_lower + (e6_upper-e6_lower)/(upper_Tp-lower_Tp) * (Tp - lower_Tp);
                            e7 = e7_lower + (e7_upper-e7_lower)/(upper_Tp-lower_Tp) * (Tp - lower_Tp);
                            Mh = Mh_lower + (Mh_upper-Mh_lower)/(upper_Tp-lower_Tp) * (Tp - lower_Tp);
                            sigma = sigma_lower + (sigma_upper-sigma_lower)/(upper_Tp-lower_Tp) * (Tp - lower_Tp);
                            if isnan(blin) == 1
                                blin = blin_lower;
                            b1 = b1_lower;
                            b2 = b2_lower;
                            c1 = c1_lower;
                            c2 = c2_lower;
                            c3 = c3_lower;
                            h1 = h_lower;
                            e1 = e1_lower;
                            e2 = e2_lower;
                            e3 = e3_lower;
                            e4 = e4_lower;
                            e5 = e5_lower;
                            e6 = e6_lower;
                            e7 = e7_lower;
                            Mh = Mh_lower;
                            sigma = sigma_lower;
                            end
                        end
                    Mref = 4.5;
                    Rref = 1;
                    a1 = 0.03;
                    pga_low = 0.06;
                    a2 = 0.09;
                    v1 = 180;
                    v2 = 300;
                    vref = 760;  
                    prompt = 'What is the type of fault? (SS 0, RS 1)';
                    fault_type = input(prompt);
                    if fault_type == 0
                        SS = 1;
                        RS = 0;
                    else
                        RS = 1;
                        SS = 0;
                    end
prompt = 'What is the shear wave velocity in m/s ?';
                    Vs = input(prompt);
                    BA08_check = BA08_check + 1;
                    end
%ii = 1;
R = sqrt(unified_distance(i)^2+h1^2);
Rp = sqrt(unified_distance(i)^2+1.35^2);
FD = (c1+c2*(unified_magnitude(j)-Mref))*log(R/Rref)+c3*(R-Rref);
FDp = (-0.66050+0.11970*(unified_magnitude(j)-Mref))*log(Rp/5)-0.01151*(Rp-5);
if unified_magnitude(j)<=Mh
FM = e2*SS+e4*RS+e5*(unified_magnitude(j)-Mh)+e6*(unified_magnitude(j)-Mh)^2;
else
    FM = e2*SS+e4*RS+e7*(unified_magnitude(j)-Mh);
end
if unified_magnitude(j)<=6.75
FMp = -0.50350*SS-0.50970*RS+0.28805*(unified_magnitude(j)-6.75)-0.10164*(unified_magnitude(j)-6.75).^2;
else
    FMp = -0.50350*SS-0.50970*RS;
end
pga4nl = exp(FMp+FDp);
Flin = blin*log(Vs/vref);
if Vs<=v1
bnl = b1;
elseif Vs>v1 && Vs<=v2
bnl = (b1-b2)*log(Vs/v2)/log(v1/v2)+b2;
elseif Vs>v2 && Vs<vref
    bnl = b2*log(Vs/vref)/log(v2/vref);
else
    bnl = 0;
end
dx = log(a2/a1);
dy = bnl*log(a2/pga_low);
C = (3*dy-bnl*dx)/dx^2;
D = -(2*dy-bnl*dx)/dx^3;
if pga4nl <= a1
    Fnl = bnl*log(pga_low/0.1);
elseif pga4nl>a1 && pga4nl<=a2
    Fnl = bnl*log(pga_low/0.1)+C*(log(pga4nl/a1))^2+D*(log(pga4nl/a1))^3;
else
    Fnl = bnl*log(pga4nl/0.1);
end
PGA = (FM+FD+Flin+Fnl);
prob_PGA(j,i) = 1-normcdf(log(Y),PGA,sigma);
                end
         end
         if specord_check == 1
                if serial_spectral == 1001
                    period_Number = 46; %% VERY IMPORTANT %%
                    T = [0.100000000000000,0.110000000000000,0.120000000000000,0.130000000000000,0.140000000000000,0.150000000000000,0.160000000000000,0.170000000000000,0.180000000000000,0.190000000000000,0.200000000000000,0.220000000000000,0.240000000000000,0.260000000000000,0.280000000000000,0.300000000000000,0.320000000000000,0.340000000000000,0.360000000000000,0.380000000000000,0.400000000000000,0.420000000000000,0.440000000000000,0.460000000000000,0.480000000000000,0.500000000000000,0.550000000000000,0.600000000000000,0.650000000000000,0.700000000000000,0.750000000000000,0.800000000000000,0.850000000000000,0.900000000000000,0.950000000000000,1,1.10000000000000,1.20000000000000,1.30000000000000,1.40000000000000,1.50000000000000,1.60000000000000,1.70000000000000,1.80000000000000,1.90000000000000,2];
                    c1 = [-0.840000000000000,-0.860000000000000,-0.870000000000000,-0.870000000000000,-0.940000000000000,-0.980000000000000,-1.05000000000000,-1.08000000000000,-1.13000000000000,-1.19000000000000,-1.21000000000000,-1.28000000000000,-1.37000000000000,-1.40000000000000,-1.46000000000000,-1.55000000000000,-1.63000000000000,-1.65000000000000,-1.69000000000000,-1.82000000000000,-1.94000000000000,-1.99000000000000,-2.05000000000000,-2.11000000000000,-2.17000000000000,-2.25000000000000,-2.38000000000000,-2.49000000000000,-2.58000000000000,-2.67000000000000,-2.75000000000000,-2.86000000000000,-2.93000000000000,-3.03000000000000,-3.10000000000000,-3.17000000000000,-3.30000000000000,-3.38000000000000,-3.43000000000000,-3.52000000000000,-3.61000000000000,-3.68000000000000,-3.74000000000000,-3.79000000000000,-3.80000000000000,-3.79000000000000];
                    c2 = [0.210000000000000,0.221000000000000,0.231000000000000,0.238000000000000,0.244000000000000,0.247000000000000,0.252000000000000,0.258000000000000,0.268000000000000,0.278000000000000,0.284000000000000,0.295000000000000,0.308000000000000,0.318000000000000,0.326000000000000,0.338000000000000,0.349000000000000,0.351000000000000,0.354000000000000,0.364000000000000,0.377000000000000,0.384000000000000,0.393000000000000,0.401000000000000,0.410000000000000,0.420000000000000,0.434000000000000,0.438000000000000,0.451000000000000,0.463000000000000,0.477000000000000,0.485000000000000,0.492000000000000,0.502000000000000,0.503000000000000,0.508000000000000,0.513000000000000,0.513000000000000,0.514000000000000,0.522000000000000,0.524000000000000,0.520000000000000,0.517000000000000,0.514000000000000,0.508000000000000,0.503000000000000];
                    ho = [4.50000000000000,4.50000000000000,4.70000000000000,5.30000000000000,4.90000000000000,4.70000000000000,4.40000000000000,4.30000000000000,4,3.90000000000000,4.20000000000000,4.10000000000000,3.90000000000000,4.30000000000000,4.40000000000000,4.2,4.20000000000000,4.40000000000000,4.5,3.90000000000000,3.60000000000000,3.70000000000000,3.90000000000000,3.70000000000000,3.50000000000000,3.30000000000000,3.10000000000000,2.50000000000000,2.80000000000000,3.10000000000000,3.50000000000000,3.70000000000000,3.90000000000000,4,4,4.30000000000000,4,3.60000000000000,3.60000000000000,3.40000000000000,3,2.50000000000000,2.50000000000000,2.40000000000000,2.80000000000000,3.20000000000000];
                    c4 = [-0.954000000000000,-0.945000000000000,-0.960000000000000,-0.981000000000000,-0.955000000000000,-0.938000000000000,-0.907000000000000,-0.896000000000000,-0.901000000000000,-0.907000000000000,-0.922000000000000,-0.911000000000000,-0.916000000000000,-0.942000000000000,-0.946000000000000,-0.933000000000000,-0.932000000000000,-0.939000000000000,-0.936000000000000,-0.900000000000000,-0.888000000000000,-0.897000000000000,-0.908000000000000,-0.911000000000000,-0.920000000000000,-0.913000000000000,-0.911000000000000,-0.881000000000000,-0.901000000000000,-0.914000000000000,-0.942000000000000,-0.925000000000000,-0.920000000000000,-0.920000000000000,-0.892000000000000,-0.885000000000000,-0.857000000000000,-0.851000000000000,-0.848000000000000,-0.839000000000000,-0.817000000000000,-0.781000000000000,-0.759000000000000,-0.730000000000000,-0.724000000000000,-0.728000000000000];
                    ca = [0.0780000000000000,0.0980000000000000,0.111000000000000,0.131000000000000,0.136000000000000,0.143000000000000,0.152000000000000,0.140000000000000,0.129000000000000,0.133000000000000,0.135000000000000,0.120000000000000,0.124000000000000,0.134000000000000,0.134000000000000,0.133000000000000,0.125000000000000,0.118000000000000,0.124000000000000,0.132000000000000,0.139000000000000,0.147000000000000,0.153000000000000,0.149000000000000,0.100000000000000,0.147000000000000,0.134000000000000,0.124000000000000,0.122000000000000,0.116000000000000,0.113000000000000,0.127000000000000,0.124000000000000,0.124000000000000,0.121000000000000,0.128000000000000,0.123000000000000,0.128000000000000,0.115000000000000,0.109000000000000,0.109000000000000,0.108000000000000,0.105000000000000,0.104000000000000,0.103000000000000,0.101000000000000];
                    cs = [0.0270000000000000,0.0360000000000000,0.0520000000000000,0.0680000000000000,0.0770000000000000,0.0850000000000000,0.101000000000000,0.102000000000000,0.107000000000000,0.130000000000000,0.142000000000000,0.143000000000000,0.155000000000000,0.163000000000000,0.158000000000000,0.148000000000000,0.161000000000000,0.163000000000000,0.160000000000000,0.164000000000000,0.172000000000000,0.180000000000000,0.187000000000000,0.191000000000000,0.197000000000000,0.201000000000000,0.203000000000000,0.212000000000000,0.215000000000000,0.214000000000000,0.212000000000000,0.218000000000000,0.218000000000000,0.225000000000000,0.217000000000000,0.219000000000000,0.206000000000000,0.214000000000000,0.200000000000000,0.197000000000000,0.204000000000000,0.206000000000000,0.206000000000000,0.204000000000000,0.194000000000000,0.182000000000000];
                    sig = [0.270000000000000,0.270000000000000,0.270000000000000,0.270000000000000,0.270000000000000,0.270000000000000,0.270000000000000,0.270000000000000,0.270000000000000,0.280000000000000,0.270000000000000,0.280000000000000,0.280000000000000,0.280000000000000,0.290000000000000,0.300000000000000,0.310000000000000,0.310000000000000,0.310000000000000,0.310000000000000,0.310000000000000,0.320000000000000,0.320000000000000,0.320000000000000,0.320000000000000,0.320000000000000,0.320000000000000,0.320000000000000,0.320000000000000,0.330000000000000,0.320000000000000,0.320000000000000,0.320000000000000,0.320000000000000,0.320000000000000,0.320000000000000,0.320000000000000,0.310000000000000,0.310000000000000,0.310000000000000,0.310000000000000,0.310000000000000,0.310000000000000,0.320000000000000,0.320000000000000,0.320000000000000];
                    for sp = 1:46
                        spec_PGA = c1(sp)+c2(sp)*unified_magnitude(j)+c4(sp)*log10(sqrt(unified_distance(i)^2+ho(sp)^2))+ca(sp)*SA+cs(sp)*SS;
                        sigma = sig(sp);
                        specprob_PGA(j,i,sp) = 1-normcdf(log10(Y),spec_PGA,sigma);
                    end
                end
         end
             if specord_check == 1
                if serial_spectral == 1002
                    period_Number = 28; %% VERY IMPORTANT %%
                    ABS_periods = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.075, 0.09, 0.1, 0.12, 0.15, 0.17, 0.2, 0.24, 0.3, 0.36, 0.4, 0.46, 0.5, 0.6, 0.75, 0.85, 1, 1.5, 2, 3, 4, 5];
            ABS_c1  = 6.4;
                 ABS_n  = 2;
                 ABS_a1 = [1.64, 1.64, 1.69, 1.78, 1.87, 1.94, 2.037, 2.1, 2.16, 2.272, 2.407, 2.43, 2.406, 2.293, 2.114, 1.955, 1.86, 1.717, 1.615, 1.428, 1.16, 1.02, 0.828, 0.26, -0.15, -0.69, -1.13, -1.46];
                 ABS_a2  = 0.512;
                 ABS_a12  = [0, 0, 0.0143, 0.0245, 0.028, 0.03, 0.03, 0.03, 0.028, 0.018, 0.005, -0.004, -0.0138, -0.0238, -0.036, -0.046, -0.0518, -0.0594, -0.0635, -0.074, -0.0862, -0.0927, -0.102, -0.12, -0.14, -0.1726, -0.1956, -0.215];
                 ABS_a3 = [-1.145, -1.145, -1.145, -1.145, -1.145, -1.145, -1.145, -1.145, -1.145, -1.145, -1.145, -1.135, -1.115, -1.079, -1.035, -1.0052, -0.988, -0.9652, -0.9515, -0.9218, -0.8852, -0.8648, -0.8383, -0.7721, -0.725, -0.725, -0.725, -0.725];
                 ABS_a13  = 0.17;
                 ABS_a4  = -0.144;
                ABS_c4 = [5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.58, 5.54, 5.5, 5.39, 5.27, 5.19, 5.1, 4.97, 4.8, 4.62, 4.52, 4.38, 4.3, 4.12, 3.9, 3.81, 3.7, 3.55, 3.5, 3.5, 3.5, 3.5];
                 ABS_b5  = [0.7, 0.7, 0.7, 0.71, 0.71, 0.72, 0.73, 0.74, 0.74, 0.75, 0.75, 0.76, 0.77, 0.77, 0.78, 0.79, 0.79, 0.8, 0.8, 0.81, 0.81, 0.82, 0.83, 0.84, 0.85, 0.87, 0.88, 0.89];
                 ABS_b6 = [0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.132, 0.13, 0.127, 0.123, 0.121, 0.118, 0.11, 0.105, 0.097, 0.092, 0.087];
                T = ABS_periods;
                 for sp = 1:28
                  c1 = ABS_c1;
                        n1 = ABS_n;
                        a1 = ABS_a1(sp);
                        a2 = ABS_a2;
                        a12 = ABS_a12(sp);
                        a3 = ABS_a3(sp);
                        a13 = ABS_a13;
                        a4 = ABS_a4;
                        c4 = ABS_c4(sp);
                        b5 = ABS_b5(sp);
                        b6 = ABS_b6(sp);
                        if unified_magnitude(j) <= 5 
                            sigma = b5;
                        elseif unified_magnitude(j) > 5 && unified_magnitude(j) < 7 
                                sigma = b5 - b6 * (unified_magnitude(j) - 5);
                            else
                                sigma = b5 - 2 * b6;
                        end
                        if unified_magnitude(j)<= c1
                            spec_PGA = a1 + a2 * (unified_magnitude(j) - c1) + (a3 + a13 * (unified_magnitude(j) - c1)) * log(sqrt(unified_distance(i) * unified_distance(i) + c4 * c4)) + a12 * (8.5 - unified_magnitude(j)) ^ n1;
                        else
                            spec_PGA = a1 + a4 * (unified_magnitude(j) - c1) + (a3 + a13 * (unified_magnitude(j) - c1)) * log(sqrt(unified_distance(i) * unified_distance(i) + c4 * c4)) + a12 * (8.5 - unified_magnitude(j)) ^ n1;
                        end
                        specprob_PGA(j,i,sp) = 1-normcdf(log(Y),spec_PGA,sigma);
                    end
                end
             end
                
     end
 end
end
if check == 0
    for d = 1:size3(2)
        for i = 1:size1(2)
            for j = 1:size2(2)
                if serial_unified(d) == 5
                    PGA = 6.74+0.859*unified_magnitude(j)-1.8*log(unified_distance(i)+25);
                    sigma = 0.57;
                    prob_PGA(j,i) = 1-normcdf(log(Y*981),PGA,sigma);
                end
                if serial_unified(d) == 2
                    PGA = 0.46+0.63*unified_magnitude(j)-1.1*log10(unified_distance(i));
                    sigma = 0.32;
                    prob_PGA(j,i) = 1-normcdf(log10(Y*981),PGA,sigma);
                end
                if serial_unified(d) == 4
                    PGA = 3.4+0.89*unified_magnitude(j)-1.17*log(unified_distance(i))-0.2*S;
                    sigma = 0.62;
                    prob_PGA(j,i) = 1-normcdf(log(Y*981),PGA,sigma);
                end
                if serial_unified(d) == 3
                    PGA = 0.0159*exp(0.868*unified_magnitude(j))*(unified_distance(i)+0.0606*exp(0.7*unified_magnitude(j)))^(-1.09);
                    sigma = 0.372;
                    prob_PGA(j,i) = 1-normcdf(log(Y*981),log(PGA*981),sigma);
                end
                if serial_unified(d) == 6
                    PGA = 0.0185*exp(1.28*unified_magnitude(j))*(unified_distance(i)+0.147*exp(0.732*unified_magnitude(j)))^(-1.75);
                    sigma = 0.384;
                    prob_PGA(j,i) = 1-normcdf(log(Y),log(PGA),sigma);
                end
                if serial_unified(d) == 7
                    R = sqrt(unified_distance(i)^2+7.3^2);
                    PGA = (-1.02+0.249*unified_magnitude(j)-log10(R)-0.00255*R);
                    sigma = 0.26;
                    prob_PGA(j,i) = 1-normcdf(log10(Y),(PGA),sigma);
                end
                if serial_unified(d) == 8
                    R = sqrt(unified_distance(i)^2+8^2);
                    PGA = 0.43+0.23*(unified_magnitude(j)-6)-log10(R)-0.0027*R;
                    sigma = 0.28;
                    prob_PGA(j,i) = 1-normcdf(log10(Y),(PGA),sigma);
                end
                if serial_unified(d) == 9
                    PGA = 0.41*unified_magnitude(j)-log10(unified_distance(i)+0.032*10^(0.41*unified_magnitude(j)))-0.0034*unified_distance(i)+1.3;
                    sigma = 0.21;
                    prob_PGA(j,i) = 1-normcdf(log10(Y*981),PGA,sigma);
                end
                if serial_unified(d) == 10
                    PGA = -2.501+0.623*unified_magnitude(j)-log(unified_distance(i)+7.28);
                    sigma = 0.506;
                    prob_PGA(j,i) = 1-normcdf(log(Y),PGA,sigma);
                end
            end
        end
        f = f+Nintervals;
    end
end
if Npoint>0
for i = 1:Npoint
    eval(sprintf('a = aP%d',i));
    eval(sprintf('b = bP%d',i));
    eval(sprintf('Mmin = Mminp%d',i));
    eval(sprintf('probMAG = Mprob_point%d',i));
    eval(sprintf('probDIST = Dprob_point%d',i));
    lambda = 10^(a-b*Mmin);
    probPGA = prob_PGA(1+(i-1)*Nintervals:i*Nintervals,1+(i-1)*Nintervals:i*Nintervals);
    recurr = lambda*probMAG*probPGA(1:Nintervals,1);                              %% THE CARE HAS BEEN TAKEN HERE %%
    one_probP(dummy,i) = sum(recurr(:));
    KP = one_probP(:,i);
    eval(sprintf('recurrence_point%d = KP',i));
    if specord_check == 1
        for sp = 1:period_Number
        SpecPGA = specprob_PGA(1+(i-1)*Nintervals:i*Nintervals,1+(i-1)*Nintervals:i*Nintervals,sp);
        recurr_spec = lambda*probMAG*SpecPGA(1:Nintervals,1);
        one_probPspec(dummy,i,sp) = sum(recurr_spec(:));
        KP1 = one_probPspec(:,i,sp);
        eval(sprintf('recurrence_point%d_per%d = KP1',i,sp));
        end
    end
end
else
    one_probP(Nacc,1) = 0;
    if specord_check == 1
    one_probPspec(Nacc,1,1:period_Number) = 0;
    end
end
if Nline>0
for i = 1:Nline
    eval(sprintf('a = aL%d',i));
    eval(sprintf('b = bL%d',i));
    eval(sprintf('Mmin = Mminl%d',i));
    eval(sprintf('probMAG = Mprob_line%d',i));
    eval(sprintf('probDIST = Dprob_line%d',i));
    lambda = 10^(a-b*Mmin);
    probPGA = prob_PGA(1+(i-1)*Nintervals+Npoint*Nintervals:i*Nintervals+Npoint*Nintervals,1+(i-1)*Nintervals+Npoint*Nintervals:i*Nintervals+Npoint*Nintervals);
    recurr = lambda*probMAG*probPGA*probDIST';
    one_probL(dummy,i) = sum(recurr(:));
    KL = one_probL(:,i);
    eval(sprintf('recurrence_line%d = KL',i));
    if specord_check == 1
        for sp = 1:period_Number
        SpecPGA = specprob_PGA(1+(i-1)*Nintervals+Npoint*Nintervals:i*Nintervals+Npoint*Nintervals,1+(i-1)*Nintervals+Npoint*Nintervals:i*Nintervals+Npoint*Nintervals,sp);
        recurr_spec = lambda*probMAG*SpecPGA*probDIST';
        one_probLspec(dummy,i,sp) = sum(recurr_spec(:));
        KL1 = one_probLspec(:,i,sp);
        eval(sprintf('recurrence_line%d_per%d = KL1',i,sp));
        end
    end
end
else
    one_probL(Nacc,1) = 0;
    if specord_check == 1
    one_probLspec(Nacc,1,1:period_Number) = 0;
    end
end
if Narea>0
for i = 1:Narea
    eval(sprintf('a = aA%d',i));
    eval(sprintf('b = bA%d',i));
    eval(sprintf('Mmin = Mmina%d',i));
    eval(sprintf('probMAG = Mprob_area%d',i));
    eval(sprintf('probDIST = Dprob_area%d',i));
    lambda = 10^(a-b*Mmin);
    probPGA = prob_PGA(1+(i-1)*Nintervals+(Npoint+Nline)*Nintervals:i*Nintervals+(Npoint+Nline)*Nintervals,1+(i-1)*Nintervals+(Npoint+Nline)*Nintervals:i*Nintervals+(Npoint+Nline)*Nintervals);
    recurr = lambda*probMAG*probPGA*probDIST';
    one_probA(dummy,i) = sum(recurr(:));
    KA = one_probA(:,i);
    eval(sprintf('recurrence_area%d = KA',i));
    if specord_check == 1
        for sp = 1:period_Number
        SpecPGA = specprob_PGA(1+(i-1)*Nintervals+(Npoint+Nline)*Nintervals:i*Nintervals+(Npoint+Nline)*Nintervals,1+(i-1)*Nintervals+(Npoint+Nline)*Nintervals:i*Nintervals+(Npoint+Nline)*Nintervals,sp);
        recurr_spec = lambda*probMAG*SpecPGA*probDIST';
        one_probAspec(dummy,i,sp) = sum(recurr_spec(:));
        KA1 = one_probAspec(:,i,sp);
        eval(sprintf('recurrence_area%d_per%d = KA1',i,sp));
        end
    end
end
else
    one_probA(Nacc,1) = 0;
    if specord_check == 1
    one_probAspec(Nacc,1,1:period_Number) = 0;
    end
end
dummy = dummy+1;
%%%%%%%%%%%%%%%%%%%%%%%%         %%% GARNERING NECESSARY DATA FOR
%%%%%%%%%%%%%%%%%%%%%%%%         DEAGGREGATION %%%%
if Deagg_check == 1
 for m = 1:agg
     dx = Deagg(m);
     if abs(Y-dx)<0.00001
         if Npoint>0
             for n = 1:Npoint
             Deagg_PGAP(:,n) = prob_PGA(:,n*Nintervals);
             end
         else
             Deagg_PGAP = zeros(Nintervals*Nsources,1);
         end
         if (Nsources-Npoint) > 0
             Deagg_PGAO = prob_PGA(:,(Npoint)*Nintervals+1:Nsources*Nintervals);
         else
           Deagg_PGAO = zeros(Nintervals*Nsources,1);  
         end
         Deagg_PGA = horzcat(Deagg_PGAP,Deagg_PGAO);
         Deagg_PGA(:,~any(Deagg_PGA,1)) = []
         eval(sprintf('Deagg_PGA%d = Deagg_PGA',m));
     end
 end
end
 %%%%%%%%%%%%%%%%%%%%%%
end
Y = Ylow:0.01:Yhigh;
for i = 1:Npoint
    eval(sprintf('recurrence = recurrence_point%d',i));
    figure()
    loglog(Y,recurrence)
    grid on
    grid minor
    xlabel('PGA IN g')
    ylabel('RECUURENCE OR AVERAGE ANNUAL OCCURANCES')
    title(sprintf('RECURRENCE PLOT OF POINT SOURCE %d',i))
    figure()
    loglog(Y,1./recurrence)
    grid on
    grid minor
    xlabel('PGA IN g')
    ylabel('RETURN PERIOD')
    title(sprintf('RETURN PERIOD PLOT OF POINT SOURCE %d',i))
    figure()
    loglog(Y,1-exp(-recurrence*Time_Period))
    grid on
    grid minor
    xlabel('PGA IN g')
    ylabel('PROBABILITY OF ATLEAST ON OCCURANCE IN THE SPECIFIED TIME PERIOD')
    title(sprintf('PROBABILITY OF ATLEAST ONE OCCURANCE (POISSON) PLOT OF POINT SOURCE %d',i))
end
for i = 1:Nline
    eval(sprintf('recurrence = recurrence_line%d',i));
    figure()
    loglog(Y,recurrence)
    grid on
    grid minor
    xlabel('PGA IN g')
    ylabel('RECUURENCE OR AVERAGE ANNUAL OCCURANCES')
    title(sprintf('RECURRENCE PLOT OF LINE SOURCE %d',i))
    figure()
    loglog(Y,1./recurrence)
    grid on
    grid minor
    xlabel('PGA IN g')
    ylabel('RETURN PERIOD')
    title(sprintf('RETURN PERIOD PLOT OF LINE SOURCE %d',i))
    figure()
    loglog(Y,1-exp(-recurrence*Time_Period))
    grid on
    grid minor
    xlabel('PGA IN g')
    ylabel('PROBABILITY OF ATLEAST ON OCCURANCE IN THE SPECIFIED TIME PERIOD')
    title(sprintf('PROBABILITY OF ATLEAST ONE OCCURANCE (POISSON) PLOT OF LINE SOURCE %d',i))
end
for i = 1:Narea
    eval(sprintf('recurrence = recurrence_area%d',i));
    figure()
    loglog(Y,recurrence)
    grid on
    grid minor
    xlabel('PGA IN g')
    ylabel('RECUURENCE OR AVERAGE ANNUAL OCCURANCES')
    title(sprintf('RECURRENCE PLOT OF AREA SOURCE %d',i))
    figure()
    loglog(Y,1./recurrence)
    grid on
    grid minor
    xlabel('PGA IN g')
    ylabel('RETURN PERIOD')
    title(sprintf('RETURN PERIOD PLOT OF AREA SOURCE %d',i))
    figure()
    loglog(Y,1-exp(-recurrence*Time_Period))
    grid on
    grid minor
    xlabel('PGA IN g')
    ylabel('PROBABILITY OF ATLEAST ON OCCURANCE IN THE SPECIFIED TIME PERIOD')
    title(sprintf('PROBABILITY OF ATLEAST ONE OCCURANCE (POISSON) PLOT OF AREA SOURCE %d',i))
end
total = horzcat(one_probP,one_probL,one_probA);
if specord_check == 1
total_spec = horzcat(one_probPspec,one_probLspec,one_probAspec);
end
total = total';
total_recurrence = sum(total);
if specord_check == 1
 for sp = 1:period_Number
     temp_last = total_spec(:,:,sp);
     temp_last = sum(temp_last');
     eval(sprintf('specord_recurr%d = temp_last',sp));
 end
end
figure()
    loglog(Y,total_recurrence)
    grid on
    grid minor
    xlabel('PGA IN g')
    ylabel('RECUURENCE OR AVERAGE ANNUAL OCCURANCES')
    title('TOTAL RECURRENCE OF ALL SOURCES')
    figure()
    loglog(Y,1./total_recurrence)
    grid on
    grid minor
    xlabel('PGA IN g')
    ylabel('RETURN PERIOD')
    title('RETURN PERIOD PLOT OF ALL SOURCES')
    figure()
    loglog(Y,1-exp(-total_recurrence*Time_Period))
    grid on
    grid minor
    xlabel('PGA IN g')
    ylabel('PROBABILITY OF ATLEAST ON OCCURANCE IN THE SPECIFIED TIME PERIOD')
    title('PROBABILITY OF ATLEAST ONE OCCURANCE (POISSON) PLOT OF ALL SOURCES')
    
    %%% DEAGGREGATION %%%
    
    if Deagg_check == 1
    if Npoint>0
    for i = 1:Npoint
        eval(sprintf('m1 = Mminp%d',i));
        eval(sprintf('m2 = Mmaxp%d',i));
        Deagg_distexp = dist1;
        Deagg_magex1(i) = m1;
        Deagg_magex2(i) = m2;
    end
    else
        Deagg_distexp = 0;
        Deagg_magex1 = 0;
        Deagg_magex2 = 0;
    end
    if Nline>0
    for i  = 1:Nline
        eval(sprintf('d1 = Rminline%d',i));
        eval(sprintf('d2 = Rmaxline%d',i));
        eval(sprintf('m1 = Mminl%d',i));
        eval(sprintf('m2 = Mmaxl%d',i));
        Deagg_distex1(i) = d1;
        Deagg_distex2(i) = d2;
        Deagg_magex3(i) = m1;
        Deagg_magex4(i) = m2;
    end
    else
        Deagg_distex1 = 0;
        Deagg_distex2 = 0;
        Deagg_magex3 = 0;
        Deagg_magex4 = 0;
    end
    if Narea>0
    for i = 1:Narea
        eval(sprintf('d1 = Rminarea%d',i));
        eval(sprintf('d2 = Rmaxarea%d',i));
        eval(sprintf('m1 = Mmina%d',i));
        eval(sprintf('m2 = Mmaxa%d',i));
        Deagg_distex3(i) = d1;
        Deagg_distex4(i) = d2;
        Deagg_magex5(i) = m1;
        Deagg_magex6(i) = m2;
    end
    else
    Deagg_magex5 = 0;
    Deagg_magex6 = 0;
    Deagg_distex3 = 0;
    Deagg_distex4 = 0;
    end
    Deagg_distex = horzcat(Deagg_distexp,Deagg_distex1,Deagg_distex2,Deagg_distex3,Deagg_distex4);
    Deagg_magex = horzcat(Deagg_magex1,Deagg_magex2,Deagg_magex3,Deagg_magex4,Deagg_magex5,Deagg_magex6);
    Deagg_distex(Deagg_distex == 0) = [];
    Deagg_magex(Deagg_magex == 0) = [];
    
    for i = 1:agg
        PGA_Deagg = eval(sprintf('Deagg_PGA%d',i));
        for j = 1:Npoint
            PGA_final(1+(j-1)*Nintervals:j*Nintervals,j) = PGA_Deagg(1+(j-1)*Nintervals:j*Nintervals,j);
        end
        for k = 1:(Nsources-Npoint)
            PGA_final(((Npoint*Nintervals+1)+(k-1)*Nintervals):(Npoint+k)*Nintervals,((1+Npoint)+(k-1)*Nintervals):((1+Npoint)+(k-1)*Nintervals+Nintervals-1)) = PGA_Deagg(((Npoint*Nintervals+1)+(k-1)*Nintervals):(Npoint+k)*Nintervals,((1+Npoint)+(k-1)*Nintervals):((1+Npoint)+(k-1)*Nintervals+Nintervals-1));
        end
        eval(sprintf('Deagg_PGAA%d = PGA_final',i));
    end
    if Npoint>0
        for i = 1:Npoint
        recurrence_matrix1(1,i,1) = eval(sprintf('aP%d',i));
        recurrence_matrix1(1,i,2) = eval(sprintf('bP%d',i));
        recurrence_matrix1(1,i,3) = eval(sprintf('Mminp%d',i));
        end
    else
        for i = 1:3
        recurrence_matrix1(1,1,i) = 0;
        end
    end
    if Nline>0
        for i = 1:Nline
        recurrence_matrix2(1,i,1) = eval(sprintf('aL%d',i));
        recurrence_matrix2(1,i,2) = eval(sprintf('bL%d',i));
        recurrence_matrix2(1,i,3) = eval(sprintf('Mminl%d',i));
        end
    else
        for i = 1:3
        recurrence_matrix2(1,1,i) = 0;
        end
    end
    if Narea>0
        for i = 1:Narea
        recurrence_matrix3(1,i,1) = eval(sprintf('aA%d',i));
        recurrence_matrix3(1,i,2) = eval(sprintf('bA%d',i));
        recurrence_matrix3(1,i,3) = eval(sprintf('Mmina%d',i));
        end
    else
        for i = 1:3
        recurrence_matrix3(1,1,i) = 0;
        end
    end
    recurrence_matrix = horzcat(recurrence_matrix1,recurrence_matrix2,recurrence_matrix3);
    recurrence_matrix(recurrence_matrix == 0) = []
    for i = 1:3
        matrix(1,1:Nsources,1) = recurrence_matrix(1,1:Nsources);
        matrix(1,1:Nsources,2) = recurrence_matrix(1,(Nsources+1):2*Nsources);
        matrix(1,1:Nsources,3) = recurrence_matrix(1,(2*Nsources+1):3*Nsources);
    end
    recurrence_matrix = matrix;
    DRmaxN = max(Deagg_distex);
    DRminN = min(Deagg_distex);
    DMmaxN = max(Deagg_magex);
    DMminN = min(Deagg_magex);
    DRmax = max(Deagg_distance);
    DRmin = min(Deagg_distance);
    DMmax = max(Deagg_magnitude);
    DMmin = min(Deagg_magnitude);
    v1 = size(Deagg_distance);
    v2 = size(Deagg_magnitude);
    DdR = (DRmax-DRmin)/Nintervals;
    DdM = (DMmax-DMmin)/Nintervals;
    DdRN = (DRmaxN-DRminN)/Nintervals;
    DdMN = (DMmaxN-DMminN)/Nintervals;
    concatR = find(Deagg_distance == DRmin);
    concatM = find(Deagg_magnitude == DMmin);
    for i = 1:Npoint
        source_control(i) = i;
    end
    for i = 1:(Nsources-Npoint)
        source_control(((1+Npoint)+(i-1)*Nintervals):((1+Npoint)+(i-1)*Nintervals+Nintervals-1)) = i+Npoint;
    end
    o = 1;
    
    for k = 1:agg
        eval(sprintf('PGA_values = Deagg_PGAA%d',k));
        YD = Deagg(k);
        for i = 1:Nintervals%%%%%%%%%%
            
        Rlow = DRmin+(i-1)*DdR;
        Rhigh = DRmin+(i)*DdR;
        Deagg_Ravg(i) = (Rlow+Rhigh)/2;
        RlowN = DRminN+(i-1)*DdRN;
        RhighN = DRminN+(i)*DdRN;
        Deagg_RavgN(i) = (RlowN+RhighN)/2;
        indicesR = find(Deagg_distance<=Rhigh&Deagg_distance>Rlow);
        if i == 1
            indicesR = horzcat(concatR,indicesR);
        end
        for j = 1:Nintervals%%%%%%%%%%
            Mlow = DMmin+(j-1)*DdM;
            Mhigh = DMmin+j*DdM;
            Deagg_Mavg(j) = (Mlow+Mhigh)/2;
            MlowN = DMminN+(j-1)*DdMN;
            MhighN = DMminN+j*DdMN;
            Deagg_MavgN(j) = (MlowN+MhighN)/2;
            indicesM = find(Deagg_magnitude<=Mhigh&Deagg_magnitude>Mlow);
            if j == 1
                indicesM = horzcat(concatM,indicesM);
            end
            controlR = size(indicesR);
            controlM = size(indicesM);
            multiply = 0;
            for h = 1:controlR(2)
                for l = 1:controlM(2)
                    multiply = multiply+PGA_values(indicesM(l),indicesR(h))*Deagg_distprob(indicesR(h))*Deagg_magnprob(indicesM(l))*10.^(recurrence_matrix(1,source_control(indicesR(h)),1)-recurrence_matrix(1,source_control(indicesR(h)),2)*recurrence_matrix(1,source_control(indicesR(h)),3));
                end
            end
            intermediate_values(j,i,k) = multiply;
            o = o+1;
        end
        end
    end
    for i = 1:agg
        Yspec = Deagg(i);
        
        tt = ((Yhigh-Ylow)/0.01+1)-(Yhigh-Yspec)/0.01;
        tt = round(tt);
        sample_space = total_recurrence(tt);
        prob_final = intermediate_values(:,:,i)/sample_space;
        k1 = sum(prob_final,1);
        k2 = sum(k1);
        abs(k2-1)
        eval(sprintf('percentage_contribution%d = prob_final',i));
        figure()
        bar3(prob_final)
        set(gca,'YTickLabel',Deagg_MavgN,'YTick',1:numel(Deagg_MavgN))
        set(gca,'XTickLabel',Deagg_RavgN,'XTick',1:numel(Deagg_RavgN))
        zlabel('probability')
        ylabel('Magnitude')
        xlabel('Distance(Km)')
        title(sprintf('Probability contribution of various distances and magnitudes of acceleration value %d',i))
        if abs(k2-1)>0.2
            i
            display('Error in Deaagregation')
            beep
        else
            display('SUCCESS')
        end
    end
  %  end
 
    %%% UNIFORM HAZARD SPECTRUM %%%
    if specord_check == 1
    ss = size(Hazard_spectrum);
    for j = 1:ss(2)
        for i = 1:period_Number
        eval(sprintf('temp_spec = specord_recurr%d',i))
            unif_matrix(i,j) = interp1(temp_spec,Y,Hazard_spectrum(j));
        end
    end
    for i = 1:ss(2)
        figure()
        plot(T,unif_matrix(:,i));
        xlabel('Time period (Seconds)');
        ylabel('Absolute Spectral accleration in terms of g');
        title(sprintf('Uniform Hazard Spectrum for hazard value %d',i))
        grid minor
    end
    end
    end

    
    
    
    
        
            
            
    
            
    
    
    
        
    
    
    





    
    




    
    
    



             

     
 
 
 
 
 
 
 
 
     
      
 
 
     
     
 
 
 

    












    
    
        
        
    
    

        
  
        
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    



