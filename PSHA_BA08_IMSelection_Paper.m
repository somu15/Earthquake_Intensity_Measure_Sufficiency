function [Deagg_PGA1,recurr] =  PSHA_BA08_IMSelection_Paper(Y,sno,Tp,fault_type,Vs)
format long
unified_distance = [18.4469817918902,23.8294870613328,29.2119923307754,34.5944976002181,39.9770028696607,45.3595081391033,50.7420134085459,56.1245186779886,61.5070239474312,66.8895292168738,72.2720344863164,77.6545397557591,83.0370450252017,88.4195502946443,93.8020555640869,99.1845608335296,104.567066102972,109.949571372415,115.332076641857,120.714581911300];
prob_dist = [0.165384232234926,0.0849849061411171,0.0752046374992244,0.0710276342728635,0.0687671431863026,0.0673827200014182,0.0486311667340931,0.0329116884732186,0.0326774325357362,0.0325009246625283,0.0323644360723405,0.0322566128096217,0.0321698906316319,0.0320990625647778,0.0320404440673310,0.0319913655975675,0.0319498530073355,0.0319144199832119,0.0318839296184171,0.0318574999063374];
unified_magnitude = [3.12500000000000,3.37500000000000,3.62500000000000,3.87500000000000,4.12500000000000,4.37500000000000,4.62500000000000,4.87500000000000,5.12500000000000,5.37500000000000,5.62500000000000,5.87500000000000,6.12500000000000,6.37500000000000,6.62500000000000,6.87500000000000,7.12500000000000,7.37500000000000,7.62500000000000,7.87500000000000];
prob_mag = [0.365888933346237,0.230841153441669,0.145638835356232,0.0918842678079041,0.0579702429639990,0.0365737154953527,0.0230745395696497,0.0145578421317075,0.00918461522891983,0.00579461956930919,0.00365585439521673,0.00230649677673643,0.00145518032338926,0.000918080525816079,0.000579221583974965,0.000365433787024558,0.000230553999356629,0.000145457668411390,9.17699686785742e-05,5.78981310730747e-05];
%Y = 5.05;
lambda = 10^(2-0.8*3);
intervals = 20;
eps_max = 4;
eps_min = -eps_max;
global BA08_check
BA08_check = 1;
deps = (eps_max-eps_min)/(intervals-1);
eps = eps_min:deps:eps_max;
prob_eps = normpdf(eps);
req_pga = 1:5:1000;%0.01:0.01:15;
number = 0;
for k = 1:max(size(Y))
    for l = 1:intervals
for j = 1:intervals
    for i = 1:intervals
        if BA08_check == 1
%fprintf(   'S.no    Ground motion intensity measure\n')
                    %names = {'1'  'PGV';'2'  'PGA';'3'  'Sa'}; 
                    %disp(names)
                    %prompt = 'Enter the appropriate serial number';
                    %sno = input(prompt);
                    if sno == 3
                        %prompt = 'Enter the time period for Sa between 0.01 - 10 seconds';
                        %Tp = input(prompt);
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
                    %prompt = 'What is the type of fault? (SS 0, RS 1)';
                    %fault_type = input(prompt);
                    if fault_type == 0
                        SS = 1;
                        RS = 0;
                    else
                        RS = 1;
                        SS = 0;
                    end
%prompt = 'What is the shear wave velocity in m/s ?';
                    %Vs = input(prompt);
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
prob_PGA(j,i,l) = 1-normcdf(log(Y(k)),PGA+eps(l)*sigma,sigma);
matrix(j,i,l) = prob_mag(j)*prob_dist(i)*prob_eps(l);
    end
end
    end
recurr(k) = lambda*sum(sum(sum(prob_PGA.*matrix)))*deps;
if any(abs(Y-(req_pga))<1e-10),
    number = number + 1;
    temp = lambda*prob_PGA.*matrix*deps/recurr(k);
    Deagg_PGA1 = temp;
    %eval(sprintf('Deagg_PGA%d = temp',number));
end
end
recurr = recurr';
end