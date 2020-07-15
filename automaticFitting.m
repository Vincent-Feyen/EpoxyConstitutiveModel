close all
clear
clc

%define sample size (area and volume)
A0 = 0.006^2*pi();
V0 = 0.012*A0;

%choose the material

%for the NTPT material
% CSRdataNames = {'NTPT_CSR_1','NTPT_CSR_2','NTPT_CSR_3'}; %CSR 0.1 0.01 0.001
% HALdataNames = {'NTPT_HAL_1','NTPT_HAL_2','NTPT_HAL_3','NTPT_HAL_4','NTPT_HAL_5','NTPT_HAL_6','NTPT_HAL_7','NTPT_HAL_8','NTPT_HAL_9'};
% CSRrates = [0.1;0.01;0.001];
% timeLimitLow = 0;
% timeLimitHigh = 2000;
% filename = 'NTPT_e0_manual.xlsx';

%for the RTM material
CSRdataNames = {'RTM_CSR_1','RTM_CSR_2','RTM_CSR_3','RTM_CSR_4'}; %CSR 0.1 0.01 0.001 0.0001
CSRrates = [0.1;0.01;0.001;0.0001];
HALdataNames = {'RTM_HAL_1','RTM_HAL_2','RTM_HAL_3','RTM_HAL_4','RTM_HAL_5','RTM_HAL_6','RTM_HAL_7','RTM_HAL_8','RTM_HAL_9'};
timeLimitLow = 0; %startTime of the test
timeLimitHigh = 2000; %stopTime of the test
filename = 'RTM_e0_manual.xlsx';

%load start strains and loads for the HAL data
sheet = 1;
xlRange1 = 'E2:E10';%contains the start strains
xlRange2 = 'K2:K10';%contains the applied load
startStrains = xlsread(filename,sheet,xlRange1);
loads = xlsread(filename,sheet,xlRange2);

%load data for CSR
for i=1:length(CSRdataNames)  
    filename = [CSRdataNames{i} '.csv'];
    CSRdata{i} = csvread(filename);
    CSRdata{i} = [CSRdata{i};[0 0]]; 
    CSRdata{i} = sortrows(CSRdata{i});
    indices = CSRdata{i}(:,1)<0;
    CSRdata{i}(indices,1) = 0.00001;
    dataPoints = size(CSRdata{i},1);
    CSRdata{i} = [CSRdata{i} ones(dataPoints,1).*CSRrates(i)]; %add strain rates
end

%visually check the CSR data
figure
hold on 
for i=1:length(CSRdataNames) 
plot(CSRdata{i}(:,1),CSRdata{i}(:,2))
end
ylim([0 180])
%%
%load data for HAL
for i=1:length(HALdataNames)  
    filename = [HALdataNames{i} '.csv'];
    HALdata{i} = csvread(filename);
    HALdata{i} = sortrows(HALdata{i});
    indices = HALdata{i}(:,1)<0;
    HALdata{i}(indices,1) = 0.01;
end
%visually check the HAL data
figure
hold on 
for i=1:length(HALdataNames) 
plot(HALdata{i}(:,1),HALdata{i}(:,2))
end

%filter HAL data if the fit is not good (reduce the selected time)
for i=1:length(HALdataNames) 
    startRow(i) = 1;
    dataRows = size(HALdata{i},1);
    for r=1:dataRows 
        if(HALdata{i}(r,1) < timeLimitLow && r > startRow(i))
        startRow(i) = r;
        end
    end
end

%fit the individual HAL curves to obtain the local strain rate,
%using a simple power law, visually check the result

for i=1:length(HALdataNames)
    [xData, yData] = prepareCurveData(HALdata{i}(startRow(i):end,1),HALdata{i}(startRow(i):end,2));
    ft = fittype( 'a*x^b', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Robust = 'Bisquare';
    opts.MaxFunEvals = 60000;
    opts.MaxIter = 40000;
    opts.Lower = [0 0];
    opts.Upper = [Inf Inf];
    opts.StartPoint = [0.0001 0.2];
    [fitresult, gof] = fit( xData, yData, ft, opts );
    coeffvals = coeffvalues(fitresult);
    alpha(i) = coeffvals(1);
    beta(i) = coeffvals(2);
    figure
    h = plot( fitresult, xData, yData);
    h(1).LineWidth = 2;
    h(2).LineWidth = 2;
    h(1).MarkerSize = 20;
    h(2).MarkerSize = 20;
end
%store values for time-dependent creep law
test_alpha = alpha;
test_beta = beta;
%Calculate and add the strain rates to the HAL data
for i=1:length(HALdataNames)
    dataPoints = size(HALdata{i},1);
    totalStrain = zeros(dataPoints,1);
    trueStress = zeros(dataPoints,1);
    strainRate = zeros(dataPoints,1);
    for r=1:dataPoints
        strainRate(r) = alpha(i)*beta(i)*HALdata{i}(r,1)^(beta(i)-1);
        totalStrain(r) = alpha(i)*HALdata{i}(r,1)^beta(i)+startStrains(i);
        A = A0*exp(totalStrain(r));%note that the true strain is compressive
        trueStress(r) = loads(i)/A/1000;
    end
    HALdata{i} = [HALdata{i} totalStrain trueStress strainRate];
end
%%keep the first three for powerlaw fit
t_p = [HALdata{1}(:,1);HALdata{2}(:,1);HALdata{3}(:,1);HALdata{4}(:,1)];
s_p = [HALdata{1}(:,4);HALdata{2}(:,4);HALdata{3}(:,4);HALdata{4}(:,4)];
z_p = [HALdata{1}(:,2);HALdata{2}(:,2);HALdata{3}(:,2);HALdata{4}(:,2)];
%%
%visually check the strain rate range in which we can use the HAL data
figure
hold on
for i=1:length(HALdataNames)
    scatter3(HALdata{i}(:,3),log10(HALdata{i}(:,5)),HALdata{i}(:,4))
end
for i=1:length(CSRdataNames)
    scatter3(CSRdata{i}(:,1),log10(CSRdata{i}(:,3)),CSRdata{i}(:,2))
end
ylim([-6 0])
zlim([0 160])

%select the range for the HAL data
eLevels = [0.001;0.0001;0.00001];

%create stress strain curves from the HAL data at the selected strain rates
HALcurves = cell(1,length(eLevels));
for i=1:length(eLevels)
    strainRate = eLevels(i); %[1/s]
    for r=1:length(HALdataNames)
        time(r) = (strainRate/alpha(r)/beta(r))^(1/(beta(r)-1));
        if(time(r) > timeLimitLow && time(r) < timeLimitHigh)
           totalStrain = alpha(r)*time(r)^beta(r)+startStrains(r);
           A = A0*exp(totalStrain);%note that the true strain is compressive
           trueStress = loads(r)/A/1000;
           HALcurves{i} = [HALcurves{i};[totalStrain trueStress strainRate]];
        end
    end  
end
%visually check the obtained curves
figure
hold on
for i=1:length(HALcurves)
    scatter3(HALcurves{i}(:,1),log10(HALcurves{i}(:,3)),HALcurves{i}(:,2))
end

%define the data to fit the coefficients
strainRateDataSize = length(CSRdata)+length(eLevels);
strainRateData = cell(1,strainRateDataSize);

%add CSR data
for i=1:length(CSRdata)
    strainRateData{1,i} = CSRdata{1,i};
end

%add HAL data
for i=1:length(HALcurves)
    p = i+length(CSRdata);
    strainRateData{1,p} = HALcurves{1,i};
end

%fit the prefactors
finalCoefficients = zeros(1,9);%a,b,c,d,e,f,g,h,j
testCoefficients = zeros(strainRateDataSize,6);%reused every fitting

%fit c,f and j individually
for i = 1:1:strainRateDataSize
    xfit_psr = strainRateData{i}(:,1);
    yfit_psr = strainRateData{i}(:,2);
    [xData, yData] = prepareCurveData( xfit_psr, yfit_psr );
    funct = "abs(a)*exp(-abs(b)*x)-abs(a+e)*exp(-abs(d)*x)+abs(e)*exp(abs(f)*x);";
    ft = fittype(char(funct), 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Robust = 'Bisquare';
    opts.Lower = [0 0 0 0 0];
    opts.MaxFunEvals = 600000;
    opts.MaxIter = 400000;
    opts.StartPoint = [0.01 0.01 0.01 0.01 0.01];
    opts.Upper = [Inf Inf Inf Inf Inf];
    [fitresult, gof] = fit( xData, yData, ft, opts );
    coeffvals = coeffvalues(fitresult);
    figure
    h = plot( fitresult, xData, yData );
    h(1).LineWidth = 2;
    h(2).LineWidth = 2;
    h(1).MarkerSize = 20;
    h(2).MarkerSize = 20;
    testCoefficients(i,1) = coeffvals(1);
    testCoefficients(i,2) = coeffvals(2);
    testCoefficients(i,3) = coeffvals(1)+coeffvals(5);
    testCoefficients(i,4) = coeffvals(3);
    testCoefficients(i,5) = coeffvals(4);
    testCoefficients(i,6) = coeffvals(5);
end

%the exact shape of the curve is more pronounced by the CSR curves,
%as a result we will only use this data to fit the c,f and j coefficient

%visually check the spread of these coefficients 
figure
plot(testCoefficients(1:length(CSRdata),2),'o')
figure
plot(testCoefficients(1:length(CSRdata),4),'o')
figure
plot(testCoefficients(1:length(CSRdata),6),'o')

%calculate the average and refit everything
finalCoefficients(1,3) = mean(testCoefficients(1:length(CSRdata),2));
finalCoefficients(1,6) = mean(testCoefficients(1:length(CSRdata),4));
finalCoefficients(1,9) = mean(testCoefficients(1:length(CSRdata),6));
%%
%refit pre-factors using the average value coefficients of c,f and j
%and store the strain rate as well
for i = 1:1:strainRateDataSize
    xfit_psr = strainRateData{i}(:,1);
    yfit_psr = strainRateData{i}(:,2);
    [xData, yData] = prepareCurveData( xfit_psr, yfit_psr );
    funct =  strcat("abs(a)*exp(-abs(",string(finalCoefficients(1,3)),")*x)-abs(a+g)*exp(-abs(",string(finalCoefficients(1,6)),")*x)+abs(g)*exp(abs(",string(finalCoefficients(1,9)),")*x);");
    ft = fittype(char(funct), 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Robust = 'Bisquare';
    opts.Lower = [0 0];
    opts.MaxFunEvals = 600000;
    opts.MaxIter = 400000;
    opts.StartPoint = [0.01 0.01];
    opts.Upper = [Inf Inf];
    [fitresult, gof] = fit( xData, yData, ft, opts );
    coeffvals = coeffvalues(fitresult);
    figure
    h = plot( fitresult, xData, yData );
    h(1).LineWidth = 2;
    h(2).LineWidth = 2;
    h(1).MarkerSize = 20;
    h(2).MarkerSize = 20;
    testCoefficients(i,1) = coeffvals(1);
    testCoefficients(i,2) = coeffvals(2);
    testCoefficients(i,3) = log10(strainRateData{i}(1,3)); %this variable stores the corresponding strain rate
end
%%
%fit the linear dependence of the prefactor
for i=1:2
    [xData, yData] = prepareCurveData(testCoefficients(:,3),testCoefficients(:,i));
    ft = fittype( 'a*x+b', 'independent', 'x', 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
    opts.Display = 'Off';
    opts.Robust = 'Bisquare';
    opts.MaxFunEvals = 60000;
    opts.MaxIter = 40000;
    opts.Lower = [-Inf -Inf];
    opts.Upper = [Inf Inf];
    opts.StartPoint = [0.01 0.01];
    [fitresult, gof] = fit( xData, yData, ft, opts );
    coeffvals = coeffvalues(fitresult);
    a_values(i) = coeffvals(1);
    b_values(i) = coeffvals(2);
    figure
    h = plot( fitresult, xData, yData );
    h(1).LineWidth = 2;
    h(2).LineWidth = 2;
    h(1).MarkerSize = 20;
    h(2).MarkerSize = 20;
end
%%
%assign the values to the corresponding coefficients
finalCoefficients(1,1) = a_values(1);
finalCoefficients(1,2) = b_values(1);
finalCoefficients(1,4) = a_values(1)+a_values(2);
finalCoefficients(1,5) = b_values(1)+b_values(2);
finalCoefficients(1,7) = a_values(2);
finalCoefficients(1,8) = b_values(2);
%%
%transfer the coefficients for sake of brevity
a = finalCoefficients(1,1);
b = finalCoefficients(1,2);
c = finalCoefficients(1,3);
d = finalCoefficients(1,4);
e = finalCoefficients(1,5);
f = finalCoefficients(1,6);
g = finalCoefficients(1,7);
h = finalCoefficients(1,8);
j = finalCoefficients(1,9);
%%
%construct fitted stress strain curves, and overplot with the scatter
maxStrain = max(strainRateData{1}(:,1)); %the maximum strain of the data
colorMatrix = ['r';'b';'c';'k';'g']; %colors for the plots
figure
hold on
for i = 1:1:strainRateDataSize
    xfit_psr = strainRateData{i}(:,1);
    yfit_psr = strainRateData{i}(:,2);
    e_v = [0:0.01:maxStrain]'; %strain vector
    e_r = log10(strainRateData{i}(1,3)); %strain rate
    stress = (a.*e_r+b).*exp(-c.*e_v)-(d.*e_r+e).*exp(-f.*e_v)+(g.*e_r+h).*exp(j.*e_v);
    if(i <= (strainRateDataSize-length(eLevels))) %CSR data
        scatter(xfit_psr,yfit_psr,[],colorMatrix(-e_r),'filled');
    else %HAL data
        scatter(xfit_psr,yfit_psr,[],colorMatrix(-e_r));
    end
    plot(e_v,stress,'color',colorMatrix(-e_r));
end
xlabel("Total true strain [/]")
ylabel("True stress [MPa]")
title("NTPT - CSR and HAL tests")
xlim([0 0.8])
ylim([0 150])
%
%reconstruct the fitted time, strain curves
dt = 1;%[s]
time = [timeLimitLow:1:timeLimitHigh];
RECdata = cell(1,length(HALdata)); %size of the cell
colorMatrix = rand(length(HALdata),3)./1.5; %colors for the plots
figure
hold on
for r = 1:1:length(HALdata)
    RECdata{r}(1,1) = timeLimitLow;%time
    RECdata{r}(1,2) = startStrains(r);%true total plastic strain
    A = A0*exp(RECdata{r}(1,2));
    RECdata{r}(1,3) = loads(r)/A/1000;%true stress
    for i = 1:1:(length(time)-1)
        t_temp = RECdata{r}(i,1);%local total time
        e_temp = RECdata{r}(i,2);%local total strain
        s_temp = RECdata{r}(i,3);%local total stress
        e_r_temp = 10^((s_temp - b*exp(-c*e_temp)+e*exp(-f*e_temp)-h*exp(j*e_temp))/(a*exp(-c*e_temp)-d*exp(-f*e_temp)+g*exp(j*e_temp)));%local strain rate
        if(e_r_temp > 0.01)
           e_r_temp = 0.01;
        end
        t_new = t_temp+dt;
        e_new = e_temp+e_r_temp*dt;
        s_new = loads(r)/(A0*exp(e_new)*1000);
        RECdata{r} = [RECdata{r};[t_new e_new s_new]];
    end
    scatter(HALdata{r}(:,1),HALdata{r}(:,3),[],colorMatrix(r,:));
    plot(RECdata{r}(:,1),RECdata{r}(:,2),'color',colorMatrix(r,:));
end
xlim([0 2000])
xlabel("time [s]")
ylabel("Total true strain [/]")
title("NTPT - HAL tests")

figure
hold on
for r = 1:1:length(HALdata)
    scatter(HALdata{r}(:,3),HALdata{r}(:,4),[],colorMatrix(r,:));
    plot(RECdata{r}(:,2),RECdata{r}(:,3),'color',colorMatrix(r,:));
end
xlim([0 0.6])
ylim([0 150])
xlabel("Total true strain [/]")
ylabel("Total true stress [MPa]")
title("NTPT - HAL tests")
