
%%%%%%%%% Script in "A unified representation of the temperature dependences of soil microbial growth and respiration, Communications Earth & Environment, 2025".
%%%%%%%%% Code to generate model results and figures by Albert C. Brangarí, University of Amsterdam, IBED, 2025. E-mail: a.carlesbrangari@uva.nl

%%%%%%%%% Important: To run this script, external data need to be downloaded and saved in the same folder where this code is run.
%%%%%%%%% Please download the 'SoilTemp’ database (Hoogen et al., 2022, ﻿Global Soil Bioclimatic variables at 30 arc second resolution (Version 2) [Data set]): 'SBIO1_0_5cm_Annual_Mean_Temperature.tif', 'SBIO2_0_5cm_mean_diurnal_range.tif' and 'SBIO7_0_5cm_annual_range.tif'

clearvars;
close all;

runX = 0; %0: run calibration across the European dataset without model comparison mapping (necessary if 'SoilTemp' dataset has not been dowloaded), %1: run the whole scriot (needs 'SoilTemp')

posTdep = 1;
posSite = 1;
nSites = 70;
ttLL = linspace(0,45,10);
regT = linspace(-20,60,500);
regT_Zwi = linspace(0,50,500); 
regT_Hib = linspace(-5,25,500);
regT_Fan = linspace(-5,50,500);
regT_Cor = linspace(-5,50,500);
sizeF = 18;
TolS = 1e-12;
MxIt = 500000;
t_days = linspace(1,365,365);
t_hours = linspace(1,24,24);
t_time = linspace(0,365*24,365*24);

opts = optimset('MaxIter',MxIt,'MaxFunEval',MxIt,'TolFun',TolS,'TolX',TolS,'TolConSQP',TolS,'Display','off'); % solver options definition
optsCO = optimoptions('fmincon','MaxIterations',MxIt,'MaxFunctionEvaluations',MxIt,'OptimalityTolerance',TolS,'StepTolerance',TolS,'ConstraintTolerance',TolS,'Display','off');
optsS = optimoptions(@fsolve,'Algorithm','trust-region-dogleg','MaxIter',MxIt,'MaxFunctionEvaluations',MxIt,'TolFun',TolS,'TolX',TolS,'Display','off');

% Model definition
funLin = @(x_,P) max(P(1)*(x_+P(2)),0);
funRatkowskyL = @(x_,P) max(P(1)*(x_+P(2)).*(1-exp(P(3)*(x_-P(4)))),0).^2;
funRatkowskyLSQ = @(x_,P) max(P(1)*(x_+P(2)).*(1-exp(P(3)*(x_-P(4)))),0);
funRatkowskyModL = @(x_,P) max(P(1)*(x_+P(2)).*(1+exp(P(3)*(x_-P(4)))),0).^2;
k = 1.38064852e-23; h = 6.62607015e-34; R = 8.314463; Tr = 273.15;
funArrhenius = @(x_,P) P(1)*exp(-P(2)/R./(x_+Tr));
funMMRT = @(x_,P) k/h.*x_.*exp((-P(1)+P(3).*(x_-P(4)))./(R.*x_)+(P(2)-P(3).*(log(x_)-log(P(4))))/R);
funCwX = @(x_,P) (P(3)+P(4)*exp(-P(5)*(1-x_/P(6))/R./x_)) ./ (1+exp(-P(5)*(1-x_/P(6))/R./x_));
funMMRTmod = @(x_,P) (k/h*x_) .* exp( -(P(1) - x_*P(2) + funCwX(x_,P).*(x_-P(7)-x_.*log(x_/P(7))) ) /R./x_);

% Model calibration and plotting
dataTdep = readtable('EuropeanGradientData.xlsx','Sheet',1,'Range','A1:K761','Format','auto','VariableNamingRule','preserve'); %Data from Cruz-Paredes et al., 2023, ﻿Applied and environmental microbiology.
dataSite = readtable('EuropeanGradientData.xlsx','Sheet',2,'Range','A1:AM73','Format','auto','VariableNamingRule','preserve');
if runX == 1
    [MAST_World,R_World] = readgeoraster('SBIO1_0_5cm_Annual_Mean_Temperature.tif');
    MAST_World = double(MAST_World);
    [DiurnalRange_Matrix,~] = readgeoraster('SBIO2_0_5cm_mean_diurnal_range.tif');
    [SeasonalRange_Matrix,~] = readgeoraster('SBIO7_0_5cm_annual_range.tif');
    [MAST_Europe,R_Europe] = geocrop(MAST_World,R_World,[34,72],[-13,30]);
end

for nn = 1:nSites
    if nn == 3
       posTdep = posTdep + 10;
    end
    if nn == 25
       posTdep = posTdep + 10;
       posSite = posSite + 1;
    end
    [nn dataTdep(posTdep,1).Variables dataSite(posSite,1).Variables]

    vecT = dataTdep(posTdep:posTdep+9,2).Variables;
    vecV(1,:) = dataTdep(posTdep:posTdep+9,7).Variables + dataTdep(posTdep:posTdep+9,8).Variables;  %Microbial growth
    vecV(2,:) = dataTdep(posTdep:posTdep+9,9).Variables;  %Respiration
    MAT(nn) = dataSite(posSite,14).Variables;
    lon(nn) = dataSite(posSite,9).Variables;
    lat(nn) = dataSite(posSite,10).Variables;
    soc(nn) = dataSite(posSite,5).Variables*0.55;

    if runX == 1
        [row(nn),col(nn)] = geographicToDiscrete(R_World,lat(nn),lon(nn));
        [rowX(nn),colX(nn)] = geographicToDiscrete(R_Europe,lat(nn),lon(nn));
        MAST(nn) = MAST_Europe(rowX(nn),colX(nn));
        A_daily(nn) = double(DiurnalRange_Matrix(row(nn),col(nn))/2);
        A_seasonal(nn) = double(SeasonalRange_Matrix(row(nn),col(nn))/2);
        seasonal_lag = 365/12*8;  % Seasonal phase lag (days)
        daily_lag = 14;  % Daily phase lag (hours) (ex: 14)
        Temp = [];
        for dd = 1:length(t_days)
            season_comp = A_seasonal(nn) * cos((2*pi/365) * (t_days(dd)-seasonal_lag));
            day_comp = A_daily(nn) * cos((2*pi/24)*(t_hours-daily_lag));

            Temp = [Temp MAST(nn) + season_comp + day_comp];
        end
    end
    
    posTdep = posTdep + 10;  
    posSite = posSite + 1;
    
    for ii = 1:2
        nanIn = isnan(vecV(ii,1:6));
        yy = vecV(ii,~nanIn);
        tt = vecT(~nanIn);
        [regX_aux,PsolX(ii,1:2),~] = Regression_Tmodel(tt,sqrt(yy),regT,funLin,[0.01 10],opts);
        regX(ii,:) = regX_aux.^2;

        nanIn = isnan(vecV(ii,:));
        yyL = vecV(ii,~nanIn);
        ttL = vecT(~nanIn); 

        posOpt = find(ttL<25);   
        yyOpt = yyL(posOpt);
        ttOpt = ttL(posOpt);

        if ii == 1
            figure(99+nn);  set(gcf,'Position',get(0,'Screensize')); set(gcf,'color','w'); subplot(3,2,[1 3]); plot(ttL,sqrt(yyL),'o','Color',[0.6 0.6 0.6],'LineWidth',2,'MarkerSize',10); hold on;
        
            ParMG_Rat(nn,:) = [PsolX(1,:) NaN NaN];
            funRatkowskySQ = @(x_,P) max(PsolX(1,1)*(x_+PsolX(1,2)).*(1-exp(P(1)*(x_-P(2)))),0);
            [regMG_aux,ParMG_Rat(nn,3:4),err] = Regression_Tmodel(ttL,sqrt(yyL),regT,funRatkowskySQ,[0.1 50],opts);
            regMG_Rat = regMG_aux.^2;
            [regMG_RatX,ParMG_RatX(nn,1:4),err] = Regression_Tmodel(ttL,yyL,regT,funRatkowskyL,ParMG_Rat(nn,:),opts);
            ParMG_RatM(nn,:) = ParMG_Rat(nn,:); regMG_RatM = regMG_Rat;
            [regMG_Arr,ParMG_Arr(nn,:),errA] = Regression_Tmodel(ttOpt,yyOpt,regT,funArrhenius,[30 10],opts);
            [regMG_MMRT,ParMG_MMRT(nn,:),errMMRT] = Regression_Tmodel(ttL+Tr,yyL,regT+Tr,funMMRT,[1e5 10 1e4 290],opts);
            [regMG_MMRTmod,ParMG_MMRTmod(nn,:),errMMRTmod] = Regression_Tmodel(ttL+Tr,yyL,regT+Tr,funMMRTmod,[1e5 -100 -2e3 -2e4 5e5 310 300],opts);
            ParMG_RatMX(nn,:) = ParMG_RatX(nn,:); regMG_RatMX = regMG_RatX;
            
            [MG_R2(nn,1),MG_RMSE(nn,1),MG_AIC(nn,1),MG_BIC(nn,1)] = funModEval(yyL/max(yyL),funArrhenius(ttL',ParMG_Arr(nn,:))/max(yyL),2);
            [MG_R2(nn,2),MG_RMSE(nn,2),MG_AIC(nn,2),MG_BIC(nn,2)] = funModEval(yyL/max(yyL),funMMRT(ttL'+Tr,ParMG_MMRT(nn,:))/max(yyL),4);
            [MG_R2(nn,3),MG_RMSE(nn,3),MG_AIC(nn,3),MG_BIC(nn,3)] = funModEval(yyL/max(yyL),funRatkowskyL(ttL',ParMG_Rat(nn,:))/max(yyL),4);
            [MG_R2(nn,4),MG_RMSE(nn,4),MG_AIC(nn,4),MG_BIC(nn,4)] = funModEval(yyL/max(yyL),funRatkowskyL(ttL',ParMG_RatX(nn,:))/max(yyL),4);
            [MG_R2(nn,5),MG_RMSE(nn,5),MG_AIC(nn,5),MG_BIC(nn,5)] = funModEval(yyL/max(yyL),funRatkowskyL(ttL',ParMG_RatM(nn,:))/max(yyL),4);
            [MG_R2(nn,6),MG_RMSE(nn,6),MG_AIC(nn,6),MG_BIC(nn,6)] = funModEval(yyL/max(yyL),funRatkowskyL(ttL',ParMG_RatMX(nn,:))/max(yyL),4);
            [MG_R2(nn,7),MG_RMSE(nn,7),MG_AIC(nn,7),MG_BIC(nn,7)] = funModEval(yyL/max(yyL),funMMRTmod(ttL'+Tr,ParMG_MMRTmod(nn,:))/max(yyL),7);

            yy_MG = yyL/max(yyL);
            val1 = funArrhenius(ttL',ParMG_Arr(nn,:))/max(yyL);
            val2 = funMMRT(ttL'+Tr,ParMG_MMRT(nn,:))/max(yyL);
            val3 = funRatkowskyL(ttL',ParMG_Rat(nn,:))/max(yyL);
            val4 = funRatkowskyL(ttL',ParMG_RatX(nn,:))/max(yyL);
            val5 = funRatkowskyL(ttL',ParMG_RatM(nn,:))/max(yyL);
            val6 = funRatkowskyL(ttL',ParMG_RatMX(nn,:))/max(yyL);
            val7 = funMMRTmod(ttL'+Tr,ParMG_MMRTmod(nn,:))/max(yyL);
            
            GbiasArr(nn,:) = (funArrhenius(vecT',ParMG_Arr(nn,:))-vecV(ii,:))/std(yyL);
            GbiasMMRT(nn,:) = (funMMRT(vecT'+Tr,ParMG_MMRT(nn,:))-vecV(ii,:))/std(yyL);
            GbiasMMRTmod(nn,:) = (funMMRTmod(vecT'+Tr,ParMG_MMRTmod(nn,:))-vecV(ii,:))/std(yyL);
            GbiasRat(nn,:) = (funRatkowskyL(vecT',ParMG_Rat(nn,:))-vecV(ii,:))/std(yyL);
            GbiasRatX(nn,:) = (funRatkowskyL(vecT',ParMG_RatX(nn,:))-vecV(ii,:))/std(yyL);
            GbiasRatMX(nn,:) = (funRatkowskyL(vecT',ParMG_RatMX(nn,:))-vecV(ii,:))/std(yyL);
            
            funTopt = @(x_) x_ - ParMG_RatM(nn,4) + 1/ParMG_RatM(nn,3)*log(1+ParMG_RatM(nn,3)*(x_+ParMG_RatM(nn,2)));        
            [ToptMG(nn),~,~,~]  = fsolve(funTopt,35,optsS);
            funTopt = @(x_) x_ - ParMG_RatMX(nn,4) + 1/ParMG_RatMX(nn,3)*log(1+ParMG_RatMX(nn,3)*(x_+ParMG_RatMX(nn,2)));        
            [ToptMGX(nn),~,~,~]  = fsolve(funTopt,35,optsS);
            
            if runX == 1
                if ttL(end) > ParMG_RatMX(nn,4)
                    intD = interp1([-ParMG_RatMX(nn,2) ttL(1:end-1)' ParMG_RatMX(nn,4)],[0 yyL(1:end-1) 0],regT,'pchip');
                else
                    intD = interp1([-ParMG_RatMX(nn,2) ttL' ParMG_RatMX(nn,4)],[0 yyL 0],regT,'pchip');
                end
                intD((regT<=-ParMG_RatMX(nn,2)) | (regT>=ParMG_RatMX(nn,4))) = 0;
                y_interp = smooth(regT,intD,0.2,'loess');
                yy_interp = y_interp;
                yy_interp((regT<=-ParMG_RatMX(nn,2)) | (regT>=ParMG_RatMX(nn,4))) = 0;

                MG_data = interp1(regT,yy_interp,Temp,'linear');
                CumMG_data(nn) = trapz(t_time,MG_data);
                TclimCh = 5.8;
                MG_data_fut = interp1(regT,yy_interp,Temp+TclimCh,'linear');
                CumMG_data_fut(nn) = trapz(t_time,MG_data_fut);

                CumMG_Arr(nn) = trapz(t_time,funArrhenius(Temp,ParMG_Arr(nn,:)));
                CumMG_MMRT(nn) = trapz(t_time,funMMRT(Temp+Tr,ParMG_MMRT(nn,:)));
                CumMG_MMRTmod(nn) = trapz(t_time,funMMRTmod(Temp+Tr,ParMG_MMRTmod(nn,:)));
                CumMG_RatX(nn) = trapz(t_time,funRatkowskyL(Temp,ParMG_RatX(nn,:)));
                CumMG_RatMX(nn) = trapz(t_time,funRatkowskyL(Temp,ParMG_RatMX(nn,:)));
                CumMG_Arr_fut(nn) = trapz(t_time,funArrhenius(Temp+TclimCh,ParMG_Arr(nn,:)));
                CumMG_MMRT_fut(nn) = trapz(t_time,funMMRT(Temp+Tr+TclimCh,ParMG_MMRT(nn,:)));
                CumMG_MMRTmod_fut(nn) = trapz(t_time,funMMRTmod(Temp+Tr+TclimCh,ParMG_MMRTmod(nn,:)));
                CumMG_RatX_fut(nn) = trapz(t_time,funRatkowskyL(Temp+TclimCh,ParMG_RatX(nn,:)));
                CumMG_RatMX_fut(nn) = trapz(t_time,funRatkowskyL(Temp+TclimCh,ParMG_RatMX(nn,:)));

                overT(nn) = sum(Temp>ToptMGX(nn))/length(Temp)*100;
                overT_fut(nn) = sum(Temp+5>ToptMGX(nn))/length(Temp)*100;
            end
    
        else
            figure(99+nn);  set(gcf,'Position',get(0,'Screensize')); set(gcf,'color','w'); subplot(3,2,[2 4]); plot(ttL,sqrt(yyL),'o','Color',[0.6 0.6 0.6],'LineWidth',2,'MarkerSize',10); hold on;
        
            ParR_Rat(nn,:) = PsolX(2,:);
            [regX_aux,ParR_RatX(nn,1:2),~] = Regression_Tmodel(ttL,sqrt(yyL),regT,funLin,ParR_Rat(nn,:),opts);
            regR_RatX = regX_aux.^2;

            ParR_RatM(nn,:) = [PsolX(2,:) NaN NaN];
            funRatkowskyModSQ = @(x_,P) max(PsolX(2,1)*(x_+PsolX(2,2)).*(1+exp(P(1)*(x_-P(2)))),0);
            [regR_aux,ParR_RatM(nn,3:4),errM] = Regression_Tmodel(ttL,sqrt(yyL),regT,funRatkowskyModSQ,[0.1 50],opts);
            regR_RatM = regR_aux.^2;
            [regR_Arr,ParR_Arr(nn,:),errA] = Regression_Tmodel(ttL,yyL,regT,funArrhenius,[30 10],opts);
            [regR_MMRT,ParR_MMRT(nn,:),errMMRT] = Regression_Tmodel(ttL+Tr,yyL,regT+Tr,funMMRT,[2.2e4 1 8e3 330],opts);
            [regR_MMRTmod,ParR_MMRTmod(nn,:),errMMRTmod] = Regression_Tmodel(ttL+Tr,yyL,regT+Tr,funMMRTmod,[-1e5 -100 -500 -2e4 1e6 310 300],opts);
            [regR_RatMX,ParR_RatMX(nn,:),errMX] = Regression_TmodelCON(ttL,yyL,regT,funRatkowskyModL,ParR_RatM(nn,:),optsCO);

            [R_R2(nn,1),R_RMSE(nn,1),R_AIC(nn,1),R_BIC(nn,1)] = funModEval(yyL/max(yyL),funArrhenius(ttL',ParR_Arr(nn,:))/max(yyL),2);
            [R_R2(nn,2),R_RMSE(nn,2),R_AIC(nn,2),R_BIC(nn,2)] = funModEval(yyL/max(yyL),funMMRT(ttL'+Tr,ParR_MMRT(nn,:))/max(yyL),4);
            [R_R2(nn,3),R_RMSE(nn,3),R_AIC(nn,3),R_BIC(nn,3)] = funModEval(yyL/max(yyL),funLin(ttL',ParR_Rat(nn,:)).^2/max(yyL),2);
            [R_R2(nn,4),R_RMSE(nn,4),R_AIC(nn,4),R_BIC(nn,4)] = funModEval(yyL/max(yyL),funLin(ttL',ParR_RatX(nn,:)).^2/max(yyL),4);
            [R_R2(nn,5),R_RMSE(nn,5),R_AIC(nn,5),R_BIC(nn,5)] = funModEval(yyL/max(yyL),funRatkowskyModL(ttL',ParR_RatM(nn,:))/max(yyL),4);
            [R_R2(nn,6),R_RMSE(nn,6),R_AIC(nn,6),R_BIC(nn,6)] = funModEval(yyL/max(yyL),funRatkowskyModL(ttL',ParR_RatMX(nn,:))/max(yyL),4);
            [R_R2(nn,7),R_RMSE(nn,7),R_AIC(nn,7),R_BIC(nn,7)] = funModEval(yyL/max(yyL),funMMRTmod(ttL'+Tr,ParR_MMRTmod(nn,:))/max(yyL),7);
            
            [T_R2(nn,1),T_RMSE(nn,1),T_AIC(nn,1),T_BIC(nn,1)] = funModEval([yyL/max(yyL) yy_MG],[funArrhenius(ttL',ParR_Arr(nn,:))/max(yyL) val1],4);
            [T_R2(nn,2),T_RMSE(nn,2),T_AIC(nn,2),T_BIC(nn,2)] = funModEval([yyL/max(yyL) yy_MG],[funMMRT(ttL'+Tr,ParR_MMRT(nn,:))/max(yyL) val2],8);
            [T_R2(nn,3),T_RMSE(nn,3),T_AIC(nn,3),T_BIC(nn,3)] = funModEval([yyL/max(yyL) yy_MG],[funLin(ttL',ParR_Rat(nn,:)).^2/max(yyL) val3],6);
            [T_R2(nn,4),T_RMSE(nn,4),T_AIC(nn,4),T_BIC(nn,4)] = funModEval([yyL/max(yyL) yy_MG],[funLin(ttL',ParR_RatX(nn,:)).^2/max(yyL) val4],8);
            [T_R2(nn,5),T_RMSE(nn,5),T_AIC(nn,5),T_BIC(nn,5)] = funModEval([yyL/max(yyL) yy_MG],[funRatkowskyModL(ttL',ParR_RatM(nn,:))/max(yyL) val5],8);
            [T_R2(nn,6),T_RMSE(nn,6),T_AIC(nn,6),T_BIC(nn,6)] = funModEval([yyL/max(yyL) yy_MG],[funRatkowskyModL(ttL',ParR_RatMX(nn,:))/max(yyL) val6],8);
            [T_R2(nn,7),T_RMSE(nn,7),T_AIC(nn,7),T_BIC(nn,7)] = funModEval([yyL/max(yyL) yy_MG],[funMMRTmod(ttL'+Tr,ParR_MMRTmod(nn,:))/max(yyL) val7],14);
         
            RbiasArr(nn,:) = (funArrhenius(vecT',ParR_Arr(nn,:))-vecV(ii,:))/std(yyL);
            RbiasMMRT(nn,:) = (funMMRT(vecT'+Tr,ParR_MMRT(nn,:))-vecV(ii,:))/std(yyL);
            RbiasMMRTmod(nn,:) = (funMMRTmod(vecT'+Tr,ParR_MMRTmod(nn,:))-vecV(ii,:))/std(yyL);
            RbiasRat(nn,:) = (funLin(vecT',ParR_Rat(nn,:)).^2-vecV(ii,:))/std(yyL);
            RbiasRatX(nn,:) = (funLin(vecT',ParR_RatX(nn,:)).^2-vecV(ii,:))/std(yyL);
            RbiasRatM(nn,:) = (funRatkowskyModL(vecT',ParR_RatM(nn,:))-vecV(ii,:))/std(yyL);
            RbiasRatMX(nn,:) = (funRatkowskyModL(vecT',ParR_RatMX(nn,:))-vecV(ii,:))/std(yyL);
            
            funTopt = @(x_) x_ - ParR_RatM(nn,4) + 1/ParR_RatM(nn,3)*log(1+ParR_RatM(nn,3)*(x_+ParR_RatM(nn,2)));        
            [ToptR(nn),~,~,~]  = fsolve(funTopt,35,optsS);
            funTopt = @(x_) x_ - ParR_RatMX(nn,4) + 1/ParR_RatMX(nn,3)*log(1+ParR_RatMX(nn,3)*(x_+ParR_RatMX(nn,2)));        
            [ToptRX(nn),~,~,~]  = fsolve(funTopt,35,optsS);

            if runX == 1
                intD = interp1([-ParR_RatMX(nn,2) ttL' 60],[0 yyL yyL(end)*10], regT, 'pchip');
                intD(regT<=-ParR_RatMX(nn,2)) = 0;
                y_interp = smooth(regT, intD, 0.2, 'loess');
                yy_interp = y_interp;
                yy_interp(regT<=-ParR_RatMX(nn,2)) = 0;

                R_data = interp1(regT,yy_interp,Temp,'linear');
                CumR_data(nn) = trapz(t_time,R_data);
                R_data_fut = interp1(regT,yy_interp,Temp+TclimCh,'linear');
                CumR_data_fut(nn) = trapz(t_time,R_data_fut);

                CumR_Arr(nn) = trapz(t_time,funArrhenius(Temp,ParR_Arr(nn,:)));
                CumR_MMRT(nn) = trapz(t_time,funMMRT(Temp+Tr,ParR_MMRT(nn,:)));
                CumR_MMRTmod(nn) = trapz(t_time,funMMRTmod(Temp+Tr,ParR_MMRTmod(nn,:)));
                CumR_Rat(nn) = trapz(t_time,funLin(Temp,ParR_Rat(nn,:)).^2);
                CumR_RatX(nn) = trapz(t_time,funLin(Temp,ParR_RatX(nn,:)).^2);
                CumR_RatMX(nn) = trapz(t_time,funRatkowskyModL(Temp,ParR_RatMX(nn,:)));
                CumR_Arr_fut(nn) = trapz(t_time,funArrhenius(Temp+TclimCh,ParR_Arr(nn,:)));
                CumR_MMRT_fut(nn) = trapz(t_time,funMMRT(Temp+Tr+TclimCh,ParR_MMRT(nn,:)));
                CumR_MMRTmod_fut(nn) = trapz(t_time,funMMRTmod(Temp+Tr+TclimCh,ParR_MMRTmod(nn,:)));
                CumR_Rat_fut(nn) = trapz(t_time,funLin(Temp+TclimCh,ParR_Rat(nn,:)).^2);
                CumR_RatX_fut(nn) = trapz(t_time,funLin(Temp+TclimCh,ParR_RatX(nn,:)).^2);
                CumR_RatMX_fut(nn) = trapz(t_time,funRatkowskyModL(Temp+TclimCh,ParR_RatMX(nn,:)));
            end
        end 
    end
        
    figure(99+nn);
    subplot(3,2,[1 3]); plot(regT,sqrt(regMG_Arr),'Color',[0.466 0.674 0.188],'LineWidth',2,'LineStyle','-');
                        plot(regT,sqrt(regMG_MMRT),'Color',[1 0.84 0],'LineWidth',2,'LineStyle','-');
                        plot(regT,sqrt(regMG_MMRTmod),'Color',[0 0.447 0.741],'LineWidth',2,'LineStyle','-');
                        plot(regT,sqrt(regMG_Rat(1,:)),'Color',[0.93 0.66 0.55],'LineWidth',2,'LineStyle','-');  
                        plot(regT,sqrt(regMG_RatX),'Color',[0.85 0.325 0.098],'LineWidth',2,'LineStyle','-'); 
                        plot(regT,sqrt(regMG_RatM(1,:)),'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle',':');
                        plot(regT,sqrt(regMG_RatMX(1,:)),'Color','k','LineWidth',2,'LineStyle',':');
                        ylim([0 1.1*max(sqrt(vecV(1,:)))]); xlabel('Temperature [°C]'); ylabel('\surd Microbial growth [\mugC/g/h]'); set(gca,'fontsize',sizeF); 
    subplot(3,2,[2 4]); hhd(1) = plot(NaN,NaN,'o','Color',[0.6 0.6 0.6],'LineWidth',2,'MarkerSize',10);
                        hhd(2) = plot(regT,sqrt(regR_Arr),'Color',[0.466 0.674 0.188],'LineWidth',2,'LineStyle','-');
                        hhd(3) = plot(regT,sqrt(regR_MMRT),'Color',[1 0.84 0],'LineWidth',2,'LineStyle','-');
                        hhd(4) = plot(regT,sqrt(regR_MMRTmod),'Color',[0 0.447 0.741],'LineWidth',2,'LineStyle','-');
                        hhd(5) = plot(regT,sqrt(regX(2,:)),'Color',[0.93 0.66 0.55],'LineWidth',2,'LineStyle','-'); 
                        hhd(6) = plot(regT,sqrt(regR_RatX),'Color',[0.85 0.325 0.098],'LineWidth',2,'LineStyle','-'); 
                        hhd(7) = plot(regT,sqrt(regR_RatM),'Color',[0.5 0.5 0.5],'LineWidth',2,'LineStyle','-');
                        hhd(8) = plot(regT,sqrt(regR_RatMX),'Color','k','LineWidth',2,'LineStyle','-');
                        ylim([0 1.1*max(sqrt(vecV(2,:)))]);  xlabel('Temperature [°C]'); ylabel('\surd Respiration [\mugC/g/h]'); legend(hhd,{'Data','Arrhenius','MMRT','MMRT-2S','Ratkowsky-2','Ratkowsky-1','Ratkowsky DK-2','Ratkowsky DK-1'},'Orientation','horizontal'); set(gca,'fontsize',sizeF);

    if ~ismember(nn,41)
        close(99+nn) %comment this line to keep all graphs (the 70 sites)
    else
        figure(99+nn); 
            subplot(3,2,[1 3]); xlabel('Temperature [°C]'); ylabel('\surd Growth [\muC/g/h]');
                                plot([-ParMG_Rat(nn,2) -ParMG_Rat(nn,2)],[0 10],'k:');
                                plot([ParMG_Rat(nn,4) ParMG_Rat(nn,4)],[0 10],'k:');
                                plot([ToptMG(nn) ToptMG(nn)],[0 10],'k:');
            subplot(3,2,[2 4]); xlabel('Temperature [°C]'); ylabel('\surd Respiration [\muC/g/h]');
                                plot([-ParR_RatM(nn,2) -ParR_RatM(nn,2)],[0 10],'k:');
                                plot([ParR_RatM(nn,4) ParR_RatM(nn,4)],[0 10],'k:');
                                plot([ToptR(nn) ToptR(nn)],[0 10],'k:');                             
    end
end

figure(4); set(gcf,'Position',get(0,'Screensize')); set(gcf,'color','w'); set(findall(gcf,'-property','Interpreter'),'Interpreter','tex');
        subplot(2,4,1); plot([-5 50],[0 0],'k:'); hold on; ylabel('Growth bias [-]'); set(gca,'fontsize',sizeF); ylim([-0.7 2]);
                        plot(ttLL,nanmean(GbiasArr),'-','LineWidth',2,'Color',[0.466 0.674 0.188]);
                        plot(ttLL,GbiasArr,'o','Color',[0.466 0.674 0.188]);      
                        hh3(1) = plot(NaN,NaN,'-','Color',[0.466 0.674 0.188],'LineWidth',2); 
                        hh3(2) = plot(NaN,NaN,'-','Color',[1 0.84 0],'LineWidth',2); 
                        hh3(3) = plot(NaN,NaN,'-','Color',[0 0.447 0.741],'LineWidth',2);
                        hh3(4) = plot(NaN,NaN,'-','Color',[0.93 0.66 0.55],'LineWidth',2);
                        hh3(5) = plot(NaN,NaN,'-','Color',[0.85 0.325 0.098],'LineWidth',2);
                        hh3(6) = plot(NaN,NaN,'-','Color',[0.5 0.5 0.5],'LineWidth',2);
                        hh3(7) = plot(NaN,NaN,'-','Color','k','LineWidth',2); 
                        legend(hh3,{'Arrhenius','MMRT','MMRT-2S','Ratkowsky-2','Ratkowsky-1','Ratkowsky DK-2','Ratkowsky DK-1'},'Orientation','horizontal'); set(gca,'fontsize',sizeF);
        subplot(2,4,2); plot([-5 50],[0 0],'k:'); hold on; set(gca,'fontsize',sizeF); ylim([-0.7 2]);
                        plot(ttLL,nanmean(GbiasMMRT),'-','LineWidth',2,'Color',[1 0.84 0]);
                        plot(ttLL,GbiasMMRT,'o','Color',[1 0.84 0]);
                        plot(ttLL,nanmean(GbiasMMRTmod),'-','LineWidth',2,'Color',[0 0.447 0.741]);
                        plot(ttLL,GbiasMMRTmod,'o','Color',[0 0.447 0.741]);              
        subplot(2,4,3); plot([-5 50],[0 0],'k:'); hold on; set(gca,'fontsize',sizeF); ylim([-0.7 2]);
                        plot(ttLL,nanmean(GbiasRat),'-','LineWidth',2,'Color',[0.93 0.66 0.55]);
                        plot(ttLL,nanmean(GbiasRatX),'-','LineWidth',2,'Color',[0.85 0.325 0.098]);
                        plot(ttLL,GbiasRat,'o','Color',[0.93 0.66 0.55]); 
                        plot(ttLL,GbiasRatX,'+','Color',[0.85 0.325 0.098]);                            
        subplot(2,4,4); plot([-5 50],[0 0],'k:'); hold on; set(gca,'fontsize',sizeF); ylim([-0.7 2]);
                        plot(ttLL,nanmean(GbiasRat),'Color',[0.5 0.5 0.5],'LineWidth',2);
                        plot(ttLL,nanmean(GbiasRatMX),'-','LineWidth',2,'Color','k');
                        plot(ttLL,GbiasRat,'o','Color',[0.5 0.5 0.5]');
                        plot(ttLL,GbiasRatMX,'+','Color','k');
        subplot(2,4,5); plot([-5 50],[0 0],'k:'); hold on; ylabel('Respiration bias [-]'); xlabel('Temperature [°C]'); set(gca,'fontsize',sizeF); ylim([-2 0.5]);
                        plot(ttLL,nanmean(RbiasArr),'-','LineWidth',2,'Color',[0.466 0.674 0.188]);
                        plot(ttLL,RbiasArr,'o','Color',[0.466 0.674 0.188]);                     
        subplot(2,4,6); plot([-5 50],[0 0],'k:'); hold on; set(gca,'fontsize',sizeF); xlabel('Temperature [°C]'); ylim([-2 0.5]);
                        plot(ttLL,nanmean(RbiasMMRT),'-','Color',[1 0.84 0],'LineWidth',2);
                        plot(ttLL,RbiasMMRT,'o','Color',[1 0.84 0]);      
                        plot(ttLL,nanmean(RbiasMMRTmod),'-','Color',[0 0.447 0.741],'LineWidth',2);
                        plot(ttLL,RbiasMMRTmod,'o','Color',[0 0.447 0.741]);                            
        subplot(2,4,7); plot([-5 50],[0 0],'k:'); hold on; set(gca,'fontsize',sizeF); xlabel('Temperature [°C]'); ylim([-2 0.5]);
                        plot(ttLL,nanmean(RbiasRat),'-','LineWidth',2,'Color',[0.93 0.66 0.55]);
                        plot(ttLL,nanmean(RbiasRatX),'-','LineWidth',2,'Color',[0.85 0.325 0.098]);
                        plot(ttLL,RbiasRat,'o','Color',[0.93 0.66 0.55]); 
                        plot(ttLL,RbiasRatX,'+','Color',[0.85 0.325 0.098]);                            
        subplot(2,4,8); plot([-5 50],[0 0],'k:'); hold on; set(gca,'fontsize',sizeF); xlabel('Temperature [°C]'); ylim([-2 0.5]);
                        plot(ttLL,nanmean(RbiasRatM),'-','LineWidth',2,'Color',[0.5 0.5 0.5]);
                        plot(ttLL,nanmean(RbiasRatMX),'-k','LineWidth',2);
                        plot(ttLL,RbiasRatM,'o','Color',[0.5 0.5 0.5]); 
                        plot(ttLL,RbiasRatMX,'+','Color','k'); 

posZ = find(isinf(MG_AIC(:,2))); % to avoid problems related to few data points, only N>6 used
MG_AICs = MG_AIC;
MG_AICs(posZ,:) = [];
posZ = find(isinf(R_AIC(:,2)));
R_AICs = R_AIC;
R_AICs(posZ,:) = [];
posZ = find(isinf(T_AIC(:,2)));
T_AICs = T_AIC;
T_AICs(posZ,:) = [];
posZ = find(isinf(MG_BIC(:,2)));
MG_BICs = MG_BIC;
MG_BICs(posZ,:) = [];
posZ = find(isinf(R_BIC(:,2)));
R_BICs = R_BIC;
R_BICs(posZ,:) = [];
posZ = find(isinf(T_BIC(:,2)));
T_BICs = T_BIC;
T_BICs(posZ,:) = [];

funLinMAT1 = @(x_,P) P(1)*x_-P(2);
funLinMAT2 = @(x_,P) P(1)*x_+P(2);
xxMAST = linspace(0,22,10);

figure(5); set(gcf,'Position',get(0,'Screensize')); set(gcf,'color','w'); 
    subplot(2,3,1); hold on; set(gca,'fontsize',sizeF); ylabel('T_m_i_n [°C]'); xlabel('MAST [°C]'); axis([0 22 -15 0]); box on;
    subplot(2,3,3); hold on; set(gca,'fontsize',sizeF); ylabel('T_m_a_x [°C]'); xlabel('MAST [°C]'); axis([0 22 40 60]); box on;
    subplot(2,3,2); hold on; set(gca,'fontsize',sizeF); ylabel('T_t_p [°C]'); xlabel('MAST [°C]'); axis([0 22 20 45]); box on;
    subplot(2,3,4); hold on; set(gca,'fontsize',sizeF); ylabel('T_m_i_n^G [°C]'); xlabel('T_m_i_n^R [°C]'); axis([-8 -3 -18 0]); box on;
    subplot(2,3,6); hold on; set(gca,'fontsize',sizeF); ylabel('T_m_a_x^G [°C]'); xlabel('T_m_a_x^R [°C]'); axis([42 60 40 55]); box on;
    subplot(2,3,5); hold on; set(gca,'fontsize',sizeF); ylabel('T_t_p^G [°C]'); xlabel('T_t_p^R [°C]'); axis([27 41 25 42]); box on;

for nn = 1:nSites
   subplot(2,3,1); plot(MAST(nn),-ParMG_RatM(nn,2),'s','Color',[0.494 0.184 0.556]);
                   plot(MAST(nn),-ParR_RatM(nn,2),'d','Color',[0.369 0.529 0.161]);
   subplot(2,3,3); plot(MAST(nn),ParMG_RatM(nn,4),'s','Color',[0.494 0.184 0.556]);
                   plot(MAST(nn),ParR_RatM(nn,4),'d','Color',[0.369 0.529 0.161]);
   subplot(2,3,2); plot(MAST(nn),ToptMG(nn),'s','Color',[0.494 0.184 0.556]);
                   plot(MAST(nn),ToptR(nn),'d','Color',[0.369 0.529 0.161]);
   
   subplot(2,3,4); plot(-ParR_RatM(nn,2),-ParMG_RatM(nn,2),'ok');
   subplot(2,3,6); plot(ParR_RatM(nn,4),ParMG_RatM(nn,4),'ok');
   subplot(2,3,5); plot(ToptR(nn),ToptMG(nn),'ok');
end
[matTmin(1,:),Pmin(1,:),~] = Regression_Tmodel(MAST,-ParMG_RatM(:,2),xxMAST,funLinMAT1,[1 -10],opts);
[matTmin(2,:),Pmin(2,:),~] = Regression_Tmodel(MAST,-ParR_RatM(:,2),xxMAST,funLinMAT1,[1 -10],opts);
[matTmax(1,:),Pmax(1,:),~] = Regression_Tmodel(MAST,ParMG_RatM(:,4),xxMAST,funLinMAT1,[1 40],opts);
[matTmax(2,:),Pmax(2,:),~] = Regression_Tmodel(MAST(ParR_RatM(:,4)<65),ParR_RatM(ParR_RatM(:,4)<65,4),xxMAST,funLinMAT1,[1 40],opts);
[matTopt(1,:),Popt(1,:),~] = Regression_Tmodel(MAST,ToptMG,xxMAST,funLinMAT1,[1 40],opts);
[matTopt(2,:),Popt(2,:),~] = Regression_Tmodel(MAST(ToptR<45),ToptR(ToptR<45),xxMAST,funLinMAT1,[1 40],opts);
subplot(2,3,1); hhx(1) = plot(xxMAST,matTmin(1,:),'-','Color',[0.494 0.184 0.556],'LineWidth',2); text(7,-13,['T_m_i_n = ',num2str(sprintf('%.2f',Pmin(1,1))),'*MAST - ',num2str(sprintf('%.2f',Pmin(1,2)))],'Color',[0.494 0.184 0.556],'FontSize',15);
                hhx(2) = plot(xxMAST,matTmin(2,:),'-','Color',[0.369 0.529 0.161],'LineWidth',2); text(2,-3,['T_m_i_n = ',num2str(sprintf('%.2f',Pmin(2,1))),'*MAST - ',num2str(sprintf('%.2f',Pmin(2,2)))],'Color',[0.369 0.529 0.161],'FontSize',15);
                legend(hhx,{'Growth','Respiration'},'Orientation','horizontal'); set(gca,'fontsize',sizeF);
subplot(2,3,3); plot(xxMAST,matTmax(1,:),'-','Color',[0.494 0.184 0.556],'LineWidth',2); text(7,43,['T_m_a_x = ',num2str(sprintf('%.2f',Pmax(1,1))),'*MAST + ',num2str(sprintf('%.2f',-Pmax(1,2)))],'Color',[0.494 0.184 0.556],'FontSize',15);
                plot(xxMAST,matTmax(2,:),'-','Color',[0.369 0.529 0.161],'LineWidth',2); text(2,55,['T_m_a_x = ',num2str(sprintf('%.2f',Pmax(2,1))),'*MAST + ',num2str(sprintf('%.2f',-Pmax(2,2)))],'Color',[0.369 0.529 0.161],'FontSize',15);
subplot(2,3,2); plot(xxMAST,matTopt(1,:),'-','Color',[0.494 0.184 0.556],'LineWidth',2); text(7,23,['T_t_p = ',num2str(sprintf('%.2f',Popt(1,1))),'*MAST + ',num2str(sprintf('%.2f',-Popt(1,2)))],'Color',[0.494 0.184 0.556],'FontSize',15);
                plot(xxMAST,matTopt(2,:),'-','Color',[0.369 0.529 0.161],'LineWidth',2); text(2,40,['T_t_p = ',num2str(sprintf('%.2f',Popt(2,1))),'*MAST + ',num2str(sprintf('%.2f',-Popt(2,2)))],'Color',[0.369 0.529 0.161],'FontSize',15);

[matTminGR,PminGR,~] = Regression_Tmodel(-ParR_RatM(:,2),-ParMG_RatM(:,2),[-10 0],funLinMAT2,[1 0],opts);
[matTmaxGR,PmaxGR,~] = Regression_Tmodel(ParR_RatM(:,4),ParMG_RatM(:,4),[40 60],funLinMAT2,[0 50],opts);
[matToptGR,PoptGR,~] = Regression_Tmodel(ToptR,ToptMG,[20 45],funLinMAT2,[0 30],opts);
subplot(2,3,4); plot([-10 0],matTminGR,'-k','LineWidth',2); text(-8,-2,['T_m_i_n^G = ',num2str(sprintf('%.2f',PminGR(1))),'*T_m_i_n^R - ',num2str(sprintf('%.2f',-PminGR(2)))],'Color','k','FontSize',15);
subplot(2,3,6); plot([40 60],matTmaxGR,'-k','LineWidth',2); text(45,52,['T_m_a_x^G = ',num2str(sprintf('%.2f',PmaxGR(1))),'*T_m_a_x^R + ',num2str(sprintf('%.2f',PmaxGR(2)))],'Color','k','FontSize',15);
subplot(2,3,5); plot([20 45],matToptGR,'-k','LineWidth',2); text(28,40,['T_t_p^G = ',num2str(sprintf('%.2f',PoptGR(1))),'*T_t_p^R + ',num2str(sprintf('%.2f',PoptGR(2)))],'Color','k','FontSize',15);

cf_data = CumR_data./CumMG_data;
flu_Arr = -(cf_data.*CumMG_Arr - CumR_Arr)./CumR_data*100;
flu_MMRT = -(cf_data.*CumMG_MMRT - CumR_MMRT)./CumR_data*100;
flu_MMRTmod = -(cf_data.*CumMG_MMRTmod - CumR_MMRTmod)./CumR_data*100;
flu_Rat = -(cf_data.*CumMG_RatX - CumR_RatX)./CumR_data*100;
flu_RatMX = -(cf_data.*CumMG_RatMX - CumR_RatMX)./CumR_data*100;
cf_data_fut = CumR_data_fut./CumMG_data_fut;
flu_Arr_fut = -(cf_data_fut.*CumMG_Arr_fut - CumR_Arr_fut)./CumR_data_fut*100;
flu_MMRT_fut = -(cf_data_fut.*CumMG_MMRT_fut - CumR_MMRT_fut)./CumR_data_fut*100;
flu_MMRTmod_fut = -(cf_data_fut.*CumMG_MMRTmod_fut - CumR_MMRTmod_fut)./CumR_data_fut*100;
flu_Rat_fut = -(cf_data_fut.*CumMG_RatX_fut - CumR_RatX_fut)./CumR_data_fut*100;
flu_RatMX_fut = -(cf_data_fut.*CumMG_RatMX_fut - CumR_RatMX_fut)./CumR_data_fut*100;

funEmi = @(x_,P) P(1).*x_-P(2);
[~,pEmi_Arr,~] = Regression_Tmodel(MAST,flu_Arr,xxMAST,funEmi,[-10 -100],opts);
[~,pEmi_MMRT,~] = Regression_Tmodel(MAST,flu_MMRT,xxMAST,funEmi,[-10 -100],opts);
[~,pEmi_MMRTmod,~] = Regression_Tmodel(MAST,flu_MMRTmod,xxMAST,funEmi,[-10 -100],opts);
[~,pEmi_Rat,~] = Regression_Tmodel(MAST,flu_Rat,xxMAST,funEmi,[-10 -100],opts);
[~,pEmi_RatMX,~] = Regression_Tmodel(MAST,flu_RatMX,xxMAST,funEmi,[-10 -100],opts);
[~,pEmi_Arr_fut,~] = Regression_Tmodel(MAST,flu_Arr_fut,xxMAST,funEmi,[-10 -100],opts);
[~,pEmi_MMRT_fut,~] = Regression_Tmodel(MAST,flu_MMRT_fut,xxMAST,funEmi,[-10 -100],opts);
[~,pEmi_MMRTmod_fut,~] = Regression_Tmodel(MAST,flu_MMRTmod_fut,xxMAST,funEmi,[-10 -100],opts);
[~,pEmi_Rat_fut,~] = Regression_Tmodel(MAST,flu_Rat_fut,xxMAST,funEmi,[-10 -100],opts);
[~,pEmi_RatMX_fut,~] = Regression_Tmodel(MAST,flu_RatMX_fut,xxMAST,funEmi,[-10 -100],opts);
 
figure(8); set(gcf,'Position',get(0,'Screensize')); set(gcf,'color','w');
    subplot(2,1,1); 
        hha(1) = plot(MAST,flu_Arr,'o','Color',[0.466 0.674 0.188]); hold on;
        plot(xxMAST,funEmi(xxMAST,pEmi_Arr),'Color',[0.466 0.674 0.188],'Linewidth',2);
        hha(2) = plot(MAST,flu_MMRT,'o','Color',[1 0.84 0]); 
        plot(xxMAST,funEmi(xxMAST,pEmi_MMRT),'Color',[1 0.84 0],'Linewidth',2);
        hha(3) = plot(MAST,flu_MMRTmod,'o','Color',[0 0.447 0.741]); 
        plot(xxMAST,funEmi(xxMAST,pEmi_MMRTmod),'Color',[0 0.447 0.741],'Linewidth',2);
        hha(4) = plot(MAST,flu_Rat,'o','Color',[0.85 0.325 0.098]);
        plot(xxMAST,funEmi(xxMAST,pEmi_Rat),'Color',[0.85 0.325 0.098],'Linewidth',2);
        hha(5) = plot(MAST,flu_RatMX,'o','Color','k');
        plot(xxMAST,funEmi(xxMAST,pEmi_RatMX),'Color','k','Linewidth',2); 
        legend(hha,{'Arrhenius','MMRT','MMRT-2S','Ratkowsky-1','Ratkowsky DK-1'},'Orientation','horizontal'); axis([0 22 -150 150]); xlabel('MAST [°C]'); ylabel({'Dev [%]'}); title('Current scenario'); set(gca,'fontsize',sizeF); 
    subplot(2,1,2); 
        plot(MAST+TclimCh,flu_Arr_fut,'o','Color',[0.466 0.674 0.188]); hold on;
        plot(xxMAST+TclimCh,funEmi(xxMAST,pEmi_Arr_fut),'Color',[0.466 0.674 0.188],'Linewidth',2);
        plot(MAST+TclimCh,flu_MMRT_fut,'o','Color',[1 0.84 0]); 
        plot(xxMAST+TclimCh,funEmi(xxMAST,pEmi_MMRT_fut),'Color',[1 0.84 0],'Linewidth',2);
        plot(MAST+TclimCh,flu_MMRTmod_fut,'o','Color',[0 0.447 0.741]); 
        plot(xxMAST+TclimCh,funEmi(xxMAST,pEmi_MMRTmod_fut),'Color',[0 0.447 0.741],'Linewidth',2);
        plot(MAST+TclimCh,flu_Rat_fut,'o','Color',[0.85 0.325 0.098]);
        plot(xxMAST+TclimCh,funEmi(xxMAST,pEmi_Rat_fut),'Color',[0.85 0.325 0.098],'Linewidth',2); 
        plot(MAST+TclimCh,flu_RatMX_fut,'o','Color','k');
        plot(xxMAST+TclimCh,funEmi(xxMAST,pEmi_RatMX_fut),'Color','k','Linewidth',2); 
        xlabel('MAST [°C]'); ylabel({'Dev [%]'}); title('Future scenario +5.8°C'); axis([0+TclimCh 22+TclimCh -150 150]); set(gca,'fontsize',sizeF); 

predF_Arr = funEmi(MAST_Europe,pEmi_Arr);
predF_MMRT = funEmi(MAST_Europe,pEmi_MMRT);
predF_MMRTmod = funEmi(MAST_Europe,pEmi_MMRTmod);
predF_Rat = funEmi(MAST_Europe,pEmi_Rat);
predF_RatM = funEmi(MAST_Europe,pEmi_RatMX);
predF_Arr_fut = funEmi(MAST_Europe,pEmi_Arr_fut);
predF_MMRT_fut = funEmi(MAST_Europe,pEmi_MMRT_fut);
predF_MMRTmod_fut = funEmi(MAST_Europe,pEmi_MMRTmod_fut);
predF_Rat_fut = funEmi(MAST_Europe,pEmi_Rat_fut);
predF_RatM_fut = funEmi(MAST_Europe,pEmi_RatMX_fut);

meanpdF_Arr = nanmean(abs(flu_Arr(:)));
meanpdF_MMRT = nanmean(abs(flu_MMRT(:)));
meanpdF_MMRTmod = nanmean(abs(flu_MMRTmod(:)));
meanpdF_Rat = nanmean(abs(flu_Rat(:)));
meanpdF_RatM = nanmean(abs(flu_RatMX(:)));
meanpdF_Arr_fut = nanmean(abs(flu_Arr_fut(:)));
meanpdF_MMRT_fut = nanmean(abs(flu_MMRT_fut(:)));
meanpdF_MMRTmod_fut = nanmean(abs(flu_MMRTmod_fut(:)));
meanpdF_Rat_fut = nanmean(abs(flu_Rat_fut(:)));
meanpdF_RatM_fut = nanmean(abs(flu_RatMX_fut(:)));

nancolor = [1,1,1]; % RGB color for NaN (white)
numSteps = 100; % Half of the colormap will be for negative, half for positive
blueGradient = [linspace(1,0,numSteps)',linspace(1,0,numSteps)',linspace(1,1,numSteps)'];
redGradient = [linspace(1,1,numSteps)',linspace(1,0,numSteps)',linspace(1,0,numSteps)'];
customCmap = [flipud(blueGradient); ones(12,3); redGradient];
customCmap = [nancolor; customCmap];
maxV = 100;
minV = -maxV;

coast = shaperead('landareas.shp','UseGeoCoords',true);  % installed with MATLAB
figure(6); set(gcf,'Position',get(0,'Screensize')); set(gcf,'color','w');
subplot(2,5,1); axesm('MapProjection', 'mercator');
    predF_Arr_nan = predF_Arr;
    predF_Arr_nan(isnan(predF_Arr)) = -999;
    hgA = geoshow(predF_Arr,R_Europe,'DisplayType','texturemap');
    colormap(customCmap);
    caxis([minV,maxV]);
    set(hgA,'CData',predF_Arr_nan);
    title('Arrhenius');
    colorbar;
    grid on;
    axis tight;
    hold on;
    axLimits = axis;
    geoshow([coast.Lat],[coast.Lon],'DisplayType','polygon','FaceColor','none','EdgeColor','k','LineWidth',0.5);
    axis(axLimits);
    for nn = 1:nSites
        if flu_Arr(nn) > maxV
            markerColor_Arr = [1 0 0];
        elseif flu_Arr(nn) < minV
            markerColor_Arr = [0 0 1];
        else
            markerColor_Arr = interp1(linspace(minV,maxV,size(customCmap,1)),customCmap,flu_Arr(nn));
        end
        geoshow(lat(nn),lon(nn),'DisplayType','Point','Marker','o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',markerColor_Arr);
    end
    text(-0.2,1.7,sprintf('%.1f%%',meanpdF_Arr),'FontSize',25);
subplot(2,5,4); axesm('MapProjection', 'mercator');
    predF_MMRTmod_nan = predF_MMRTmod;
    predF_MMRTmod_nan(isnan(predF_MMRTmod)) = -999;
    hgA = geoshow(predF_MMRTmod,R_Europe,'DisplayType','texturemap');
    colormap(customCmap);
    caxis([minV,maxV]);
    set(hgA,'CData',predF_MMRTmod_nan);
    title('MMRT-2S');
    colorbar;
    grid on;
    axis tight;
    hold on;
    axLimits = axis;
    geoshow([coast.Lat],[coast.Lon],'DisplayType','polygon','FaceColor','none','EdgeColor','k','LineWidth',0.5);
    axis(axLimits);
    for nn = 1:nSites
        if flu_MMRTmod > maxV
            markerColor_MMRT = [1 0 0];
        elseif flu_MMRTmod(nn) < minV
            markerColor_MMRT = [0 0 1];
        else
            markerColor_MMRT = interp1(linspace(minV,maxV,size(customCmap,1)),customCmap,flu_MMRTmod(nn));
        end
        geoshow(lat(nn),lon(nn),'DisplayType','Point','Marker','o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',markerColor_MMRT);
    end
    text(-0.2,1.7,sprintf('%.1f%%',meanpdF_MMRTmod),'FontSize',25);
subplot(2,5,3); axesm('MapProjection', 'mercator');
    predF_MMRT_nan = predF_MMRT;
    predF_MMRT_nan(isnan(predF_MMRT)) = -999;
    hgA = geoshow(predF_MMRT,R_Europe,'DisplayType','texturemap');
    colormap(customCmap);
    caxis([minV,maxV]);
    set(hgA,'CData',predF_MMRT_nan);
    title('MMRT');
    colorbar;
    grid on;
    axis tight;
    hold on;
    axLimits = axis;
    geoshow([coast.Lat],[coast.Lon],'DisplayType','polygon','FaceColor','none','EdgeColor','k','LineWidth',0.5);
    axis(axLimits);
    for nn = 1:nSites
        if flu_MMRT(nn) > maxV
            markerColor_MMRT = [1 0 0];
        elseif flu_MMRT(nn) < minV
            markerColor_MMRT = [0 0 1];
        else
            markerColor_MMRT = interp1(linspace(minV,maxV,size(customCmap,1)),customCmap,flu_MMRT(nn));
        end
        geoshow(lat(nn),lon(nn),'DisplayType','Point','Marker','o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',markerColor_MMRT);
    end
    text(-0.2,1.7,sprintf('%.1f%%',meanpdF_MMRT),'FontSize',25);
subplot(2,5,2); axesm('MapProjection', 'mercator');
    predF_Rat_nan = predF_Rat;
    predF_Rat_nan(isnan(predF_Rat)) = -999;
    hgA = geoshow(predF_Rat,R_Europe,'DisplayType','texturemap');
    colormap(customCmap);
    caxis([minV,maxV]);
    set(hgA,'CData',predF_Rat_nan);
    title('Ratkowsky');
    colorbar;
    grid on;
    axis tight;
    hold on;
    axLimits = axis;
    geoshow([coast.Lat],[coast.Lon],'DisplayType','polygon','FaceColor','none','EdgeColor','k','LineWidth',0.5);
    axis(axLimits);
    for nn = 1:nSites
        if flu_Rat(nn) > maxV
            markerColor_Rat = [1 0 0];
        elseif flu_Rat(nn) < minV
            markerColor_Rat = [0 0 1];
        else
            markerColor_Rat = interp1(linspace(minV,maxV,size(customCmap,1)),customCmap,flu_Rat(nn));
        end
        geoshow(lat(nn),lon(nn),'DisplayType','Point','Marker','o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',markerColor_Rat);
    end
    text(-0.2,1.7,sprintf('%.1f%%',meanpdF_Rat),'FontSize',25);
subplot(2,5,5); axesm('MapProjection', 'mercator');
    predF_RatM_nan = predF_RatM;
    predF_RatM_nan(isnan(predF_RatM)) = -999;
    hgA = geoshow(predF_RatM,R_Europe,'DisplayType','texturemap');
    colormap(customCmap);
    caxis([minV,maxV]);
    set(hgA,'CData',predF_RatM_nan);
    title('Ratkowsky DK');
    colorbar;
    grid on;
    axis tight;
    hold on;
    axLimits = axis;
    geoshow([coast.Lat],[coast.Lon],'DisplayType','polygon','FaceColor','none','EdgeColor','k','LineWidth',0.5);
    axis(axLimits);
    for nn = 1:nSites
        if flu_RatMX(nn) > maxV
            markerColor_RatM = [1 0 0];
        elseif flu_RatMX(nn) < minV
            markerColor_RatM = [0 0 1];
        else
            markerColor_RatM = interp1(linspace(minV,maxV,size(customCmap,1)),customCmap,flu_RatMX(nn));
        end
        geoshow(lat(nn),lon(nn),'DisplayType','Point','Marker','o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',markerColor_RatM);
    end
    text(-0.2,1.7,sprintf('%.1f%%',meanpdF_RatM),'FontSize',25);
subplot(2,5,6); axesm('MapProjection', 'mercator');
    predF_Arr_nan_fut = predF_Arr_fut;
    predF_Arr_nan_fut(isnan(predF_Arr_fut)) = -999;
    hgA = geoshow(predF_Arr_fut,R_Europe,'DisplayType','texturemap');
    colormap(customCmap);
    caxis([minV,maxV]);
    set(hgA,'CData',predF_Arr_nan_fut);
    title('Arrhenius');
    colorbar;
    grid on;
    axis tight;
    hold on;
    axLimits = axis;
    geoshow([coast.Lat],[coast.Lon],'DisplayType','polygon','FaceColor','none','EdgeColor','k','LineWidth',0.5);
    axis(axLimits);
    for nn = 1:nSites
        if flu_Arr_fut(nn) > maxV
            markerColor_Arr_fut = [1 0 0];
        elseif flu_Arr_fut(nn) < minV
            markerColor_Arr_fut = [0 0 1];
        else
            markerColor_Arr_fut = interp1(linspace(minV,maxV,size(customCmap,1)),customCmap,flu_Arr_fut(nn));
        end
        geoshow(lat(nn),lon(nn),'DisplayType','Point','Marker','o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',markerColor_Arr_fut);
    end
    text(-0.2,1.7,sprintf('%.1f%%',meanpdF_Arr_fut),'FontSize',25);
subplot(2,5,8); axesm('MapProjection', 'mercator');
    predF_MMRT_nan_fut = predF_MMRT_fut;
    predF_MMRT_nan_fut(isnan(predF_MMRT_fut)) = -999;
    hgA = geoshow(predF_MMRT_fut,R_Europe,'DisplayType','texturemap');
    colormap(customCmap);
    caxis([minV,maxV]);
    set(hgA,'CData',predF_MMRT_nan_fut);
    title('MMRT');
    colorbar;
    grid on;
    axis tight;
    hold on;
    axLimits = axis;
    geoshow([coast.Lat],[coast.Lon],'DisplayType','polygon','FaceColor','none','EdgeColor','k','LineWidth',0.5);
    axis(axLimits);
    for nn = 1:nSites
        if flu_MMRT_fut(nn) > maxV
            markerColor_MMRT_fut = [1 0 0];
        elseif flu_MMRT_fut(nn) < minV
            markerColor_MMRT_fut = [0 0 1];
        else
            markerColor_MMRT_fut = interp1(linspace(minV,maxV,size(customCmap,1)),customCmap,flu_MMRT_fut(nn));
        end
        geoshow(lat(nn),lon(nn),'DisplayType','Point','Marker','o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',markerColor_MMRT_fut);
    end
    text(-0.2,1.7,sprintf('%.1f%%',meanpdF_MMRT_fut),'FontSize',25);
subplot(2,5,9); axesm('MapProjection', 'mercator');
    predF_MMRTmod_nan_fut = predF_MMRTmod_fut;
    predF_MMRTmod_nan_fut(isnan(predF_MMRTmod_fut)) = -999;
    hgA = geoshow(predF_MMRTmod_fut,R_Europe,'DisplayType','texturemap');
    colormap(customCmap);
    caxis([minV,maxV]);
    set(hgA,'CData',predF_MMRTmod_nan_fut);
    title('MMRT-2S');
    colorbar;
    grid on;
    axis tight;
    hold on;
    axLimits = axis;
    geoshow([coast.Lat],[coast.Lon],'DisplayType','polygon','FaceColor','none','EdgeColor','k','LineWidth',0.5);
    axis(axLimits);
    for nn = 1:nSites
        if flu_MMRTmod_fut(nn) > maxV
            markerColor_MMRT_fut = [1 0 0];
        elseif flu_MMRTmod_fut(nn) < minV
            markerColor_MMRT_fut = [0 0 1];
        else
            markerColor_MMRT_fut = interp1(linspace(minV,maxV,size(customCmap,1)),customCmap,flu_MMRTmod_fut(nn));
        end
        geoshow(lat(nn),lon(nn),'DisplayType','Point','Marker','o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',markerColor_MMRT_fut);
    end
    text(-0.2,1.7,sprintf('%.1f%%',meanpdF_MMRTmod_fut),'FontSize',25);
subplot(2,5,7); axesm('MapProjection', 'mercator');
    predF_Rat_nan_fut = predF_Rat_fut;
    predF_Rat_nan_fut(isnan(predF_Rat_fut)) = -999;
    hgA = geoshow(predF_Rat_fut,R_Europe,'DisplayType','texturemap');
    colormap(customCmap);
    caxis([minV,maxV]);
    set(hgA,'CData',predF_Rat_nan_fut);
    title('Ratkowsky');
    colorbar;
    grid on;
    axis tight;
    hold on;
    axLimits = axis;
    geoshow([coast.Lat],[coast.Lon],'DisplayType','polygon','FaceColor','none','EdgeColor','k','LineWidth',0.5);
    axis(axLimits);
    for nn = 1:nSites
        if flu_Rat_fut(nn) > maxV
            markerColor_Rat_fut = [1 0 0];
        elseif flu_Rat_fut(nn) < minV
            markerColor_Rat_fut = [0 0 1];
        else
            markerColor_Rat_fut = interp1(linspace(minV,maxV,size(customCmap,1)),customCmap,flu_Rat_fut(nn));
        end
        geoshow(lat(nn),lon(nn),'DisplayType','Point','Marker','o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',markerColor_Rat_fut);
    end
    text(-0.2,1.7,sprintf('%.1f%%',meanpdF_Rat_fut),'FontSize',25);
subplot(2,5,10); axesm('MapProjection', 'mercator');
    predF_RatM_nan_fut = predF_RatM_fut;
    predF_RatM_nan_fut(isnan(predF_RatM_fut)) = -999;
    hgA = geoshow(predF_RatM_fut,R_Europe,'DisplayType','texturemap');
    colormap(customCmap);
    caxis([minV,maxV]);
    set(hgA,'CData',predF_RatM_nan_fut);
    title('Ratkowsky DK');
    colorbar;
    grid on;
    axis tight;
    hold on;
    axLimits = axis;
    geoshow([coast.Lat],[coast.Lon],'DisplayType','polygon','FaceColor','none','EdgeColor','k','LineWidth',0.5);
    axis(axLimits);
    for nn = 1:nSites
        if flu_RatMX_fut(nn) > maxV
            markerColor_RatMX_fut = [1 0 0];
        elseif flu_RatMX_fut(nn) < minV
            markerColor_RatMX_fut = [0 0 1];
        else
            markerColor_RatMX_fut = interp1(linspace(minV,maxV,size(customCmap,1)),customCmap,flu_RatMX_fut(nn));
        end
        geoshow(lat(nn),lon(nn),'DisplayType','Point','Marker','o','MarkerSize',10,'MarkerEdgeColor','k','MarkerFaceColor',markerColor_RatMX_fut);
    end
    text(-0.2,1.7,sprintf('%.1f%%',meanpdF_RatM_fut),'FontSize',25);

 beep

function [regY,P,err] = Regression_TmodelCON(Tp,Yp,regT,eq,P0,opts)

    errfh = @(P,x,z) sum((z(:)-eq(x(:),P)).^2);
    [P,err,a,b,c,d,hessian] = fmincon(errfh,P0,[],[],[],[],zeros(1,length(P0)),[],[],opts,Tp,Yp);
    regY = eq(regT,P);

end

function [regY,P,err] = Regression_Tmodel(Tp,Yp,regT,eq,P0,opts)

    errfh = @(P,x,z) sum((z(:)-eq(x(:),P)).^2);
    [P,err,a,b] = fminsearch(errfh,P0,opts,Tp,Yp);
    regY = eq(regT,P);

end

function [R2,RMSE,AICc,BIC] = funModEval(Yp,Ymod,nPar)
    
    RSS = sum((Yp-Ymod).^2);
    R2 = 1 - RSS/sum((Yp-mean(Yp)).^2);
    nData = length(Yp);
    RMSE = sqrt(RSS/nData);
    AIC = nData*log(RSS/nData)+2*(nPar+1);
    AICc = AIC + 2*nPar*(nPar+1)/(nData-nPar-1);
    if nData-nPar-1<=0
        pa = -1/(nData-nPar-2);
        AICc = AIC + 2*nPar*(nPar+1)/pa;
    end
    BIC = nData*log(RSS/nData)+nPar*log(nData);
    
end