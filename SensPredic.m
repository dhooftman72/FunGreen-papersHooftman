function Sensi_predic = SensPredic(xDOrg,ToTestyDMarg,Sense_factor,trait,Beta,yPredic_SensBase)
xSens = xDOrg;
basemode = 1;
if trait == [3,4,13:16,21,22] % so a categorical trait
    basemode = 0;
end
Sensi_predic = zeros(1,5);
% Distance (X)
factor = 1;
chance = ((10.^xSens(:,1))-1).*(1-Sense_factor);
xSens(:,1) =  log10(chance+1);
Sensi_predic = SensPredicFunc(Sensi_predic,factor,ToTestyDMarg,Sense_factor,Beta,xSens,basemode,yPredic_SensBase);
xSens = xDOrg; % set varables back

%TotalGI (Tot)
factor = 2;
chance  = ((10.^xSens(:,factor))-1).*(1+Sense_factor);
xSens(:,factor) =  log10(chance+1);
Sensi_predic = SensPredicFunc(Sensi_predic,factor,ToTestyDMarg,Sense_factor,Beta,xSens,basemode,yPredic_SensBase);
xSens = xDOrg; % set varables back

%RoadGI (Road)
factor = 3;
xSens = AddGI(xSens,factor,Sense_factor);
Sensi_predic = SensPredicFunc(Sensi_predic,factor,ToTestyDMarg,Sense_factor,Beta,xSens,basemode,yPredic_SensBase);
xSens = xDOrg; % set varables back
%
%Other Linear GI (Other)
factor = 4;
xSens = AddGI(xSens,factor,Sense_factor);
Sensi_predic = SensPredicFunc(Sensi_predic,factor,ToTestyDMarg,Sense_factor,Beta,xSens,basemode,yPredic_SensBase);
xSens = xDOrg; % set varables back

%AreaGI (Area)
factor = 5;
xSens = AddGI(xSens,factor,Sense_factor);
Sensi_predic = SensPredicFunc(Sensi_predic,factor,ToTestyDMarg,Sense_factor,Beta,xSens,basemode,yPredic_SensBase);
xSens = xDOrg; % set varables back
end % Funtion SensiPredic

function Sensi_predic = SensPredicFunc(Sensi_predic,factor,ToTestyDMarg,Sense_factor,Beta,xCalc,basemode,yPredic_SensBase)
for i = 1:1:length(ToTestyDMarg)
    yPredict(i) = (Beta(1).*xCalc(i,1)) + (Beta(2).*xCalc(i,2)) +...
        (Beta(3).*xCalc(i,3)) + (Beta(4).*xCalc(i,4)) + (Beta(5).*xCalc(i,5)) + Beta(6);
end
yPredic = nanmean(reshape(yPredict,[],1)); % so this is 1 value
Sensi_predic(factor) = (-(((abs(basemode-yPredic))./(abs(basemode-yPredic_SensBase)))-1))./Sense_factor;
end % SensiPredicFunc

function xSens = AddGI(xSens,factor,Sense_factor)
Total_GIClean = ((10.^xSens(:,2))-1) - ((10.^xSens(:,factor))-1);
chance  = ((10.^xSens(:,factor))-1).*(1+Sense_factor);
xSens(:,factor) =  log10(chance+1);
xSens(:,2) = log10((Total_GIClean +chance)+1);% make sure aso total GI is elevated
end