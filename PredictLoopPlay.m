function PredictLoop(loop,List,PropRadi,trait,traittype)
rng('shuffle')
Landscape_toGet_1 = randperm(12,6);
Landscape_toGet_2 = randperm(12,6)+12;
Landscape_toGet_3 = randperm(12,6)+24;
Training_landscapes = sort([Landscape_toGet_1,Landscape_toGet_2,Landscape_toGet_3]);
Testing_landscapes = [1:36];
Testing_landscapes(Training_landscapes) = [];

Numbers_training = [];
for i = 1:length(Training_landscapes)
    number = find(PropRadi.SpeciesSegmentIDs.Landscape == Training_landscapes(i));
    number =  reshape(number,1,[]);
    Numbers_training = [Numbers_training, number]; %#ok<AGROW>
end
Numbers_testing = [];
for i = 1:length(Testing_landscapes)
    number = find(PropRadi.SpeciesSegmentIDs.Landscape == Testing_landscapes(i));
    number =  reshape(number,1,[]);
    Numbers_testing = [Numbers_testing, number]; %#ok<AGROW>;
end


%% Training model
Ins = MakeIns(PropRadi,Numbers_training,List,trait,traittype);
Outs = RegressFunc(Ins);% Here are the Marginal values calculated
xD = [log10(Outs.XC+1),log10(Outs.total_GIC +1),log10(Outs.LinRoadC+1),log10(Outs.LinOtherC+1),log10(Outs.AreaC+1)];
xD =[xD,Outs.CountryC];
MaX = max(Outs.yD);
%% Do the multiple variable fitting
clear prior
% prior = [Outs.Disrico,Outs.Totalrico,Outs.LinRoadrico,Outs.LinOtherrico,Outs.Arearico,Outs.constant,0];
% series = [1 2 3 4 5 6];
% seriesRand = [6];
% modelFun = @(b,xD,MaX) (b(1).*xD(:,1)) + (b(2).*xD(:,2)) + (b(3).*xD(:,3)) + (b(4).*xD(:,4)) + (b(5).*xD(:,5))+ b(6);
% modelFunw = @(b,xD)modelFun(b,xD,MaX);
% Options = statset('FunValCheck','off','Display','off','MaxIter',100,...
%     'TolFun',1.0000e-4, 'TolX',1.0e-4, 'Jacobian','off', 'DerivStep', 6.0555e-05, 'OutputFcn',' ');
% Group = grp2idx(Outs.CountryC);%[1:length(yD)];
%     [bFit,~,stats] = nlmefit(xD,Outs.yD,Group,[],modelFunw,prior,'RefineBeta0','off','Options',Options,'FEParamsSelect',series,'REParamsSelect',seriesRand); %#ok<*NASGU>

prior = [Outs.Disrico,Outs.Totalrico,Outs.constant];
series = [1 2 3];
seriesRand = [3];
modelFun = @(b,xD,MaX) (b(1).*xD(:,1)) + (b(2).*xD(:,2)) + b(3);
modelFunw = @(b,xD)modelFun(b,xD,MaX);
Options = statset('FunValCheck','off','Display','off','MaxIter',100,...
    'TolFun',1.0000e-4, 'TolX',1.0e-4, 'Jacobian','off', 'DerivStep', 6.0555e-05, 'OutputFcn',' ');
Group = grp2idx(Outs.CountryC);%[1:length(yD)];
    [bFitOne,~,statsOne] = nlmefit(xD,Outs.yD,Group,[],modelFunw,prior,'RefineBeta0','off','Options',Options,'FEParamsSelect',series,'REParamsSelect',seriesRand); %#ok<*NASGU>    
    
prior = [Outs.LinRoadrico,Outs.LinOtherrico,Outs.Arearico];
series = [1 2 3];
seriesRand = [3];
modelFun = @(bdue,xD,MaX) (bdue(1).*xD(:,3)) + (bdue(2).*xD(:,4)) + (bdue(3).*xD(:,5));
modelFunw = @(bdue,xD)modelFun(bdue,xD,MaX);
Options = statset('FunValCheck','off','Display','off','MaxIter',100,...
     'TolFun',1.0000e-4, 'TolX',1.0e-4, 'Jacobian','off', 'DerivStep', 6.0555e-05, 'OutputFcn',' ');
 Group = grp2idx(Outs.CountryC);%[1:length(yD)];
     [bFitdue,~,stats] = nlmefit(xD,statsOne.pres,Group,[],modelFunw,prior,'RefineBeta0','off','Options',Options,'FEParamsSelect',series,'REParamsSelect',seriesRand); %#ok<*NASGU>   

series = [1 2 3 4 5 6];
bFit = [bFitOne(1:2);bFitdue;bFitOne(3)];
Stats.sebeta = [statsOne.sebeta(1:2),stats.sebeta,statsOne.sebeta(3)];
z = bFit./Stats.sebeta';
pvalu = 2*(1 - normcdf(abs(z)));
bFitw = NaN(max(series),1);
pvalue = NaN(max(series),1);
for i = 1:length(series)
    bFitw(series(i)) = bFit(i);
    pvalue(series(i)) = pvalu(i);
end
clear pvalu bFit
bFitCheck = reshape(bFitw,length(bFitw),1);
for t = 1:1:length(xD(:,1))
    modelFun = @(b,xD,MaX) (b(1).*xD(:,1)) + (b(2).*xD(:,2)) + (b(3).*xD(:,3)) + (b(4).*xD(:,4)) + (b(5).*xD(:,5))+ b(6);
    expect = modelFun(bFitCheck,xD,MaX);
end
mdl = LinearModel.fit(expect,Outs.yD);
tbl = anova(mdl);
parameters.Beta = bFitw;
%parameters.Stats = stats;
parameters.Beta_PValue = pvalue;
parameters.Rsquare = mdl.Rsquared.Adjusted;
parameters.FvalueModel = tbl.F(1);
parameters.PvalueModel = tbl.pValue(1);
parameters.DFModel = tbl.DF(2);

clear Ins Outs pvalue stats xD MaX modelFunw modelFun Group Options xDInt
clear tbl mdl bFitw
%% The Testing part
Ins = MakeIns(PropRadi,Numbers_testing,List,trait,traittype);
Outs = RegressFunc(Ins);% Here are the Marginal values calculated
ToTestyDMarg = Outs.yD;
ToTestyC = Outs.yC;
xD = [log10(Outs.XC+1),log10(Outs.total_GIC +1),log10(Outs.LinRoadC+1),log10(Outs.LinOtherC+1),log10(Outs.AreaC+1)];
xD =[xD,Outs.CountryC];
xDOrg = xD;
ToTestyDMarg =  reshape(ToTestyDMarg,[],1);
datapoints = length(ToTestyDMarg);
%% Do logaritmic predictions and STATS
Beta = parameters.Beta;
for i = 1:1:length(ToTestyDMarg)
    yPredicted(i) =  (Beta(1).*xD(i,1)) + (Beta(2).*xD(i,2)) + (Beta(3).*xD(i,3)) + (Beta(4).*xD(i,4)) + (Beta(5).*xD(i,5)) + Beta(6);
end
yPredicted = reshape(yPredicted,[],1);
yPredic_SensBase = nanmean(yPredicted); % so this is 1 value for sensitivity analyses
TransferValue_DisLog = yPredicted ;
yPredicted = TransferValue_DisLog;
clear TransferValue_DisLog
% Marginal data
%Correlation to marginal values of prediction set
[mean_double_deviation, Tran] = DeviationCalc(ToTestyDMarg,yPredicted,datapoints);
correlation.Marginal.Rho = Tran.RHO;
correlation.Marginal.Pval = Tran.PVAL;
correlation.Marginal.Deviation = mean_double_deviation ;
clear mean_double_deviation Tran
% Orginal data
[mean_double_deviation, Tran] = DeviationCalc(ToTestyC,yPredicted,datapoints);
correlation.Data.Rho = Tran.RHO;
correlation.Data.Pval = Tran.PVAL;
correlation.Data.Deviation = mean_double_deviation ;
clear mean_double_deviation Tran
%% Collate results
PredictRun.correlation = correlation;
PredictRun.parameters = parameters; % beta's and stats
PredictRun.Ynumbers = [Outs.TestingIDsC,ToTestyC,ToTestyDMarg,yPredicted];
PredictRun.Xnumbers = [Outs.TestingIDsC,xD];
PredictRun.NumbersTested = Numbers_testing; %#ok<*STRNU>
%% Sensitivity_module
Sense_factor = List.sense_factor;
Sensi_predic(:,1) = SensPredic(xDOrg,ToTestyDMarg,Sense_factor,trait,Beta,yPredic_SensBase,0);
for CIt = 2:1:4
    Lister =  find(Outs.CountryC == (CIt-1));
    Sensi_predic(:,CIt) = SensPredic(xDOrg(Lister,:),ToTestyDMarg(Lister,:),Sense_factor,trait,Beta,yPredic_SensBase,0);
end
%clearvars -except loop List PropRadi trait traittype Landscape_toGet_* PredictRun Sensi_predic Res_Country
PredictRun.SensiChanges = Sensi_predic;

%% write results
cd 'folder_out'
name_file = ['topredict','_',mat2str(trait),'_',mat2str(traittype),'_',mat2str(loop)];
save(name_file,'PredictRun')
cd ..
end % function PredictLoop

function Ins = MakeIns(PropRadi,Numbers,List,trait,traittype)
Ins.Y =  PropRadi.(genvarname(char(List.Outputs(trait))))(Numbers,traittype);
Ins.Area =  PropRadi.SpeciesSegmentIDs.AreaGI(Numbers);
Ins.LinRoad = PropRadi.SpeciesSegmentIDs.RoadGI(Numbers);
Ins.LinOther = PropRadi.SpeciesSegmentIDs.OtherLinearGI(Numbers);
Ins.Country = PropRadi.SpeciesSegmentIDs.Country(Numbers);
Ins.Landscape = PropRadi.SpeciesSegmentIDs.Landscape(Numbers);
Ins.SegmentBinar = PropRadi.SpeciesSegmentIDs.Segment(Numbers);
Ins.total_GI = PropRadi.SpeciesSegmentIDs.TotalGIAdam(Numbers);
Ins.X = PropRadi.SpeciesSegmentIDs.MeanCattle(Numbers);
Ins.TestingIDs = PropRadi.SpeciesSegmentIDs.SegmentID(Numbers);
Ins.TestingIDs = reshape(Ins.TestingIDs,[],1);
List = find(isnan(Ins.Y)==1);
if isempty(List)~= 1
Ins.Y(List) = [];
Ins.Area(List) = [];
Ins.LinRoad(List) = [];
Ins.LinOther(List) = [];
Ins.Country(List) = [];
Ins.Landscape(List) = [];
Ins.SegmentBinar(List) = [];
Ins.total_GI(List) = [];
Ins.X(List) = [];
Ins.TestingIDs(List) = [];
Ins.TestingIDs = reshape(Ins.TestingIDs,[],1);
end
end % function MakeIns

function [mean_double_deviation, Tran] = DeviationCalc(InValues,Predicted,datapoints)
[Tran.RHO,Tran.PVAL] = corr(InValues, Predicted,'type','Spearman');
x_range =InValues;
x_range_perc = prctile(x_range,95);
x_range_norm = x_range./x_range_perc;
x_range_norm(x_range_norm>1) = 1;
y_range = Predicted(:,1);
y_range_perc = prctile(y_range,95);
y_range_norm = y_range./y_range_perc;
y_range_norm(y_range_norm>1) = 1;
mean_double_deviation = 1- ((sum(abs(y_range_norm-x_range_norm)))/datapoints);
end

function Testing_landscapes = ClrTestLand(Testing_landscapes,Training_landscapes)
for i = 6:-1:1
    tst = find(Testing_landscapes == Training_landscapes(i));
    Testing_landscapes(tst) = []; %#ok<*FNDSB>
end
end

function [VarOutDis,VarOutGI] = InterRegio(Outs,VarIn)
Nr_C1 = find(Outs.CountryC ==1);
Nr_C2 = find(Outs.CountryC ==2);
Nr_C3 = find(Outs.CountryC ==3);
VarOutDis = ones(size(VarIn,1),3);
VarOutGI = ones(size(VarIn,1),3);
VarOutDis(Nr_C1,1) = VarIn(Nr_C1,1);
VarOutDis(Nr_C2,2) = VarIn(Nr_C2,1);
VarOutDis(Nr_C3,3) = VarIn(Nr_C3,1);
VarOutGI(Nr_C1,1) = VarIn(Nr_C1,2);
VarOutGI(Nr_C2,2) = VarIn(Nr_C2,2);
VarOutGI(Nr_C3,3) = VarIn(Nr_C3,2);
end

function Sensi_predic = SensPredic(xDOrg,ToTestyDMarg,Sense_factor,trait,Beta,~,Cor_factor)
for i = 1:1:length(ToTestyDMarg)
    yPredicted(i) = (Beta(1).*xDOrg(i,1)) + (Beta(2).*xDOrg(i,2)) + (Beta(3).*xDOrg(i,3)) + (Beta(4).*xDOrg(i,4)) +...
        (Beta(5).*xDOrg(i,5)) + Beta(6) - Cor_factor;
end
yPredic_SensBase = nanmean(yPredicted);

xSens = xDOrg;
basemode = 1;
if trait == [3,4,13:16,21,22] % so a categorical or shannon trait
    basemode = 0;
end
Sensi_predic = zeros(1,5);
% Distance (X)
factor = 1;
chance = ((10.^xSens(:,1))-1).*(1-Sense_factor);
xSens(:,1) =  log10(chance+1);
Sensi_predic = SensPredicFunc(Sensi_predic,factor,ToTestyDMarg,Sense_factor,Beta,xSens,basemode,yPredic_SensBase,Cor_factor);
xSens = xDOrg; % set varables back

%TotalGI (Tot)
factor = 2;
chance  = ((10.^xSens(:,factor))-1).*(1+Sense_factor);
xSens(:,factor) =  log10(chance+1);
Sensi_predic = SensPredicFunc(Sensi_predic,factor,ToTestyDMarg,Sense_factor,Beta,xSens,basemode,yPredic_SensBase,Cor_factor);
xSens = xDOrg; % set varables back

%RoadGI (Road)
factor = 3;
xSens = AddGI(xSens,factor,Sense_factor);
Sensi_predic = SensPredicFunc(Sensi_predic,factor,ToTestyDMarg,Sense_factor,Beta,xSens,basemode,yPredic_SensBase,Cor_factor);
xSens = xDOrg; % set varables back
%
%Other Linear GI (Other)
factor = 4;
xSens = AddGI(xSens,factor,Sense_factor);
Sensi_predic = SensPredicFunc(Sensi_predic,factor,ToTestyDMarg,Sense_factor,Beta,xSens,basemode,yPredic_SensBase,Cor_factor);
xSens = xDOrg; % set varables back

%AreaGI (Area)
factor = 5;
xSens = AddGI(xSens,factor,Sense_factor);
Sensi_predic = SensPredicFunc(Sensi_predic,factor,ToTestyDMarg,Sense_factor,Beta,xSens,basemode,yPredic_SensBase,Cor_factor);
xSens = xDOrg; % set varables back
end % Funtion SensiPredic

function Sensi_predic = SensPredicFunc(Sensi_predic,factor,ToTestyDMarg,Sense_factor,Beta,xCalc,basemode,yPredic_SensBase,Cor_factor)
for i = 1:1:length(ToTestyDMarg)
    yPredict(i) = (Beta(1).*xCalc(i,1)) + (Beta(2).*xCalc(i,2)) +...
        (Beta(3).*xCalc(i,3)) + (Beta(4).*xCalc(i,4)) + (Beta(5).*xCalc(i,5)) + Beta(6)- Cor_factor;
end
yPredic = nanmean(reshape(yPredict,[],1)); % so this is 1 value
Sensi_predic(factor) = abs((-(((abs(basemode-yPredic))./(abs(basemode-yPredic_SensBase)))-1)))./Sense_factor;
if isnan(Sensi_predic(factor))==1
    Sensi_predic(factor) = 0;
end
end % SensiPredicFunc

function xSens = AddGI(xSens,factor,Sense_factor)
Total_GIClean = ((10.^xSens(:,2))-1) - ((10.^xSens(:,factor))-1);
chance  = ((10.^xSens(:,factor))-1).*(1+Sense_factor);
xSens(:,factor) =  log10(chance+1);
xSens(:,2) = log10((Total_GIClean +chance)+1);% make sure aso total GI is elevated
end