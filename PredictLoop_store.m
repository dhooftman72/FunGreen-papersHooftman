function PredictLoop(loop,List,PropRadi,trait,traittype)
rng('shuffle')
Landscape_toGet_1 = randperm(12,6);
Landscape_toGet_2 = randperm(12,6)+12;
Landscape_toGet_3 = randperm(12,6)+24;

% Test for all landscapes and per country
for CIt = 1:1:4
    display(CIt)
    
if CIt == 1
    Training_landscapes = sort([Landscape_toGet_1,Landscape_toGet_2,Landscape_toGet_3]);
    Testing_landscapes = [1:36];
    Testing_landscapes(Training_landscapes) = [];
elseif CIt == 2
    Training_landscapes = sort(Landscape_toGet_1);
    Testing_landscapes = [1:12];
    Testing_landscapes = ClrTestLand(Testing_landscapes,Training_landscapes);
elseif CIt == 3
    Training_landscapes = sort(Landscape_toGet_2);
    Testing_landscapes = [13:24];
     Testing_landscapes = ClrTestLand(Testing_landscapes,Training_landscapes);
elseif CIt == 4
    Training_landscapes = sort(Landscape_toGet_3);
    Testing_landscapes = [25:36];
     Testing_landscapes = ClrTestLand(Testing_landscapes,Training_landscapes);
end

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
MaX = max(Outs.yD);
%% Do the multiple variable fitting
clear prior
prior = [Outs.Disrico,Outs.Totalrico,Outs.LinRoadrico,Outs.LinOtherrico,Outs.Arearico,Outs.constant];
series = [1 2 3 4 5 6];
modelFun = @(b,xD,MaX) (b(1).*xD(:,1)) + (b(2).*xD(:,2)) + (b(3).*xD(:,3)) + (b(4).*xD(:,4)) + (b(5).*xD(:,5))+ b(6);
modelFunw = @(b,xD)modelFun(b,xD,MaX);
Options = statset('FunValCheck','off','Display','off','MaxIter',100,...
    'TolFun',1.0000e-4, 'TolX',1.0e-4, 'Jacobian','off', 'DerivStep', 6.0555e-05, 'OutputFcn',' ');
Group = grp2idx(xD(:,1));%[1:length(yD)];
if  CIt == 1
[bFit,~,stats] = nlmefit(xD,Outs.yD,Group,[],modelFunw,prior,'RefineBeta0','off','Options',Options,'FEParamsSelect',series); %#ok<*NASGU>
else
[bFit,~,~,~,~]  = nlinfit(xD,Outs.yD,modelFunw,prior,'Options',Options);
stats.sebeta(6,1) = 1;
end
z = bFit./stats.sebeta';
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
    expect = modelFun(bFitCheck,xD,MaX);
end
mdl = LinearModel.fit(expect,Outs.yD);
display(mdl.Rsquared.Adjusted);
tbl = anova(mdl);
parameters.Beta = bFitw;
parameters.Stats = stats;
parameters.Beta_PValue = pvalue;
parameters.Rsquare = mdl.Rsquared.Adjusted;
parameters.FvalueModel = tbl.F(1);
parameters.PvalueModel = tbl.pValue(1);
parameters.DFModel = tbl.DF(2);
clear Ins Outs pvalue stats xD MaX modelFunw modelFun Group Options 
clear tbl mdl bFitw
%% The Testing part
Ins = MakeIns(PropRadi,Numbers_testing,List,trait,traittype);
Outs = RegressFunc(Ins);% Here are the Marginal values calculated
ToTestyDMarg = Outs.yD;
ToTestyC = Outs.yC;
xD = [log10(Outs.XC+1),log10(Outs.total_GIC +1),log10(Outs.LinRoadC+1),log10(Outs.LinOtherC+1),log10(Outs.AreaC+1)];
xDOrg = xD;
ToTestyDMarg =  reshape(ToTestyDMarg,[],1);
datapoints = length(ToTestyDMarg);
%% Do logaritmic predictions and STATS
Beta = parameters.Beta;
for i = 1:1:length(ToTestyDMarg)
    yPredicted(i) = (Beta(1).*xD(i,1)) + (Beta(2).*xD(i,2)) + (Beta(3).*xD(i,3)) + (Beta(4).*xD(i,4)) + (Beta(5).*xD(i,5)) + Beta(6) ;
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
if CIt == 1
    PredictRun.correlation = correlation;
    PredictRun.parameters = parameters; % beta's and stats
    PredictRun.Ynumbers = [Outs.TestingIDsC,ToTestyC,ToTestyDMarg,yPredicted];
    PredictRun.Xnumbers = [Outs.TestingIDsC,xD];
    PredictRun.NumbersTested = Numbers_testing; %#ok<*STRNU>
end
%% Sensitivity_module
Sense_factor = List.sense_factor;
Sensi_predic(:,CIt) = SensPredic(xDOrg,ToTestyDMarg,Sense_factor,trait,Beta,yPredic_SensBase);

clearvars -except loop List PropRadi trait traittype Landscape_toGet_* PredictRun Sensi_predic
end
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