function Outs = RegressFunc(Ins)
% regression_predict
% Select Segements only.
nonzeros = find(Ins.SegmentBinar);
Outs.yC = Ins.Y(nonzeros);
Outs.XC = Ins.X(nonzeros);
Outs.LinRoadC = Ins.LinRoad(nonzeros);
Outs.LinOtherC = Ins.LinOther(nonzeros);
Outs.AreaC = Ins.Area(nonzeros);
Outs.CountrOuts.yC = Ins.Country(nonzeros);
Outs.LandscapeC = Ins.Landscape(nonzeros);
Outs.total_GIC = Ins.total_GI(nonzeros);
Outs.TestingIDsC = Ins.TestingIDs(nonzeros);

VarNames_1 = {'Country', 'Landscape', 'Disrtance','Total GI'};
VarNames_2 = {'Country','Road verges','Other Linear GI','Areal GI'};
% Two tiered GLM anovan
[~,~,stats]= anovan((Outs.yC),{Outs.CountrOuts.yC,Outs.LandscapeC,(log10(Outs.XC+1)),log10(Outs.total_GIC +1)},'sstype',3,...
    'model',[1 0 0 0; 0 1 0 0; 0 0 1 0 ; 0 0 0 1; 1 0 1 0; 1 0 0 1],'nested',[0 0 0 0; 1 0 0 0; 0 0 0 0 ; 0 0 0 0 ],'display', 'off',...
    'continuous',[3,4], 'random', [2], 'varnames', VarNames_1);
[~,~,stats2]= anovan((stats.resid),{(log10(Outs.LinRoadC+1)),(log10(Outs.LinOtherC+1)),(log10(Outs.AreaC+1))},'sstype',3,...
    'model',[0 1 0 0; 0 0 1 0 ; 0 0 0 1; 1 1 0 0; 1 0 1 0; 1 0 0 1],'display', 'off',...
    'continuous',[2,3,4], 'varnames', VarNames_2);
for i = 1:length(Outs.yC)
    if isnan(Outs.yC(i)) ~= 1
        Value_Dis(i,1) = Outs.yC(i)- stats2.resid(i);    %#ok<*AGROW>
    else
        Value_Dis(i,1) = NaN;
    end
    clear landscapenr Disnr Linnr LinOthernr LinOthernr Areanr
end
Outs.yD= Value_Dis;
Outs.constant = stats.coeffs(1);
end