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

[~,~,stats]= anovan((Outs.yC),{Outs.CountrOuts.yC,Outs.LandscapeC,(log10(Outs.XC+1)),log10(Outs.total_GIC +1)},'sstype',3,...
    'model','linear','nested',[0 0 0 0; 1 0 0 0; 0 0 0 0 ; 0 0 0 0 ],'display', 'off',...
    'continuous',[3,4], 'random', [1,2], 'varnames', {'Country', 'Landscape','Distance','Total_GI'});
[~,~,stats2]= anovan((stats.resid),{(log10(Outs.LinRoadC+1)),(log10(Outs.LinOtherC+1)),(log10(Outs.AreaC+1))},'sstype',3,...
    'model','linear','display', 'off',...
    'continuous',[1,2,3], 'varnames', {'RoadGI','OtherLinearGI','AreaGI'});
test = find(isnan(Outs.yC) ==1);
LandscapeCValue = Outs.LandscapeC;
CountrOuts.yCValue = Outs.CountrOuts.yC;
LandscapeCValue(test) = NaN;
CountrOuts.yCValue(test) = NaN;
clear test
Actual_Landscapes= sort(unique(LandscapeCValue(isnan(LandscapeCValue) ~= 1)));
Actual_Countries= sort(unique(CountrOuts.yCValue(isnan(LandscapeCValue) ~= 1)));
combined_length = 1+ (length(Actual_Countries) + length(Actual_Landscapes));
clear Value_Dis
for i = 1:length(Outs.yC)
    countrynr = (find(Actual_Countries ==  Outs.CountrOuts.yC(i))) +1;
    landscapenr = 1 + (find(Actual_Landscapes == Outs.LandscapeC(i))) + length(Actual_Countries);
    Disnr = 1 + combined_length;
    Outs.Disrico = stats.coeffs(Disnr);
    Totalnr = 2 + combined_length;
    Outs.Totalrico = stats.coeffs(Totalnr);
    Outs.LinRoadrico = stats2.coeffs(2);
    Outs.LinOtherrico = stats2.coeffs(3);
    Outs.Arearico = stats2.coeffs(4);
    if isnan(Outs.yC(i)) ~= 1
        Value_Dis(i,1) = stats.coeffs(1) + stats.coeffs(countrynr) + stats.coeffs(landscapenr)+...
            (Outs.Disrico.*(log10(Outs.XC(i)+1))) + (Outs.Totalrico.*log10(Outs.total_GIC(i)+1)) + (Outs.LinRoadrico .*(log10(Outs.LinRoadC(i)+1)))+...
            (Outs.LinOtherrico .*(log10(Outs.LinOtherC(i)+1))) +(Outs.Arearico.*(log10(Outs.AreaC(i)+1)));       
    else
        Value_Dis(i,1) = NaN;
    end
    clear landscapenr Disnr Linnr LinOthernr LinOthernr Areanr
end
Outs.yD= Value_Dis;
Outs.constant = stats.coeffs(1);
end