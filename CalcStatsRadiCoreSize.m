function CalcStatsRadi
% Main GLM and differences tests
% Most of code is naming and can be replaced by other names
warning off 
load('Radiation_proportions.mat');
for trait = List.Traits_Species
    display(trait)
    %Initiate outputs
    DifferenceFocal = 0;
    BandMean = 0;
    BandSdev = 0;
    Distance.overview = 0;
    maxtrait = List.maxtraits(trait);
    for traittype = 1:maxtrait
        load('Radiation_proportions.mat');
        traittext =  Input.(genvarname([char(List.trait_txts(trait))]))(traittype,2);
        clc
        disp('Running...')
        disp(List.Outputs(trait))
        disp(traittext)        
        % Set which Y-values to use
        % correct for the species loss rate per se 
        Data.Yuncor = PropRadi.(genvarname([char(List.Outputs(trait))]))(:,traittype);
        Data.Yuncor2 = TotRadi.(genvarname([char(List.Outputs(trait))]))(:,traittype);
        catagorical = 0;
        if isempty(find(List.CatagoricalTraits == trait)) ~= 1 %#ok<*EFIND> % so a non-numeric trait
            catagorical = 1;
            display('Categorical')
            Data.y1 =  PropRadi.(genvarname([char(List.Outputs(trait))]))(:,traittype) - PropRadi.Species;
            Data.y2 = Data.Yuncor2;
        else
            Data.y1 = Data.Yuncor;
            Data.y2 = Data.Yuncor2;
        end
        if trait == 27 % correction for area
        Data.Yuncor = (bsxfun(@rdivide, Data.Yuncor',SizeProp'));
        Data.Yuncor = Data.Yuncor';
        Data.y1 = Data.Yuncor;
        Data.Yuncor2 = (bsxfun(@rdivide, Data.Yuncor2',SizeProp'));
        Data.Yuncor2 = Data.Yuncor2';
        Data.y2 = Data.Yuncor2;
        end

        Data.Y = Data.y1;
        clear  Data.Yuncor2
        % to show how much species are actually included dealing with proportional values
        Data.nrspecies = NrSpec.(genvarname([char(List.Outputs(trait))]))(:,traittype);
        Data.Species_loss = SpecLoss;
        % test for enough data
        %%
        tes = sum(isnan(Data.Y ));
        Regression.TooFewData = {'Dummy'};
        if tes < (length(Data.Y )/4) % test for enough data
            % Set data and test for NaN within Y;
            Data.total_GI = PropRadi.SpeciesSegmentIDs.TotalGIAdam;
            Data.Area =  PropRadi.SpeciesSegmentIDs.AreaGI;
            Data.LinRoad = PropRadi.SpeciesSegmentIDs.RoadGI;
            Data.LinOther = PropRadi.SpeciesSegmentIDs.OtherLinearGI;
            Data.Country = PropRadi.SpeciesSegmentIDs.Country;
            Data.Landscape = PropRadi.SpeciesSegmentIDs.Landscape;
            Data.SegmentBinar = PropRadi.SpeciesSegmentIDs.Segment;
            Data.Band = PropRadi.SpeciesSegmentIDs.Band;
            Data.IsSegment = PropRadi.SpeciesSegmentIDs.Segment;
            Data.CoreSize = PropRadi.SpeciesSegmentIDs.CoreSize;
            
            test = find(isnan(Data.Y) == 1);
            Data = removeTest(Data,test);
            to_remove = test; % save to remove within x as well
            clear test
            
            %Set in a loop which Distance data-set is the focus
            for DistanceType = 1:1:List.dis_max
                Regression.TooFewData = {'Dummy'};
                if DistanceType == 1
                    AnType = 'MeanCattle';
                    Data.Distance = PropRadi.SpeciesSegmentIDs.MeanCattle; % Mean distance
                elseif DistanceType == 2
                    AnType = 'Euclidian';
                    Data.Distance = PropRadi.SpeciesSegmentIDs.Band; % Mean distance
                elseif DistanceType == 3
                    AnType = 'MinCattle';
                    Data.Distance = PropRadi.SpeciesSegmentIDs.MinCattle; % Mean distanc
                elseif DistanceType == 4
                    AnType = 'MeanOpenness';
                    Data.Distance = PropRadi.SpeciesSegmentIDs.MeanWind; % Mean distance
                elseif DistanceType == 5
                    AnType = 'MinOpenness';
                    Data.Distance = PropRadi.SpeciesSegmentIDs.minWind; % Mean distance
                elseif DistanceType == 6
                    AnType = 'MeanHumans';
                    Data.Distance = PropRadi.SpeciesSegmentIDs.MeanHuman; % Mean distance
                elseif DistanceType == 7
                    AnType = 'MinHumans';
                    Data.Distance = PropRadi.SpeciesSegmentIDs.MinHuman; % Mean distance
                elseif DistanceType == 8
                    AnType = 'NoDistance';
                    Data.Distance = PropRadi.SpeciesSegmentIDs.NoDis; % Mean distance
               end

                Data.Distance(to_remove) = []; % removal of NaN Y values, see above
                test = find(isnan(Data.Distance) == 1);
                Data.Distance(test) = [];
                removeTest(Data,test);
                clear test
                
               %% Distinguish between Focal patches and Segments; test on
                % difference among them.
                focals = find(Data.SegmentBinar==0);
                segments = find(Data.SegmentBinar);
                specs = (Data.y2(focals));
                countries_focal = Data.Country(focals);
                countries_segments = Data.Country(segments);
                YvalueFocal(1) = nanmean(specs);
                YvalueFocal(2) = nanmean(specs(countries_focal==1));
                YvalueFocal(3) = nanmean(specs(countries_focal==2));
                YvalueFocal(4) = nanmean(specs(countries_focal==3));
                YvalueFocal(5) = nanmean([YvalueFocal(2), YvalueFocal(3), YvalueFocal(4)]);
                YvalueFocalStd(1) = nanstd(specs);
                YvalueFocalStd(2) = nanstd(specs(countries_focal==1));
                YvalueFocalStd(3) = nanstd(specs(countries_focal==2));
                YvalueFocalStd(4) = nanstd(specs(countries_focal==3));
                YvalueFocalStd(5) = nanmean([YvalueFocalStd(2), YvalueFocalStd(3), YvalueFocalStd(4)]);
                Y1 = Data.Yuncor(segments);
                Data.nrspeciesl = Data.nrspecies(segments);
                Y1(isnan(Y1)==1) = [];
                Data.Species_lossL = Data.Species_loss(segments);
                ds = dataset(Data.Yuncor,'Varnames','Yuncor');
                ds.Country = ordinal(Data.Country);
                ds.Landscape = ordinal(Data.Landscape);
                ds.CoreSize = Data.CoreSize;
                ds.Segments = logical(Data.IsSegment);
                [~,outs,stats] = anovan(ds.Yuncor,{ds.Country,ds.CoreSize,ds.Segments},'sstype',1,...
                    'model',[1 0 0; 0 1 0; 0 0 1;1 0 1],'continuous',[2], 'display', 'off',...
                    'varnames', {'Country','CoreSize','Segments'});
                save('all')
                Value_Dis = Data.Yuncor - stats.resid;

                Significance_segment(1) = outs(4,6);
                Significance_segment(2) = outs(4,7);
                Significance_segment(3) = outs(5,6);
                Significance_segment(4) = outs(5,7);
                [~,ps] = jbtest(stats.resid);
                Significance_segment(5) = {ps};
                Data.Spec_target_All(:,1) = Value_Dis(segments);
                Data.Spec_target_C1(:,1) = Value_Dis(segments(countries_segments==1));
                Data.Spec_target_C2(:,1) = Value_Dis(segments(countries_segments==2));
                Data.Spec_target_C3(:,1) = Value_Dis(segments(countries_segments==3));
                Data.Countries(1:3,1) = stats.coeffs(2:4);
                if   cell2mat(Significance_segment(2)) <= 0.05
                    if stats.coeffs(length(stats.coeffs)-6) < 0
                        targetDif = {'Lower Segment values'};
                    else
                        targetDif = {'Higher Segment values'};
                    end
                else
                    targetDif = {'Target equal to segments'};
                end
                Stats.Difference.outs = outs;
                Stats.Difference.stats = stats;
                clear stats outs ps
                
                higher_txt = {'Higher than average loss'};
                lower_txt  = {'Lower than average loss'};
                if List.species_option >= 4
                    higher_txt = {'Higher than average gain'};
                    lower_txt  = {'Lower than average gain'};
                end
                [tested,p,~,stats]   = ttest2(Data.Species_lossL,Y1);
                
                if tested == 1
                    if mean(Data.Species_lossL) > mean(Y1)
                        targetDif_2 = higher_txt;
                    else
                        targetDif_2 =   lower_txt;
                    end
                else
                    targetDif_2 = {'Equal to average'};
                end
                targetDif_2_P = p;
                targetDif_2_T = stats.tstat;
                clear Value_Dis
                
                % Save all values in cell arrays of All segments and focal
                AllValues = dataset(Data.Distance,'Varnames',char(AnType));
                AllValues.TotalGI = Data.total_GI;
                AllValues.LinRoad = Data.LinRoad;
                AllValues.LinOther = Data.LinOther;
                AllValues.AreaGI = Data.Area;
                AllValues.Country = Data.Country;
                AllValues.Landscape = Data.Landscape;
                AllValues.CoreSite = Data.CoreSize;
                AllValues.Band = Data.Band;
                AllValues.Ys = Data.Y;
                AllValues.Yuncorrected = Data.Yuncor;
                AllValues.NrSpecies = Data.nrspecies;
                Regression.(char(AnType)).AllValues = AllValues;
                clear  AllValues
                
                % Select Segements only && Save all values in cell arrays of segments only (no
                % focals)
                yC = Data.Y(segments);
                DistanceC = Data.Distance(segments);
                total_GIC = Data.total_GI(segments);
                LinRoadC = Data.LinRoad(segments);
                LinOtherC = Data.LinOther(segments);
                AreaC = Data.Area(segments);
                CountryC = Data.Country(segments);
                LandscapeC = Data.Landscape(segments);
                CoreSizeC = Data.CoreSize(segments);
                
                Segments = dataset(DistanceC,'Varnames',char(AnType));
                Segments.TotalGI = total_GIC;
                Segments.LinRoad = LinRoadC;
                Segments.LinOther = LinOtherC;
                Segments.AreaGI = AreaC;
                Segments.Country = CountryC;
                Segments.Landscape = LandscapeC;
                Segments.CoreSize = CoreSizeC;
                Segments.Ys = yC;
                Segments.Band = Data.Band(segments);
                Segments.Yuncorrected = Data.Yuncor(segments);
                Segments.NrSpecies = Data.nrspecies(segments); 
                Regression.(char(AnType)).SegmentValues =  Segments;
                clear Segments 
                   
                %% Logaritmic GLM within segments only for Distance, Linear GI and Area
                % GI
                Total_GIactual = total_GIC;
                Roads_Actual = LinRoadC;
                Linear_Other_actual = LinOtherC;
                Area_GI_Actual = AreaC;
                VarNames_1 = {'Country', 'CoreSize'};
                VarNames_1a = {'Country', AnType,'Total GI'};
                VarNames_2 = {'Country','Road verges','Other Linear GI','Areal GI'};
                Anova_Table = dataset({'Dummy'},'Varnames',char('Source'));
                % Three tiered GLM anovan
                [~,outs,stats]= anovan((yC),{CountryC,CoreSizeC},'sstype',1,...
                    'model',[1 0; 0 1],'nested',[0 0; 1 0], 'display', 'off',...
                    'continuous',[2], 'varnames', VarNames_1);
                Anova_Table.Source([1,2,13,14],1) =  outs(2:5,1);
                Anova_Table.SumSq([1,2,13,14],1) = cell2mat(outs(2:5,2));
                Anova_Table.DFs([1,2,13,14],1) = cell2mat(outs(2:5,3));
                Anova_Table.MeanSq([1,2,13],1) = cell2mat(outs(2:4,5));
                Anova_Table.F([1,2],1) = cell2mat(outs(2:3,6));
                Anova_Table.Pvalue([1,2],1) = cell2mat(outs(2:3,7));
                Stats.AnovaTableT1.outs = outs;
                Stats.AnovaTableT1.stats = stats;
                
                [~,outs1a,stats1a]= anovan((stats.resid),{CountryC,(log10(DistanceC+1)),(log10(Total_GIactual +1))},'sstype',3,...
                    'model',[0 1 0; 0 0 1; 1 1 0; 1 0 1],'display', 'off',...
                    'continuous',[2,3], 'varnames', VarNames_1a);
                Anova_Table.Source([3,4,8,9],1) =  outs1a(2:5,1);
                Anova_Table.SumSq([3,4,8,9],1) = cell2mat(outs1a(2:5,2));
                Anova_Table.DFs([3,4,8,9],1) = cell2mat(outs1a(2:5,3));
                Anova_Table.MeanSq([3,4,8,9],1) = cell2mat(outs1a(2:5,5));
                Anova_Table.F([3,4,8,9],1) = cell2mat(outs1a(2:5,6));
                Anova_Table.Pvalue([3,4,8,9],1) = cell2mat(outs1a(2:5,7));
                Stats.AnovaTableT1a.outs = outs1a;
                Stats.AnovaTableT1a.stats = stats1a;
                
                [~,outs2,stats2]= anovan((stats1a.resid),{CountryC,(log10(Roads_Actual+1)),(log10(Linear_Other_actual+1)),(log10(Area_GI_Actual+1))},'sstype',3,...
                    'model',[0 1 0 0; 0 0 1 0 ; 0 0 0 1; 1 1 0 0; 1 0 1 0; 1 0 0 1],'display', 'off',...
                    'continuous',[2,3,4], 'varnames', VarNames_2);
                Anova_Table.Source([5,6,7,10,11,12],1) = outs2(2:7,1);
                Anova_Table.SumSq([5,6,7,10,11,12],1) = cell2mat(outs2(2:7,2));
                Anova_Table.DFs([5,6,7,10,11,12],1) = cell2mat(outs2(2:7,3));
                Anova_Table.MeanSq([5,6,7,10,11,12],1) = cell2mat(outs2(2:7,5));
                Anova_Table.F([5,6,7,10,11,12],1) = cell2mat(outs2(2:7,6));
                Anova_Table.Pvalue([5,6,7,10,11,12],1) = cell2mat(outs2(2:7,7));
                Stats.AnovaTableT2.outs = outs2;
                Stats.AnovaTableT2.stats = stats2;

%                % take the right coefficient row if not all landcsapes are
%                 % included.
%                 test = find(isnan(yC) ==1);
%                 LandscapeCValue = LandscapeC;
%                  CountryCValue = CountryC;  
%                 LandscapeCValue(test) = NaN; %#ok<*FNDSB>
%                  CountryCValue(test) = NaN;
%                 Actual_Data.Landscapes= sort(unique(LandscapeCValue(isnan(LandscapeCValue) ~= 1)));
%                 Actual_Data.Countries= sort(unique(CountryCValue(isnan(LandscapeCValue) ~= 1)));
%                  clear test

                %Disnr = 5 + length(Actual_Data.Landscapes);
                Disrico = stats1a.coeffs(2);
                %Totalnr = 6 + length(Actual_Data.Landscapes);
                Totalrico = stats1a.coeffs(3);
                LinRoadrico = stats2.coeffs(2);
                LinOtherrico = stats2.coeffs(3);
                Arearico = stats2.coeffs(4);
                Anova_Table.Rico([3,4,5,6,7],1) = [Disrico;Totalrico;LinRoadrico; LinOtherrico; Arearico];
                 
                % Make marginal values
                Value_Dis_FirstPart = yC - stats.resid; % (input - residuals)
                Value_Dis = stats.resid - stats2.resid;  % (input - residuals)
                
                % Model R2 calculatiom             
                    expect(:,1) = yC - stats2.resid;
                meanY = mean(yC);
                for t = 1:1:length(yC)
                    ssres(t) = ((expect(t)-yC(t)).^2);
                    sstot(t) =  ((yC(t)-meanY).^2);
                end
                Anova_Table.RSquared([14],1) = 1- ((sum(ssres)) /(sum(sstot)));
                Anova_Table.Source(15,1) = {'NormalityResidT1'};
                [~,Anova_Table.Pvalue(15,1)] = jbtest(stats.resid);
                Anova_Table.Source(16,1) = {'NormalityResidT2'};
                [~,Anova_Table.Pvalue(16,1)] = jbtest(stats2.resid);
                
                clear expect meanY
                clear landscapenr Disnr Linnr Totalrico LinOthernr LinOthernr
                clear  Areanr Disrico LinRoadrico LinOtherrico Arearico
                
                % Allow to see how strong the correlations are in trems of
                % effects
                mdl = LinearModel.fit([(log10(DistanceC+1))],Value_Dis_FirstPart);
                DisR(1) = mdl.Rsquared.Adjusted;
                DisR(2) = mdl.Coefficients.pValue(2);
                Constants(1,1) = mdl.Coefficients.Estimate(1);
                Constants(1,2) = mdl.Coefficients.Estimate(2);
                clear mdl
                mdl = LinearModel.fit((log10(Total_GIactual +1)),Value_Dis_FirstPart);
                TotalGIR(1) = mdl.Rsquared.Adjusted;
                TotalGIR(2) = mdl.Coefficients.pValue(2);
                Constants(2,1) = mdl.Coefficients.Estimate(1);
                Constants(2,2) = mdl.Coefficients.Estimate(2);
                clear mdl
                mdl = LinearModel.fit((log10(Roads_Actual+1)),Value_Dis);
                LinRoadR(1) = mdl.Rsquared.Adjusted;
                LinRoadR(2) = mdl.Coefficients.pValue(2);
                Constants(3,1) = mdl.Coefficients.Estimate(1);
                Constants(3,2) = mdl.Coefficients.Estimate(2);
                clear mdl
                mdl = LinearModel.fit((log10(Linear_Other_actual+1)),Value_Dis);
                LinOtherR(1) = mdl.Rsquared.Adjusted;
                LinOtherR(2) = mdl.Coefficients.pValue(2);
                Constants(4,1) = mdl.Coefficients.Estimate(1);
                Constants(4,2) = mdl.Coefficients.Estimate(2);
                clear mdl
                mdl = LinearModel.fit((log10(Area_GI_Actual+1)),Value_Dis);
                AreaR(1) = mdl.Rsquared.Adjusted;
                AreaR(2) = mdl.Coefficients.pValue(2);
                Constants(5,1) = mdl.Coefficients.Estimate(1);
                Constants(5,2) = mdl.Coefficients.Estimate(2);
                clear mdl
                Anova_Table.RSquared([3,4,5,6,7],1) = [DisR(1); TotalGIR(1); LinRoadR(1);LinOtherR(1); AreaR(1)];
                Anova_Table.PvalueLine([3,4,5,6,7],1) = [DisR(2); TotalGIR(2); LinRoadR(2);LinOtherR(2); AreaR(2)];
                clear  DisR TotalGIR LinRoadR LinOtherR AreaR
                
                % Renove the Zero's
                Radjusted = [];
                Padjusted = [];
                %First tier
                Var2 = Value_Dis_FirstPart;
                Var1 = (log10(DistanceC+1));
                [Radjusted,Padjusted] = MargReg(Var1,Var2,1,Radjusted,Padjusted);
                Var1 = (log10( Total_GIactual +1));
                [Radjusted,Padjusted] = MargReg(Var1,Var2,2,Radjusted,Padjusted);
                % Secodn tier
                Var2 = Value_Dis;
                Var1 = (log10(Roads_Actual+1));
                [Radjusted,Padjusted] = MargReg(Var1,Var2,3,Radjusted,Padjusted);
                Var1 = (log10(Linear_Other_actual+1));
                [Radjusted,Padjusted] = MargReg(Var1,Var2,4,Radjusted,Padjusted);
                Var1 = (log10(Area_GI_Actual+1));
                [Radjusted,Padjusted] = MargReg(Var1,Var2,5,Radjusted,Padjusted);
                
                Anova_Table.Rnozeross([3,4,5,6,7],1) =  Radjusted;
                Anova_Table.Pnozeross([3,4,5,6,7],1) =  Padjusted;
                clear Radjusted Padjusted
                
                Distance.overview(DistanceType,traittype) =  Anova_Table.RSquared(3); %#ok<*AGROW>
                Distance.overviewP(DistanceType,traittype) =  Anova_Table.Pvalue(3);
                Distance.overviewEffect(DistanceType,traittype) = Anova_Table.Rico(3);
                Distance.overviewR2model(DistanceType,traittype) = Anova_Table.RSquared(14);
                Regression = rmfield(Regression,{'TooFewData'});
                Name(DistanceType,1) = {AnType};
                Regression.(char(AnType)).Anova_Table =  Anova_Table;
                Regression.(char(AnType)).RegLines = Constants;
                Regression.(char(AnType)).SegmentValues.FirstRegFit =   Value_Dis_FirstPart;
                Regression.(char(AnType)).SegmentValues.SecondRegFit =  Value_Dis;
                if DistanceType == 1
                 Regression.(char(AnType)).Stats = Stats;
                end
                
                if length(specs) == 36
                    for ld = 1:1:36
                        lst =  find(LandscapeC == ld);
                        focalValue = specs(ld);
                        Regression.(char(AnType)).SegmentValues.FocalValue(lst,1) = focalValue; %#ok<FNDSB>
                        clear lst focalValue
                    end
                end
                clear Anova_Table AllValues outs stats outs stats2 Stats
                clear Actual_Data.Landscapes LandscapeCValue

                %% Collate data
                % Preset Output cell array, so it shows as a data-set
                if DistanceType == 1
                    if traittype == 1
                        clear DifferenceFocal
                        DifferenceFocal = dataset({'Dummy';'Dummy'},'Varnames',char('Trait_name'));
                        DifferenceFocal.Focal_value = [NaN;NaN];
                        DifferenceFocal.StdFocal_value = [NaN;NaN];
                        DifferenceFocal.Focal_Segment_Difference = {'Dummy';'Dummy'};
                        DifferenceFocal.F_Difference = [NaN;NaN];
                        DifferenceFocal.P_Value_Difference = [NaN;NaN];
                        DifferenceFocal.SegmentMean = [NaN;NaN];
                        DifferenceFocal.SegmentStD = [NaN;NaN];
                        if  catagorical == 1;
                            DifferenceFocal.vs_Avg_Loss = {'Dummy';'Dummy'};
                            DifferenceFocal.vs_Avg_Loss_Tstat = [NaN;NaN];
                            DifferenceFocal.vs_Avg_Loss_P =[NaN;NaN];
                            DifferenceFocal.vs_Avg_Loss_Mean =[NaN;NaN];
                            DifferenceFocal.vs_Avg_Loss_std = [NaN;NaN];
                            DifferenceFocal.Avg_Spec_loss = [NaN;NaN];
                        end 
                        clear BandMean
                        clear  BandSdev
                        BandMean = dataset({'Dummy';'Dummy'},'Varnames',char('Trait_name'));
                        BandSdev = dataset({'Dummy';'Dummy'},'Varnames',char('Trait_name'));
                        BandMean.Focal = [NaN;NaN];
                        BandSdev.Focal =[NaN;NaN];
                    end
                    DifferenceFocal.Trait_name(traittype) = traittext;
                    DifferenceFocal.Focal_value(traittype) = YvalueFocal(5);
                    DifferenceFocal.StdFocal_value(traittype) = YvalueFocalStd(5);
                    DifferenceFocal.Focal_Segment_Difference(traittype) = targetDif;
                    DifferenceFocal.F_Difference(traittype) = cell2mat(Significance_segment(1));
                    DifferenceFocal.P_Value_Difference(traittype) = cell2mat(Significance_segment(2));
                    
                    %Data.SegmentDifference = Data.Spec_target_l(:,:)-1;
                    DifferenceFocal.SegmentMean(traittype) = nanmean([(nanmean(Data.Spec_target_C1-1)),(nanmean(Data.Spec_target_C2-1)),(nanmean(Data.Spec_target_C3-1))]);
                    DifferenceFocal.SegmentStD(traittype) = nanmean([(nanstd(Data.Spec_target_C1-1)),(nanstd(Data.Spec_target_C2-1)),(nanstd(Data.Spec_target_C3-1))]);
                    DifferenceFocal.P_NormalityResiduals(traittype)  = cell2mat(Significance_segment(5));
                    DifferenceFocal.F_Interaction(traittype) = cell2mat(Significance_segment(3));
                    DifferenceFocal.P_Value_InteractionDifference(traittype) = cell2mat(Significance_segment(4));
                    DifferenceFocal.Country_1_Mean(traittype,1) = DifferenceFocal.SegmentMean(traittype) - Data.Countries(1);
                    DifferenceFocal.Country_2_Mean(traittype,1) = DifferenceFocal.SegmentMean(traittype) - Data.Countries(2);
                    DifferenceFocal.Country_3_Mean(traittype,1) = DifferenceFocal.SegmentMean(traittype) - Data.Countries(3);
                    DifferenceFocal.Country_1_STD(traittype,1) = nanstd(Data.Spec_target_C1-1);
                    DifferenceFocal.Country_2_STD(traittype,1) = nanstd(Data.Spec_target_C2-1);
                    DifferenceFocal.Country_3_STD(traittype,1) = nanstd(Data.Spec_target_C3-1);
                    DifferenceFocal.Country_1_Focal_value(traittype) = YvalueFocal(2);
                    DifferenceFocal.Country_2_Focal_value(traittype) = YvalueFocal(3);
                    DifferenceFocal.Country_3_Focal_value(traittype)= YvalueFocal(4);
                    DifferenceFocal.Country_1_Focal_valueStd(traittype) = YvalueFocalStd(2);
                    DifferenceFocal.Country_2_Focal_valueStd(traittype) = YvalueFocalStd(3);
                    DifferenceFocal.Country_3_Focal_valueStD(traittype)= YvalueFocalStd(4);
                    DifferenceFocal.SegmentSpeciesMean(traittype) = nanmean(Data.nrspeciesl);
                    DifferenceFocal.SegmentSpeciesStD(traittype) = nanstd(Data.nrspeciesl);
                    if  catagorical == 1;
                        DifferenceFocal.Avg_Spec_loss(traittype) = mean(Data.Species_lossL)-1;
                        DifferenceFocal.vs_Avg_Loss(traittype) = targetDif_2;
                        DifferenceFocal.vs_Avg_Loss_Tstat(traittype) = targetDif_2_T;
                        DifferenceFocal.vs_Avg_Loss_P(traittype) = targetDif_2_P;
                        DifferenceFocal.vs_Avg_Loss_Mean(traittype) = (1+DifferenceFocal.SegmentMean(traittype))./mean(Data.Species_lossL);
                        DifferenceFocal.vs_Avg_Loss_std(traittype) = DifferenceFocal.SegmentStD(traittype);
                    end
                    if maxtrait == 1
                        DifferenceFocal(2,:) = [];
                    end
                    %% Overall segment Means
                    BandMean.Trait_name(traittype) = traittext;
                    BandSdev.Trait_name(traittype) = traittext;
                    BandMean.Focal(traittype,1) = YvalueFocal(1);
                    BandSdev.Focal(traittype,1) = YvalueFocalStd(1);
                    %display(Regression.(char(AnType)).SegmentValues)
                    for band = 300:300:1500
                        SegmentsBandID = find(Regression.(char(AnType)).SegmentValues.Band == band);
                        BandMean.(genvarname(mat2str(band)))(traittype) = mean(Value_Dis_FirstPart(SegmentsBandID));
                        BandSdev.(genvarname(mat2str(band)))(traittype) = std(Value_Dis_FirstPart(SegmentsBandID));
                    end
                    clear band
                    clear Value_Dis
                    clear Value_Dis_FirstPart
                end
            end % Distance type loops
        end % if enough data ("tes")
        if exist('Output_Anova.mat') ~= 0
            load('Output_Anova.mat')
        end
        Output_Anova.(genvarname([char(List.Outputs(trait))])).(genvarname([char(traittext)])) = Regression;
        if traittype == maxtrait
            Output_Anova.(genvarname([char(List.Outputs(trait))])).Differences = DifferenceFocal;
            BandMeans.Mean = BandMean;
            BandMeans.Std = BandSdev;
            Output_Anova.(genvarname([char(List.Outputs(trait))])).BandMeans = BandMeans;
            Distances.overview = dataset(Distance.overview(:,1),'ObsNames',Name,'VarNames',Input.(genvarname([char(List.trait_txts(trait))]))(1,2));
            Distances.overview_P  = dataset(Distance.overviewP(:,1),'ObsNames',Name,'VarNames',Input.(genvarname([char(List.trait_txts(trait))]))(1,2));
            Distances.overview_Effect = dataset(Distance.overviewEffect(:,1),'ObsNames',Name,'VarNames',Input.(genvarname([char(List.trait_txts(trait))]))(1,2));
            Distances.overview_R2Model = dataset(Distance.overviewR2model(:,1),'ObsNames',Name,'VarNames',Input.(genvarname([char(List.trait_txts(trait))]))(1,2));
            Rank = dataset((int8(tiedrank(1-(Distance.overviewR2model(1:7,1))))),'ObsNames',Name(1:7),'VarNames',Input.(genvarname([char(List.trait_txts(trait))]))(1,2));
            for i = 2:maxtrait
                Distances.overview.(genvarname([char(Input.(genvarname([char(List.trait_txts(trait))]))(i,2))])) = Distance.overview(:,i); %#ok<*NBRAK>
                Distances.overview_P.(genvarname([char(Input.(genvarname([char(List.trait_txts(trait))]))(i,2))])) = Distance.overviewP(:,i);
                Distances.overview_Effect.(genvarname([char(Input.(genvarname([char(List.trait_txts(trait))]))(i,2))])) = Distance.overviewEffect(:,i);
                Distances.overview_R2Model.(genvarname([char(Input.(genvarname([char(List.trait_txts(trait))]))(i,2))])) = Distance.overviewR2model(:,i);
                Rank.(genvarname([char(Input.(genvarname([char(List.trait_txts(trait))]))(i,2))])) = int8(tiedrank(1-(Distance.overviewR2model(1:7,i))));
            end
            clear i
            Output_Anova.(genvarname([char(List.Outputs(trait))])).DistanceOverview.R2 =Distances.overview;
            Output_Anova.(genvarname([char(List.Outputs(trait))])).DistanceOverview.Pvalue = Distances.overview_P;
            Output_Anova.(genvarname([char(List.Outputs(trait))])).DistanceOverview.EffectSize = Distances.overview_Effect;
            Output_Anova.(genvarname([char(List.Outputs(trait))])).DistanceOverview.ModelR2 = Distances.overview_R2Model;
             Output_Anova.(genvarname([char(List.Outputs(trait))])).DistanceOverview.Rankings = Rank;
        end
        save('Output_Anova','Output_Anova')
        save('trait_transfer','List','trait','traittype','maxtrait','DifferenceFocal','Distance','Name','BandMean','BandSdev','Input')       
        clear Regression
        clear all
        clear all
        load('trait_transfer.mat');
        delete('trait_transfer.mat');
    end % trait type
    clear DifferenceFocal
    clear Distance
    clear traittype
end % all traits
clear trait
clear maxtrait
end

function [Radjusted,Padjusted] = MargReg(Var1,Var2,i,Radjusted,Padjusted)
test = find(Var1);
Var1 = Var1(test);
Var2 = Var2(test);
mdl = LinearModel.fit(Var1,Var2);
Radjusted(i,1) = mdl.Rsquared.Adjusted;
Padjusted(i,1) = mdl.Coefficients.pValue(2);
end

function Data = removeTest(Data,test)
if exist('Org') ==1
    Data.Y = Org.y;
    Data.Yuncor = Org.Yuncor;
    Data.LinRoad = Org.LinRoad;
    Data.LinOther = Org.LinOther;
    Data.Area = Org.Area;
    Data.Band = Org.Band;
    Data.Country = Org.Country;
    Data.Landscape = Org.Landscape;
    Data.nrspecies = Org.nrspecies;
    Data.Species_loss = Org.Species_loss;
    Data.SegmentBinar = Org.SegmentBinar;
    Data.total_GI = Org.total_GI;
    Data.IsSegment = Org.IsSegment;
    Data.CoreSize = Org.CoreSize;
else
    Org.y = Data.Y;
    Org.Yuncor = Data.Yuncor;
    Org.LinRoad = Data.LinRoad;
    Org.LinOther = Data.LinOther;
    Org.Area = Data.Area;
    Org.Country = Data.Country;
    Org.Landscape = Data.Landscape;
    Org.nrspecies = Data.nrspecies;
    Org.Species = Data.Species_loss;
    Org.SegmentBinar = Data.SegmentBinar;
    Org.Band = Data.Band;
    Org.IsSegment = Data.IsSegment;
    Org.total_GI = Data.total_GI;
   Org.CoreSize = Data.CoreSize;
end
Data.Y(test) = [];
Data.total_GI(test) = [];
Data.Yuncor(test) = [];
Data.LinRoad(test) = [];
Data.LinOther(test) = [];
Data.Area(test) = [];
Data.Country(test) = [];
Data.Landscape(test) = [];
Data.nrspecies(test) = [];
Data.Species_loss(test) = [];
Data.SegmentBinar(test) = [];
Data.Band(test) = [];
Data.IsSegment(test) = [];
Data.CoreSize(test) = [];
end