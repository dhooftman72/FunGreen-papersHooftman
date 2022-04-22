function PredictFitMain(Input,List) 
for trait = List.Traits_Species
    load('Radiation_proportions.mat');
    %load('Output_Anova.mat');
    maxtrait = List.maxtraits(trait);
    for traittype = 1:maxtrait
        load('Radiation_proportions.mat');
        %load('Output_Anova.mat');
        traittext =  Input.(genvarname([char(List.trait_txts(trait))]))(traittype,2);
        if exist('Predictions.mat') ~= 0
            load('Predictions.mat')
        end
 %% Prediction loops     
        % Note that job are not only used for submitting paralel tasks but
        % also to be able to 1) time-out loops (hence also difficult stint
        % construction) and 2) when an error occurs it will just go on to
        % the next tasks. The non-run is dealt with in the collating runs
        % since when the swap file does not exist the loop is ignored.
        
        % Loops. Means all transfer parameters to keep should be in
        % the seperate swap files
        if List.loop_max_Predict <= 25 % Test mode with screen output of errors
            for loop = 1:1:List.loop_max_Predict
                clc
                display('Running Predictive runs')
                display(List.Outputs(trait))
                 display(traittext)
                display(loop)
                PredictLoop(loop,List,PropRadi,trait,traittype); % test if needed
            end
        else % full run mode, errors don't show
            display('Running Prediction statistics')
            display('Running loops') 
            stints_total = ceil(List.loop_max_Predict./100);
            % so each number not dividedable by ten will be scaled
            % upwards, for that the actual done loops are counted (see
            % below and used as max loops value below
            for stint = 1:1: stints_total
                job = createJob('configuration','Full');
                for withinloop = 1:100
                    loop = ((stint-1).*100)+ withinloop;
                    clc
                    display('Running Predictive runs')
                    display(List.Outputs(trait))
                     display(traittext)
                    display(loop)
                    createTask(job, @PredictLoop, 0,{loop,List,PropRadi,trait,traittype});
                    %Loops. Means all transfer parameters to keep should be in
                    %the seperate swap files
                end
                submit(job);
                time_out_time = 600;
                waitForState(job, 'finished',time_out_time);
                destroy(job)
            end
        end
        Loop_actually_done = loop; % because of the stints of 10 each
        pause(5) % allowing the computer to update the registery
 %% Collatie all runs
        count = 0;
        losefocals = find(Input.SegmentRecords.SegmentID<1000);
        SegmentID = Input.SegmentRecords.SegmentID(losefocals); %#ok<*FNDSB>
        nrsegments = length(SegmentID);
        DistanceScores = zeros(nrsegments,5);
        nrHits_log = zeros(nrsegments,1);
        DataScores = zeros(nrsegments,3);

        cd 'folder_out' % to temperary folder
        display('Collating loops')
        for loop = 1:1:Loop_actually_done
            name_file = ['topredict','_',mat2str(trait),'_',mat2str(traittype),'_',mat2str(loop),'.mat'];
            ignore = 1;
            test = exist(name_file); %#ok<*EXIST>
            if test~= 0
                ignore = 0; % File exists occured
            end
            clear test
            if ignore ~= 1
                load(name_file)
                count = count + 1;
                %delete so next reg convergence can be tested on existence
                
                
                %% means per Segment
                Segments = PredictRun.Ynumbers(:,1);
                x_scores = PredictRun.Xnumbers(:,2:6);
                y_scores = PredictRun.Ynumbers(:,2:4);
                for i = 1:1:length(Segments)
                    Segment = Segments(i);
                    num = find(SegmentID == Segment);
                    DistanceScores(num,:) = x_scores(i,:);
                    DataScores(num,:) =  DataScores(num,:)+ y_scores(i,:);
                    nrHits_log(num) = nrHits_log(num) + 1;
                    clear Segment num
                end
                clear i x_scores y_scores Segments
                                
                %% stats means
                correlation.Marginal.Rho(count) = PredictRun.correlation.Marginal.Rho;
                correlation.Marginal.Pval(count) = PredictRun.correlation.Marginal.Pval;
                correlation.Marginal.Deviation(count) = PredictRun.correlation.Marginal.Deviation;              
                correlation.Data.Rho(count) = PredictRun.correlation.Data.Rho;
                correlation.Data.Pval(count) = PredictRun.correlation.Data.Pval;
                correlation.Data.Deviation(count) = PredictRun.correlation.Data.Deviation;
                parameters.Beta(count,:) =  PredictRun.parameters.Beta;
                parameters.Beta_PValue(count,:) = PredictRun.parameters.Beta_PValue;
                parameters.BIC(count) = -PredictRun.parameters.Stats.bic;
                parameter.Rsquare(count) = PredictRun.parameters.Rsquare;
                parameter.FvalueModel(count) = PredictRun.parameters.FvalueModel;
                parameter.PvalueModel(count) = PredictRun.parameters.PvalueModel;
                parameter.DFModel(count) = PredictRun.parameters.DFModel;
                sensitivitiesAll(count,:) = PredictRun.SensiChanges(:,1);
                sensitivitiesCountries(count,:,:) = PredictRun.SensiChanges(:,2:4);        
            end % Ignore
            delete(name_file)
        end % Collating runs
        cd .. % back to main folder
 %% Summarise alll output
        if count == 0
            correlation.Marginal.Rho(1) = 0;
            correlation.Marginal.Pval(1) = 0;
            correlation.Marginal.Deviation(1) = 0;
            correlation.Data.Rho(1) = 0;
            correlation.Data.Pval(1) = 0;
            correlation.Data.Deviation(1) = 0;
            parameters.Beta(1,1:6) =  0;
            parameters.Beta_PValue(1,1:6) =  0;
            parameters.BIC(1) = 0;
            parameter.Rsquare(1) = 0;
            parameter.FvalueModel(1) = 0;
            parameter.PvalueModel(1) = 0;
            parameter.DFModel(1) = 0;
            sensitivitiesAll(1,1:5) = 0;
            sensitivitiesCountries(1,1:5,:) = 0;
        end
        Ind_yValues(:,1) = bsxfun(@rdivide, DataScores(:,1),nrHits_log);
        Ind_yValues(:,2) = bsxfun(@rdivide, DataScores(:,2),nrHits_log);
        Ind_yValues(:,3) = bsxfun(@rdivide, DataScores(:,3),nrHits_log);
        traitnames = {'Distance','Total_GI','Road_GI','OtherLinear_GI','AreaGI','DataY','MarginalY','PredictedLogaritmicY'};
        Segments_overview =  dataset((SegmentID),'VarNames', 'Segments');
        for i = 1:5
            Segments_overview.(genvarname([char(traitnames(i))])) = DistanceScores(:,i);
        end
        for i = 6:8
            Segments_overview.(genvarname([char(traitnames(i))])) = Ind_yValues(:,(i-5));
        end
        clear DataScores DistanceScores Ind_yValues
        
        start = [{'Logaritmic_marginal'};{'Logaritmic_data'}];
        correlations =  dataset((start),'VarNames', 'Type');
        correlations.Rho_mean = [nanmean(correlation.Marginal.Rho);nanmean(correlation.Data.Rho)];
        correlations.Rho_std = [nanstd(correlation.Marginal.Rho);nanstd(correlation.Data.Rho)];
        correlations.Pval = [nanmean(correlation.Marginal.Pval);nanmean(correlation.Data.Pval)];
        correlations.Deviation = [nanmean(correlation.Marginal.Deviation);nanmean(correlation.Data.Deviation)];
        correlations.Deviation_std = [nanstd(correlation.Marginal.Deviation);nanstd(correlation.Data.Deviation)];
        
        start = [(genvarname([char(List.Outputs(trait))]))];
        Parameters =  dataset({start},'VarNames', 'Trait');
        Parameters.R2_Model_mean = [nanmean(parameter.Rsquare)];
        Parameters.F_Model_mean = [nanmean(parameter.FvalueModel)];
        Parameters.P_Model_mean = [nanmean(parameter.PvalueModel)];
        Parameters.DF_Model_mean = [nanmean(parameter.DFModel)];
        Parameters.R2_Model_std = [nanstd(parameter.Rsquare)];
        Parameters.F_Model_std = [nanstd(parameter.FvalueModel)];
        Parameters.P_Model_std = [nanstd(parameter.PvalueModel)];
        Parameters.DF_Model_std = [nanstd(parameter.DFModel)];
        
        Parameters.BIC_mean = [nanmean(parameters.BIC)]; %#ok<*NBRAK>
        Parameters.BIC_std = [nanstd(parameters.BIC)];
        Parameters.B1_mean = [nanmean(parameters.Beta(:,1))];
        Parameters.B2_mean = [nanmean(parameters.Beta(:,2))];
        Parameters.B3_mean = [nanmean(parameters.Beta(:,3))];
        Parameters.B4_mean = [nanmean(parameters.Beta(:,4))];
        Parameters.B5_mean = [nanmean(parameters.Beta(:,5))];
        Parameters.B6_mean = [nanmean(parameters.Beta(:,6))];
        Parameters.B1_std = [nanstd(parameters.Beta(:,1))];
        Parameters.B2_std = [nanstd(parameters.Beta(:,2))];
        Parameters.B3_std = [nanstd(parameters.Beta(:,3))];
        Parameters.B4_std = [nanstd(parameters.Beta(:,4))];
        Parameters.B5_std = [nanstd(parameters.Beta(:,5))];
        Parameters.B6_std = [nanstd(parameters.Beta(:,6))];
        Parameters.B1_P = [nanmean(parameters.Beta_PValue(:,1))];
        Parameters.B2_P = [nanmean(parameters.Beta_PValue(:,2))];
        Parameters.B3_P = [nanmean(parameters.Beta_PValue(:,3))];
        Parameters.B4_P = [nanmean(parameters.Beta_PValue(:,4))];
        Parameters.B5_P = [nanmean(parameters.Beta_PValue(:,5))];
        Parameters.B6_P = [nanmean(parameters.Beta_PValue(:,6))];
        Parameters.B1_P_std = [nanstd(parameters.Beta_PValue(:,1))];
        Parameters.B2_P_std = [nanstd(parameters.Beta_PValue(:,2))];
        Parameters.B3_P_std = [nanstd(parameters.Beta_PValue(:,3))];
        Parameters.B4_P_std = [nanstd(parameters.Beta_PValue(:,4))];
        Parameters.B5_P_std = [nanstd(parameters.Beta_PValue(:,5))];
        Parameters.B6_P_std = [nanstd(parameters.Beta_PValue(:,6))];      
        clear parameters start
        
        
        if trait == 1 % species
            base = [0,0,0,0,0];
        else
            base = [Predictions.Species.Sensitivities.All.Distance,...
                Predictions.Species.Sensitivities.All.GI_Total,...
                Predictions.Species.Sensitivities.All.GI_Road,...
                Predictions.Species.Sensitivities.All.GI_OtherLinear,...
                Predictions.Species.Sensitivities.All.GI_Area];
        end
        nameStrucArray = {'Distance','GI_Total','GI_Road','GI_OtherLinear','GI_Area',...
                            'DistancePO','GI_TotalPO','GI_RoadPO','GI_OtherLinearPO','GI_AreaPO'};
       %Initiate, so first collumns are the means
        SensitivitiesAll = dataset(0,'VarNames',nameStrucArray(1));
        SensitivitiesCountries = dataset(0,'VarNames','DistanceP1_2');
        for factor = 2:5
        SensitivitiesAll.(genvarname(char(nameStrucArray(factor))))  = 0;
        end
        for factor = 1:5
            [SensitivitiesAll,SensitivitiesCountries] = SensColate(SensitivitiesAll,SensitivitiesCountries,sensitivitiesAll,sensitivitiesCountries,base,nameStrucArray,factor);
        end
        clear sensitivities
        
        Predictions.(genvarname([char(List.Outputs(trait))])).Sensitivities.All(traittype,:) = SensitivitiesAll;
        Predictions.(genvarname([char(List.Outputs(trait))])).Sensitivities.Countries(traittype,:) = SensitivitiesCountries;
        Predictions.(genvarname([char(List.Outputs(trait))])).(genvarname([char(traittext)])).Correlations = correlations;
        Predictions.(genvarname([char(List.Outputs(trait))])).(genvarname([char(traittext)])).Parameters = Parameters;
        Predictions.(genvarname([char(List.Outputs(trait))])).(genvarname([char(traittext)])).Segment_overview  = Segments_overview;
        Predictions.(genvarname([char(List.Outputs(trait))])).(genvarname([char(traittext)])).TimesConverged = [(count./Loop_actually_done)];
        
        % Clean sheet approach, save only what is needed and remove all remnant
        % variables
        save('Predictions','Predictions')
        save('trait_transfer','trait','traittype','List')
        clear Predictions
        clear all
        load('trait_transfer.mat');
        delete('trait_transfer.mat');
    end %traittype
end % traits
end % function

function [SensitivitiesAll,SensitivitiesCountries] = SensColate(SensitivitiesAll,SensitivitiesCountries,InVarAll,InVarCountries,base,names,factor)
nameStruc = names(factor);
sens = InVarAll(:,factor);
sens(sens<-2) = NaN;
sens(sens>2) = NaN;
SensitivitiesAll.(genvarname(char(nameStruc))) = nanmean(sens);
SensitivitiesAll.(genvarname([(char(nameStruc)),'Std'])) = nanstd(sens);
[~,p,ci] = ttest(sens(isnan(sens)~=1), base(factor));
SensitivitiesAll.(genvarname([(char(nameStruc)),'PO'])) = p;
SensitivitiesAll.(genvarname([(char(nameStruc)),'ActualRuns'])) = length(sens(isnan(sens)~=1));

[pCtry,psCtry,MeansCtry,StdsCtry] = difTest(InVarCountries(:,factor,:));
SensitivitiesAll.(genvarname([(char(nameStruc)),'RegionPdif'])) = pCtry;
SensitivitiesCountries.(genvarname([(char(nameStruc)),'P1_2'])) = psCtry(1);
SensitivitiesCountries.(genvarname([(char(nameStruc)),'P1_3'])) = psCtry(2);
SensitivitiesCountries.(genvarname([(char(nameStruc)),'P2_3'])) = psCtry(3);
SensitivitiesCountries.(genvarname([(char(nameStruc)),'Mean_1'])) = MeansCtry(1);
SensitivitiesCountries.(genvarname([(char(nameStruc)),'Mean_2'])) = MeansCtry(2);
SensitivitiesCountries.(genvarname([(char(nameStruc)),'Mean_3'])) = MeansCtry(3);
SensitivitiesCountries.(genvarname([(char(nameStruc)),'Std_1'])) = StdsCtry(1);
SensitivitiesCountries.(genvarname([(char(nameStruc)),'Std_2'])) = StdsCtry(2);
SensitivitiesCountries.(genvarname([(char(nameStruc)),'Std_3'])) = StdsCtry(3);
end

function [p,ps,Means,Stds] = difTest(InVar)
InVar(InVar<-2) = NaN;
InVar(InVar>2) = NaN;
[~,ps(1),~] = ttest(InVar(:,1), InVar(:,2));
[~,ps(2),~] = ttest(InVar(:,1), InVar(:,3));
[~,ps(3),~] = ttest(InVar(:,2), InVar(:,3));
p = min(ps);
Means =  nanmean(InVar);
Stds = nanstd(InVar);
end