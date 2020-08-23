function GIStatsMain(type)
clc
warning off
tic
cd('TheData')
load('GI_data.mat');
cd ..
List = LisFunc;
%Radiation paper
List.Traits_Species = [1,5,7,9,11,13,15,17,19,23,24];
%ES paper
%List.Traits_Species_ES = [1,2,25,27,38,40,28,29,30,31,32,18,33,34,35];
%%
List.species_option = type; % which species; 
%Type 1 = single core; Type 2= Any core; Type 3 = all species; 
% Type 4 = not present single target; Type 5: all incoming species
List.Linear_option = 0; % all segment elements
List.dis_max = 8; % counter on scenarios
List.loop_max_Predict = 300; % Prediction Loops
List.sense_factor = 0.25; % sensitivity factor
if exist('Radiation_proportions.mat') == 0
MakeMeanTraits(Input,List); % Creates outputfile: Radiation_proportions.mat
else
     testr = testRadi(type);
     if testr == 1
         MakeMeanTraits(Input,List); % Creates outputfile: Radiation_proportions.mat
     end
end
CalcStatsRadi; % Creates output file: Output_Anova.mat
namefile = ['Anova_',char(List.typeList(type)),'.mat'];
copyfile('Output_Anova.mat',namefile);
namefile = ['TraitsFile_',char(List.typeList(type)),'.mat'];
copyfile('Radiation_proportions.mat',namefile);

%% ES Paper
if exist('folder_out','dir') == 0
    mkdir('folder_out')
end
PredictFitMain(Input,List); % Creates output file:Predictions.mat
namefile = ['Predictions_',char(List.typeList(type)),'.mat'];
copyfile('Predictions.mat',namefile);
delete('Output_Anova.mat') 
delete('Radiation_proportions.mat')
delete('Predictions.mat')
clc
display('  ')
display('  ')
display('***********    READY   **********')
end

function testr = testRadi(type)
load('Radiation_proportions.mat')
testr = 0;
if List.species_option ~= type;
    testr = 1;
end
end
%%
function List = LisFunc
%% Numbering traits (presence, cover)
% 1     : Amount of Species, Numeric, single: Presence only. (Both papers)
% 2     : Species Diversity Shannon,Numeric, single (ES paper)

% "Radiation Paper" traits
% 3-4   : Grime, Categorical: translation needed for Triangle (abundance/cover)
% 5-6   : Mean Grime triangle, Numeric, 3-Axis scores (abundance/cover)
% 7-8   : Mean max dispeRsal (Tamme et al), Numeric, single (abundance/cover)
% 19-20 : Ellenberg, Numeric, multiple independent traits (abundance/cover)
% 21-22 : Habitat preference, Categorical, translation needed for Triangle (abundance/cover)
% 23-24 : Mean Habitat Preference as Triangle, Numeric, 3-Axis scores (abundance/cover)
% For Tamme et al (dispeRsal)
% 9-10  : Mean Terminal Velocity, Numeric (abundance/cover)
% 11-12 : Mean Seed Dry weight, , Numeric (abundance/cover)
% 13-14 : Dispersal syndrome for dispeRsal, Categorical (abundance/cover)
% 15-16 : Growth form for dispeRsal, Categorical (abundance/cover)
% 17-18 : Release heigh, Numeric(also ES!!) (abundance/cover)

% "Ecosystem services paper traits, not above" 
% All Numeric & Cover corrected
% 25    : Total Nectar Index, incl Cover 
% 27    : Amount of PollinatorVisits
% 28    : Mean Carbon dry weight per leaf
% 29    : Mean Dry weight per leaf
% 30    : Mean Carbon Dry weight Root
% 31    : Mean Dry weight Root, Numeric
% 32    : Mean Rooting Depth Dry weight Root
% 33    : LandScape diversity (no plant information!!)
% 34    : Life Form Shannon Diversity
% 35    : Colour Shannon diversity in # colours
%% LIST to transfer
List.trait_txts = {'Onetraittxt','Onetraittxt','Grimetxt','Grimetxt','GrimeCentroidtxt',... %5
                    'GrimeCentroidtxt','Onetraittxt','Onetraittxt','Onetraittxt','Onetraittxt',... %10
                    'Onetraittxt','Onetraittxt','DStxt','DStxt','GFtxt',... %15
                    'GFtxt','Onetraittxt','Onetraittxt','Ellenbergtxt','Ellenbergtxt',... %20
                    'Habitattxt_short','Habitattxt_short','HabitatCentroidtxt','HabitatCentroidtxt','Onetraittxt',...%25
                    'Onetraittxt','Onetraittxt','Onetraittxt','Onetraittxt','Onetraittxt',... % 30
                    'Onetraittxt','Onetraittxt','Onetraittxt','Onetraittxt','Onetraittxt','Onetraittxt'}; % 35
List.Outputs = {'Species','SpeciesShannon','Grime','Grime_cover','Grimecentroid',... %5
                 'Grimecentroid_cover','DispeRsal', 'DispeRsal_cover','TerminalVelocity','TerminalVelocity_cover',... %10
                 'DryWeight','Dryweight_cover','DispSyndR','DispSyndR_cover','GF',... %15
                 'GF_cover','Release_Height','Release_height_Cover','Ellenberg','Ellenberg_cover',... %20
                 'Habitat','Habitat_cover','HabitatCentroid','HabitatCentroid_cover','NectarIndex',... %25
                 'PollinatorVisits','PollinatorVisits_cover','CarbonLeaf','DryWeightLeaf','CarbonRoot',...%30
                 'DryWeightRoot','RootingDepth','Landscape_Diversity','Life_form_Diversity','ColourShannon'};%35
List.CatagoricalTraits = [3,4,13,14,15,16,21,22];
List.traits_to_do = length(List.Outputs);    
List.typeList = [{'single_core'},{'All_Cores'},{'All_species'},{'Incoming_SingleCore'},{'Incoming_allCores'}];
end