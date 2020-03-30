function GIStatsMain(type)
clc
warning off
tic
load('GI_data.mat');
List = LisFunc;
%Radiation paper
List.Traits_Species = [1,2,5,7,9,11,13,15,17,19,23,24];
%ES paper
%List.Traits_Species_ES = [1,2,25,27,38,40,28,29,30,31,32,18,33,34,35];
%%
Input.species_option = type; % which species
Input.Linear_option = 0; % all segment elements
Input.dis_max = 8; % counter on scenarios

MakeMeanTraits(Input,List); % Creates file
CalcStatsRadi; % Creates file
namefile = ['Anova_',char(List.typeList(type)),'.mat'];
copyfile('Output_Anova.mat',namefile);
delete('Output_Anova.mat') 
namefile = ['TraitsFile_',char(List.typeList(type)),'.mat'];
copyfile('Radiation_proportions.mat',namefile);
delete('Radiation_proportions.mat')

PredictFits; % this creates Predictions file for ES paper
namefile = ['Predictions_',char(List.typeList(type)),'.mat'];
copyfile('Predictions.mat',namefile);
delete('Predictions.mat')
display('***********    READY   **********')
end

%%
function List = LisFunc
%% Numbering traits (presence, cover)
% "Radiation paper"
% 1     : Amount of Species, Numeric, single: Presence only.
% 2     : Species Diversity Shannon,Numeric, single
% 3-4   : Grime, Categorical: translation needed for Triangle
% 5-6   : Mean Grime triangle, Numeric, 3-Axis scores
% 7-8   : Mean max dispeRsal (Tamme et al), Numeric, single
% 19-20 : Ellenberg, Numeric, multiple independent traits
% 21-22 : Habitat preference, Categorical, translation needed for Triangle
% 23-24 : Mean Habitat Preference as Triangle, Numeric, 3-Axis scores
% For Tamme et al (dispeRsal)
% 9-10  : Mean Terminal Velocity, Numeric, single
% 11-12 : Mean Seed Dry weight, , Numeric, single
% 13-14 : Dispersal syndrome for dispeRsal, Categorical
% 15-16 : Growth form for dispeRsal, Categorical
% 17-18 : Release heigh, Numeric, single
% Ecosystem services not above 
% 25    : Total Nectar Index, incl Cover, Numeric, single Cover only
% 26-27 : Amount of PollinatorVisits, Numeric, single
% 28    : Mean Carbon dry weight per leaf, Numeric, single Cover only
% 29    : Mean Dry weight per leaf, Numeric, single Cover only
% 30    : Mean Carbon Dry weight Root, Numeric, single Cover only
% 31    : Mean Dry weight Root, Numeric, single Cover only
% 32    : Mean Rooting Depth Dry weight Root, Numeric, single Cover only
% 33    : LandScape diversity Numeric, single
% 34    : Life Form Shannon Diversity, Numeric, single Cover only
% 35    : Colour Shannon diversity in # colours, Numeric, single Cover only
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