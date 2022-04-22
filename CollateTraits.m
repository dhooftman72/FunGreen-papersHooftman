function [SelectRadi,Traitmaxs] = CollateTraits(Input,List,SpeciesSegment,nrspeciesSegment,segmentnr,segments_to_do,SelectRadi)
%% number of species
SelectRadi.nrTargetspeciessegment(segmentnr,1) = nrspeciesSegment;
Trans.segmentnr = segmentnr; %Both will do in this function
Trans.nrspeciesSegment = nrspeciesSegment; %Both will do in this function

%% SPECIES calculated on dispeRsal array
Trans.Traitmax = 1;
Trans.maxi = 1;
Trans.tf = [1,3,4]; % Species list; trait list; strength nr
Trans.prop = 7; % of dispeRsal
Trans.cover = 8; % of dispeRsal
Trans.which = 'Both';
Trans.naming = 'SeedDispersalR'; %dispeRsal array
Traitmaxs(Trans.maxi) = Trans.Traitmax;
[Selection,~] = StandNum(Input,List,Trans,SpeciesSegment,SelectRadi);
SelectRadi.Species.value(Trans.segmentnr,1) =  Selection.nrspecs;
SelectRadi.Species.nrspecs(Trans.segmentnr,1) =  Selection.nrspecs;
Trans.maxi = 36;
SelectRadi.LogSpecies.value(Trans.segmentnr,1) =  log10(Selection.nrspecs+1);
SelectRadi.LogSpecies.nrspecs(Trans.segmentnr,1) =  log10(Selection.nrspecs+1);
Traitmaxs(Trans.maxi) = Trans.Traitmax;

%% Shannon Diversity
Trans.Traitmax = 1;
Trans.maxi = 2;
Traitmaxs(Trans.maxi) = Trans.Traitmax;
shannon = 0;
Total_cover = Selection.CoverTrait./Selection.Total_cover;
for i = 1:1:Selection.nrspecs
        shannon = shannon + (Total_cover(i)*log(Total_cover(i)));
end
SelectRadi.SpeciesShannon.value(Trans.segmentnr,1) =  -shannon;
SelectRadi.SpeciesShannon.nrspecs(Trans.segmentnr,1) =  Selection.nrspecs;

%% GRIME Lifeform CATEGORICAL TRAIT
Trans.maxi = [3,4];
Trans.tf = [1,2,NaN]; % Species list; trait list; strength nr
Trans.prop = Trans.maxi(1);
Trans.cover = Trans.maxi(2);
Trans.which = 'Both';
Trans.naming = 'GrimeSpecies';
Trans.TraitsToDo =[1,2,3,4,5,6,7];
Trans.Traitmax = length(Trans.TraitsToDo);
Traitmaxs(Trans.maxi) = Trans.Traitmax;
[~,SelectRadi] = StandNum(Input,List,Trans,SpeciesSegment,SelectRadi);  
RichnessTrait = zeros(Trans.nrspeciesSegment,7);
CoverTrait = zeros(Trans.nrspeciesSegment,7);
count = 0;
for traittype =  Trans.TraitsToDo
    count = count + 1;
    %Create Yes/no matrices per Life form
    for x = 1:1:nrspeciesSegment
        PresList = find(Input.GrimeSpecies(:,1) == SpeciesSegment(x,1));
        %display(length(PresList))
        if isempty(PresList) ~= 1
            Form = Input.GrimeSpecies(PresList,2);
            Pres = find(Form == traittype);
            if isempty(Pres) ~= 1
                RichnessTrait(x,count) = 1;
                CoverTrait(x,count) = SpeciesSegment(x,2);
            end
        else
            RichnessTrait(x,count) = 9999;
            CoverTrait(x,count)= 9999;
        end
    end
    clear PresList Pres Forms x
end
clear  traittype
% Remove species with no data
RemoveDisp = find(RichnessTrait(:,1) == 9999);
RichnessTrait(RemoveDisp,:) = [];
CoverTrait(RemoveDisp,:) = [];
clear  RemoveDisp
%check for zero entries
listzeros = find(mean(RichnessTrait,2) == 0);
RichnessTrait(listzeros,:) = [];
CoverTrait(listzeros,:)= [];

for i = 1:1:Trans.Traitmax
    Selects(i) = sum(RichnessTrait(:,i))./ (sum(sum(RichnessTrait)));
    Selects_Cover(i) = sum(CoverTrait(:,i))./(sum(sum(CoverTrait)));
end
SelectRadi.(genvarname(char(List.Outputs(Trans.prop)))).value(segmentnr,:)  = Selects./(sum(Selects));
SelectRadi.(genvarname(char(List.Outputs(Trans.cover)))).value(segmentnr,:)  = Selects_Cover./(sum(Selects_Cover));
SelectRadi.(genvarname(char(List.Outputs(Trans.prop)))).nrspecs(segmentnr,:) =  NaN;
SelectRadi.(genvarname(char(List.Outputs(Trans.cover)))).nrspecs(segmentnr,:) =  NaN;
clear RichnessTrait CoverTrait StrengthTrait TotalCoverSegment disp i x
clear Traitmax Selects Selects_Cover nrspecs 
%% GRIME Lifeform Version 2 MUMERICAL TRAIT
Trans.Traitmax = 3;
Trans.maxi = [5,6];
Trans.prop = Trans.maxi(1);
Trans.cover = Trans.maxi(2);
Trans.which = 'Both';
Traitmaxs(Trans.maxi) = Trans.Traitmax;
RichnessTrait = SelectRadi.(genvarname(char(List.Outputs(3)))).value(Trans.segmentnr,:);
CoverTrait = SelectRadi.(genvarname(char(List.Outputs(4)))).value(Trans.segmentnr,:);

S = Input.Grimetranslation(:,1);
C = Input.Grimetranslation(:,2);
R = Input.Grimetranslation(:,3);

Selects(1) = sum(bsxfun(@times,RichnessTrait,C'));
Selects(2) = sum(bsxfun(@times,RichnessTrait,S'));
Selects(3) = sum(bsxfun(@times,RichnessTrait,R'));
SelectRadi.(genvarname(char(List.Outputs(Trans.prop)))).value(segmentnr,:)  = Selects;
SelectRadi.(genvarname(char(List.Outputs(Trans.prop)))).nrspecs(segmentnr,:) = [Selection.nrspecs, Selection.nrspecs, Selection.nrspecs];

Selects_Cover(1) = sum(bsxfun(@times,CoverTrait,C'));
Selects_Cover(2) = sum(bsxfun(@times,CoverTrait,S'));
Selects_Cover(3) = sum(bsxfun(@times,CoverTrait,R'));
SelectRadi.(genvarname(char(List.Outputs(Trans.cover)))).value(segmentnr,:)  =  Selects_Cover;
SelectRadi.(genvarname(char(List.Outputs(Trans.cover)))).nrspecs(segmentnr,:) = [Selection.nrspecs, Selection.nrspecs, Selection.nrspecs];
clear RichnessTrait CoverTrait StrengthTrait TotalCoverSegment disp i x
clear Traitmax Selects Selects_Cover nrspecs 
%% DispeRsal tamme et al. NUMERIC TRAIT
Trans.Traitmax = 1;
Trans.maxi = [7,8];
Trans.tf = [1,3,4]; % Species list; trait list; strength nr
Trans.prop = Trans.maxi(1);
Trans.cover = Trans.maxi(2);
Traitmaxs(Trans.maxi) = Trans.Traitmax;
Trans.which = 'Both';
Trans.naming = 'SeedDispersalR';
RichnessTrait = zeros(nrspeciesSegment,Trans.Traitmax);
CoverTrait = zeros(nrspeciesSegment,Trans.Traitmax);
for x = 1:1:nrspeciesSegment
    PresList = find(Input.SeedDispersalR(:,1) == SpeciesSegment(x,1));
    total_cover = 0;
    if isempty(PresList) ~= 1
        RichnessTrait_T = (Input.SeedDispersalR(PresList,3));
        RichnessTrait_T(RichnessTrait_T<0) = 0;
        Strenght = Input.SeedDispersalR(PresList,4);
        RichnessTrait(x,:) = ((sum((bsxfun(@times, (RichnessTrait_T),Strenght))))./sum(Strenght));
        CoverTrait(x) = SpeciesSegment(x,2);
    else
        RichnessTrait(x,:) = 9999;
        CoverTrait(x) = 9999;
    end
end
clear PresList x
%Remove species without Value
RemoveLife = find(RichnessTrait(:,1) == 9999);
RichnessTrait(RemoveLife,:) = [];
CoverTrait(RemoveLife) = [];
clear RemoveLife
to_take = find(isnan(RichnessTrait) == 0); % So only existing values are taken!!!
RichnessTrait =  RichnessTrait(to_take);
RichnessTrait  = reshape(RichnessTrait ,length(RichnessTrait),1);
CoverTrait =  CoverTrait(to_take);
CoverTrait = reshape(CoverTrait,length(CoverTrait),1);
Nr_covers =  (bsxfun(@times, RichnessTrait,CoverTrait));
nrspecs = length(RichnessTrait);
SelectRadi.(genvarname(char(List.Outputs(Trans.prop)))).value(segmentnr,1)  = 10^(mean(RichnessTrait));
SelectRadi.(genvarname(char(List.Outputs(Trans.cover)))).value(segmentnr,1)  =  sum(Nr_covers)./sum(CoverTrait);
SelectRadi.(genvarname(char(List.Outputs(Trans.prop)))).nrspecs(segmentnr,1) =  nrspecs;
SelectRadi.(genvarname(char(List.Outputs(Trans.cover)))).nrspecs(segmentnr,1) =  nrspecs;

%% Terminal Velocity NUMERIC TRAIT
Trans.Traitmax = 1;
Trans.maxi = [9,10];
Trans.tf = [1,2,NaN]; % Species list; trait list; strength nr
Trans.prop = Trans.maxi(1);
Trans.cover = Trans.maxi(2);
Traitmaxs(Trans.maxi) = Trans.Traitmax;
Trans.which = 'Both';
Trans.naming = 'TerminalVelocitySpecies';
[~,SelectRadi] = StandNum(Input,List,Trans,SpeciesSegment,SelectRadi);  

%% Seed Dry Weight  NUMERIC TRAIT
Trans.Traitmax = 1;
Trans.maxi = [11,12];
Trans.tf = [1,2,NaN]; % Species list; trait list; strength nr
Trans.prop = Trans.maxi(1);
Trans.cover = Trans.maxi(2);
Traitmaxs(Trans.maxi) = Trans.Traitmax;
Trans.which = 'Both';
Trans.naming = 'SeedDryWeightSpecies';
[~,SelectRadi] = StandNum(Input,List,Trans,SpeciesSegment,SelectRadi); 

%% Dispersal type for DispeRsal CATEGORICAL TRAIT
Trans.maxi = [13,14];
Trans.tf = [1,2,NaN]; % Species list; trait list; strength nr
Trans.prop = Trans.maxi(1);
Trans.cover = Trans.maxi(2);
Trans.naming = 'SeedDispersalR';
Trans.TraitsToDo =[4,5,6,7,8];
Trans.which = 'Both';
Trans.Traitmax = length(Trans.TraitsToDo);
Traitmaxs(Trans.maxi) = Trans.Traitmax;
[~,SelectRadi] = StandNum(Input,List,Trans,SpeciesSegment,SelectRadi);  

%% GF type for DispeRsal CATEGORICAL TRAIT
Trans.maxi = [15,16];
Trans.tf = [1,2,NaN]; % Species list; trait list; strength nr
Trans.prop = Trans.maxi(1);
Trans.cover = Trans.maxi(2);
Trans.naming = 'GFType';
Trans.TraitsToDo =[1,2,3];
Trans.which = 'Both';
Trans.Traitmax = length(Trans.TraitsToDo);
Traitmaxs(Trans.maxi) = Trans.Traitmax;
[~,SelectRadi] = StandNum(Input,List,Trans,SpeciesSegment,SelectRadi); 

%% Plant/RELEASE Height, NUMERIC TRAIT
Trans.Traitmax = 1;
Trans.maxi = [17,18];
Trans.tf = [1,2,NaN]; % Species list; trait list; strength nr
Trans.prop = Trans.maxi(1);
Trans.cover = Trans.maxi(2);
Traitmaxs(Trans.maxi) = Trans.Traitmax;
Trans.which = 'Both';
Trans.naming = 'PlantHeightSpecies';
[~,SelectRadi] = StandNum(Input,List,Trans,SpeciesSegment,SelectRadi); 

%% Ellenbergs NUMERIC TRAIT; 7 in one go
Trans.Traitmax = 7;
Trans.maxi = [19,20];
Trans.prop = Trans.maxi(1);
Trans.cover = Trans.maxi(2);
Traitmaxs(Trans.maxi) = Trans.Traitmax;
Trans.which = 'Both';
RichnessTrait = zeros(nrspeciesSegment,Trans.Traitmax);
CoverTrait = zeros(nrspeciesSegment,Trans.Traitmax);
for x = 1:1:nrspeciesSegment
    PresList = find(Input.EllenbergSpecies(:,1) == SpeciesSegment(x,1));
    if isempty(PresList) ~= 1
        RichnessTrait(x,:) = Input.EllenbergSpecies(PresList,2:(Trans.Traitmax+1));
        CoverTrait(x) = SpeciesSegment(x,2);
    else
        RichnessTrait(x,:) = 9999;
        CoverTrait(x) = 9999;
    end
end
clear PresList x

%Remove species without Ellenberg
RemoveLife = find(RichnessTrait(:,1) == 9999);
RichnessTrait(RemoveLife,:) = [];
CoverTrait(RemoveLife) = [];
clear RemoveLife
for i = 1:1:Trans.Traitmax
    Nr_trait = RichnessTrait(:,i);
    to_take = find(isnan(Nr_trait) == 0); % So only existing values are taken!!!
    Nr_trait_tt =  Nr_trait(to_take);
    Nr_trait_tt  = reshape( Nr_trait_tt ,length( Nr_trait_tt),1);
    Cover_tt =  CoverTrait(to_take);
    Cover_tt = reshape(Cover_tt,length(Cover_tt),1);
    
    MeanTrait(1,i) = mean(Nr_trait_tt);
    Nr_covers =  (bsxfun(@times, Nr_trait_tt,Cover_tt));
    MeanCover(1,i) = sum(Nr_covers)./sum(Cover_tt);
    nrspecs(1,i)= length(Nr_trait_tt);
    clear to_take
end
SelectRadi.(genvarname(char(List.Outputs(Trans.prop)))).value(segmentnr,:) = MeanTrait(:) ;
SelectRadi.(genvarname(char(List.Outputs(Trans.cover)))).value(segmentnr,:) =  MeanCover(:);
SelectRadi.(genvarname(char(List.Outputs(Trans.prop)))).nrspecs(segmentnr,:) =  nrspecs(:);
SelectRadi.(genvarname(char(List.Outputs(Trans.cover)))).nrspecs(segmentnr,:) =  nrspecs(:);
clear RichnessTrait CoverTrait StrengthTrait TotalCoverSegment disp i x
clear Traitmax Selects Selects_Cover nrspecs 
%% Habitat type CATEGORICAL TRAIT 
Trans.Traitmax = 6; 
Trans.maxi = [21,22];
Trans.prop = Trans.maxi(1);
Trans.cover = Trans.maxi(2);
Traitmaxs(Trans.maxi) = Trans.Traitmax;
Trans.which = 'Both';
RichnessTrait = zeros(nrspeciesSegment,17);
CoverTrait = zeros(nrspeciesSegment,17);
HabitatTable = zeros(nrspeciesSegment,6);
HabitatCover = zeros(nrspeciesSegment,6);
count = 0;
for traittype = 1:18
    count = count + 1;
    %Create Yes/no matrices per Life form
    for x = 1:1:nrspeciesSegment
        PresList = find(Input.HabitatSpecies(:,1) == SpeciesSegment(x,1));
        if isempty(PresList) ~= 1
            Form = Input.HabitatSpecies(PresList,2);
            Pres = find(Form == traittype);
            if isempty(Pres) ~= 1
                RichnessTrait(x,count) = 1;
                CoverTrait(x,count) = SpeciesSegment(x,2);
            end
        else
            RichnessTrait(x,count) = 9999;
            CoverTrait(x,count)= 9999;
        end
    end
    clear PresList Pres Form x
end
clear traittype count

% Remove species with no data
RemoveClone = find(RichnessTrait(:,1) == 9999);
RichnessTrait(RemoveClone,:) = [];
CoverTrait(RemoveClone,:) = [];
HabitatTable(RemoveClone,:) = [];
HabitatCover(RemoveClone,:) = [];
clear  RemoveClone

%check for zero entries
listzeros = find(mean(RichnessTrait,2) == 0);
RichnessTrait(listzeros,:) = [];
CoverTrait(listzeros,:)= [];

% Join habitat_types
for x = 1:1:length(RichnessTrait(:,1))   
    %New_types HABITAT TYPE
    % 1 = Natural Grassland = 9,18,6,1,14
    % 2 = Woodland = 7, 8 & 5 (Linear and border)
    % 3 = Ruderal = 16
    % 4 = Arable/Improved = 12 & 3
    % 5 = Heathland = 11
    % 6 = Others = 2,4,10,13,17
    HabitatTable(x,1) = max([RichnessTrait(x,9), RichnessTrait(x,17),...
        RichnessTrait(x,6), RichnessTrait(x,1),...
        RichnessTrait(x,14)]);
    HabitatTable(x,2) = max([RichnessTrait(x,7), RichnessTrait(x,8),RichnessTrait(x,5)]);
    HabitatTable(x,3) = RichnessTrait(x,15);
    HabitatTable(x,4) = max([RichnessTrait(x,12), RichnessTrait(x,3)]);
    HabitatTable(x,5) = RichnessTrait(x,11);
    HabitatTable(x,6) =   max([RichnessTrait(x,2), RichnessTrait(x,4),...
        RichnessTrait(x,10),...
        RichnessTrait(x,13), RichnessTrait(x,16)]);
    % species counter
    HabitatCover(x,1) = max([CoverTrait(x,9), CoverTrait(x,17),...
        CoverTrait(x,6), CoverTrait(x,1),...
        CoverTrait(x,14)]);
    HabitatCover(x,2) = max([CoverTrait(x,7), CoverTrait(x,8)]);
    HabitatCover(x,3) = CoverTrait(x,15);
    HabitatCover(x,4) = max([CoverTrait(x,12), CoverTrait(x,3)]);
    HabitatCover(x,5) = CoverTrait(x,11);
    HabitatCover(x,6) =   max([CoverTrait(x,2), CoverTrait(x,4),...
        CoverTrait(x,5), CoverTrait(x,10),...
        CoverTrait(x,13), CoverTrait(x,16)]);
end

for i = 1:1:6
    if length(RichnessTrait(:,1)) >0
        Selects(i) = sum(HabitatTable(:,i))./ (sum(sum(HabitatTable)));
        Selects_Cover(i) = sum(HabitatCover(:,i))./(sum(sum(HabitatCover)));
    else
        Selects(i) =0;
        Selects_Cover(i) = 0;
    end
end

SelectRadi.(genvarname(char(List.Outputs(Trans.prop)))).value(segmentnr,:)  = Selects./(sum(Selects));
SelectRadi.(genvarname(char(List.Outputs(Trans.cover)))).value(segmentnr,:)  = Selects_Cover./(sum(Selects_Cover));
SelectRadi.(genvarname(char(List.Outputs(Trans.prop)))).nrspecs(segmentnr,:) =  SelectRadi.(genvarname(char(List.Outputs(Trans.prop)))).value(segmentnr,:);
SelectRadi.(genvarname(char(List.Outputs(Trans.cover)))).nrspecs(segmentnr,:) =  SelectRadi.(genvarname(char(List.Outputs(Trans.prop)))).value(segmentnr,:);
clear RichnessTrait CoverTrait StrengthTrait TotalCoverSegment disp i x
clear Traitmax Selects Selects_Cover HabitatTable HabitatTableCover
%% Habitat Prefernce Version 2 MUMERICAL TRAIT
Trans.Traitmax = 3;
Trans.maxi = [23,24];
Trans.prop = Trans.maxi(1);
Trans.cover = Trans.maxi(2);
Trans.which = 'Both';
Traitmaxs(Trans.maxi) = Trans.Traitmax;
RichnessTrait = SelectRadi.(genvarname(char(List.Outputs(21)))).value(segmentnr,:);
CoverTrait =  SelectRadi.(genvarname(char(List.Outputs(22)))).value(segmentnr,:);

Grass = cell2mat(Input.Habitattxt_short(:,3));
Wood = cell2mat(Input.Habitattxt_short(:,4));
Ruder = cell2mat(Input.Habitattxt_short(:,5));

Selects(1) = sum(bsxfun(@times,RichnessTrait,Grass'));
Selects(2) = sum(bsxfun(@times,RichnessTrait,Wood'));
Selects(3) = sum(bsxfun(@times,RichnessTrait,Ruder'));
Selects = Selects./sum(Selects);

SelectRadi.(genvarname(char(List.Outputs(Trans.prop)))).value(segmentnr,:)  = Selects ;
SelectRadi.(genvarname(char(List.Outputs(Trans.prop)))).nrspecs(segmentnr,:) = [SelectRadi.Species.nrspecs(segmentnr,1),...
                                                                          SelectRadi.Species.nrspecs(segmentnr,1),...
                                                                          SelectRadi.Species.nrspecs(segmentnr,1)];
Selects_Cover(1) = sum(bsxfun(@times,CoverTrait,Grass'));
Selects_Cover(2) = sum(bsxfun(@times,CoverTrait,Wood'));
Selects_Cover(3) = sum(bsxfun(@times,CoverTrait,Ruder'));
Selects_Cover = Selects_Cover./sum(Selects_Cover);
SelectRadi.(genvarname(char(List.Outputs(Trans.cover)))).value(segmentnr,:)  = Selects_Cover;
SelectRadi.(genvarname(char(List.Outputs(Trans.cover)))).nrspecs(segmentnr,:) = [SelectRadi.Species.nrspecs(segmentnr,1),...
                                                                          SelectRadi.Species.nrspecs(segmentnr,1),...
                                                                          SelectRadi.Species.nrspecs(segmentnr,1)];
                                                                      
clear RichnessTrait CoverTrait StrengthTrait TotalCoverSegment disp i x
clear Traitmax Selects Selects_Cover  nrspecs 
%% Nectar INDEX NUMERICAL TRAIT
Trans.Traitmax = 1;
Trans.maxi = 25;
Trans.prop = Trans.maxi(1);
Trans.which = 'Cover';
Traitmaxs(Trans.maxi) = Trans.Traitmax;
RichnessTrait = zeros(nrspeciesSegment,Trans.Traitmax);
CoverTrait = zeros(nrspeciesSegment,Trans.Traitmax);
for x = 1:1:nrspeciesSegment
    PresList = find(Input.Nectar(:,1) == SpeciesSegment(x,1));
    total_cover = 0;
    if isempty(PresList) ~= 1
        Richtmp(1) = (Input.Nectar(PresList,3)); % Duration
        Richtmp(2) = (Input.Nectar(PresList,4)); % Nectar amounts
        CoverTrait(x) = SpeciesSegment(x,2);
        if Richtmp(1) == 1 || Richtmp(1) == 2 % Relatively low flowering time 
            if  Richtmp(2) == 0 % No Nectar
                RichnessTrait(x) = 0;
            elseif  Richtmp(2) == 1 % Low Nectar production
                if CoverTrait(x) <= 0.02
                     RichnessTrait(x) = 1;
                elseif CoverTrait(x) <= 0.25
                     RichnessTrait(x) = 1;
                elseif CoverTrait(x) == 0.5
                     RichnessTrait(x) = 2;
                end
             elseif  Richtmp(2) == 2 % Medium Nectar production
                if CoverTrait(x) <= 0.02
                     RichnessTrait(x) = 1;
                elseif CoverTrait(x) <= 0.25
                     RichnessTrait(x) = 2;
                elseif CoverTrait(x) == 0.5
                     RichnessTrait(x) = 3;
                end
            elseif  Richtmp(2) == 3 % High Nectar production
                if CoverTrait(x) <= 0.02
                     RichnessTrait(x) = 2;
                elseif CoverTrait(x) <= 0.25
                     RichnessTrait(x) = 3;
                elseif CoverTrait(x) == 0.5
                     RichnessTrait(x) = 4;
                end
            elseif  Richtmp(2) == 9999 % Unknown Nectar production
                RichnessTrait(x) = 9999;
            end

        elseif  Richtmp(1) == 3 || Richtmp(1) == 4 || Richtmp(1) == 9999 % Relatively Medium flowering time and unknown
            if  Richtmp(2) == 0 % No Nectar
                RichnessTrait(x) = 0;
            elseif  Richtmp(2) == 1 % Low Nectar production
                if CoverTrait(x) <= 0.02
                     RichnessTrait(x) = 1;
                elseif CoverTrait(x) <= 0.25
                     RichnessTrait(x) = 2;
                elseif CoverTrait(x) == 0.5
                     RichnessTrait(x) = 3;
                end
             elseif  Richtmp(2) == 2 % Medium Nectar production
                if CoverTrait(x) <= 0.02
                     RichnessTrait(x) = 2;
                elseif CoverTrait(x) <= 0.25
                     RichnessTrait(x) = 3;
                elseif CoverTrait(x) == 0.5
                     RichnessTrait(x) = 4;
                end
            elseif  Richtmp(2) == 3 % High Nectar production
                if CoverTrait(x) <= 0.02
                     RichnessTrait(x) = 3;
                elseif CoverTrait(x) <= 0.25
                     RichnessTrait(x) = 4;
                elseif CoverTrait(x) == 0.5
                     RichnessTrait(x) = 5;
            elseif  Richtmp(2) == 9999 % Unknown Nectar production
                RichnessTrait(x) = 9999;
                end
            end
        elseif  Richtmp(1)> 4 && Richtmp(1) ~= 9999 % Relatively high flowering time 
        if  Richtmp(2) == 0 % No Nectar
                RichnessTrait(x) = 0;
            elseif  Richtmp(2) == 1 % Low Nectar production
                if CoverTrait(x) <= 0.02
                     RichnessTrait(x) = 2;
                elseif CoverTrait(x) <= 0.25
                     RichnessTrait(x) = 3;
                elseif CoverTrait(x) == 0.5
                     RichnessTrait(x) = 4;
                end
             elseif  Richtmp(2) == 2 % Medium Nectar production
                if CoverTrait(x) <= 0.02
                     RichnessTrait(x) = 3;
                elseif CoverTrait(x) <= 0.25
                     RichnessTrait(x) = 4;
                elseif CoverTrait(x) == 0.5
                     RichnessTrait(x) = 5;
                end
            elseif  Richtmp(2) == 3 % High Nectar production
                if CoverTrait(x) <= 0.02
                     RichnessTrait(x) = 4;
                elseif CoverTrait(x) <= 0.25
                     RichnessTrait(x) = 5;
                elseif CoverTrait(x) == 0.5
                     RichnessTrait(x) = 5;
            elseif  Richtmp(2) == 9999 % Unknown Nectar production
                RichnessTrait(x) = 9999;
                end
            end
        end
    else
        RichnessTrait(x) = 9999;
        CoverTrait(x) = 9999;
    end
end
clear PresList x CoverTrait
%Remove species without Value
RemoveLife = find(RichnessTrait(:) == 9999);
RichnessTrait(RemoveLife) = [];
to_take = find(isnan(RichnessTrait) == 0); % So only existing values are taken!!!
RichnessTrait =  RichnessTrait(to_take);
RichnessTrait  = reshape(RichnessTrait ,length(RichnessTrait),1);
clear RemoveLife to take
nrspecs= length(RichnessTrait);

SelectRadi.(genvarname([char(List.Outputs(Trans.prop))])).value(segmentnr,1)  =  sum(RichnessTrait);
SelectRadi.(genvarname([char(List.Outputs(Trans.prop))])).nrspecs(segmentnr,1) =  nrspecs;

Trans.maxi = 37;
Trans.prop = Trans.maxi(1);
Traitmaxs(Trans.maxi) = Trans.Traitmax;
SelectRadi.(genvarname([char(List.Outputs(Trans.prop))])).value(segmentnr,1)  =  log10(sum(RichnessTrait)+1);
SelectRadi.(genvarname([char(List.Outputs(Trans.prop))])).nrspecs(segmentnr,1) =  nrspecs;

clear RichnessTrait CoverTrait StrengthTrait TotalCoverSegment disp i x
clear Traitmax Selects Selects_Cover  nrspecs



%% Amount of Pollinator visits based on Germany. NUMERIC TRAIT, Note summatation over all species!!
Trans.Traitmax = 1;
Trans.maxi = [26,27];
Trans.prop = Trans.maxi(1);
Trans.cover = Trans.maxi(2);
Trans.which = 'Both';
Traitmaxs(Trans.maxi) = Trans.Traitmax;
RichnessTrait = zeros(nrspeciesSegment,Trans.Traitmax);
CoverTrait = zeros(nrspeciesSegment,Trans.Traitmax);
for x = 1:1:nrspeciesSegment
    PresList = find(Input.Pollination.GermanyPlus(:,1) == SpeciesSegment(x,1));
    total_cover = 0;
    if isempty(PresList) ~= 1
         RichnessTrait(x) = sum(Input.Pollination.GermanyPlus(PresList,2:12));
        CoverTrait(x) = SpeciesSegment(x,2);
    else
       RichnessTrait(x) = 0; % Note to include all species, no data is 0 here.
        CoverTrait(x) = 0;
    end
end
clear PresList x
%Remove species without Value
RemoveLife = find(RichnessTrait(:,1) == 9999);
RichnessTrait(RemoveLife,:) = [];
CoverTrait(RemoveLife) = [];
to_take = find(isnan(RichnessTrait) == 0); % So only existing values are taken!!!
RichnessTrait =  RichnessTrait(to_take);
RichnessTrait  = reshape(RichnessTrait ,length(RichnessTrait),1);
CoverTrait =  CoverTrait(to_take);
CoverTrait= reshape(CoverTrait,length(CoverTrait),1);
Nr_covers =  (bsxfun(@times, RichnessTrait,CoverTrait));
clear RemoveLife to take
nrspecs= length(RichnessTrait);

SelectRadi.(genvarname(char(List.Outputs(Trans.prop)))).value(segmentnr,1)  = sum(RichnessTrait);% Note SUM
SelectRadi.(genvarname(char(List.Outputs(Trans.cover)))).value(segmentnr,1)  =  sum(Nr_covers);% Note SUM
SelectRadi.(genvarname(char(List.Outputs(Trans.prop)))).nrspecs(segmentnr,1) =  nrspecs;
SelectRadi.(genvarname(char(List.Outputs(Trans.cover)))).nrspecs(segmentnr,1) =  nrspecs;
clear RichnessTrait CoverTrait StrengthTrait TotalCoverSegment disp i x
clear Traitmax Selects Selects_Cover nrspecs
%% Carbon per leaf weighted. NUMERIC TRAIT, COVER ONLY
Trans.Traitmax = 1;
Trans.maxi =28;
Trans.tf = [1,2,NaN]; % Species list; trait list; strength nr
Trans.prop = Trans.maxi(1);
Trans.cover = Trans.maxi(1);
Traitmaxs(Trans.maxi) = Trans.Traitmax;
Trans.which = 'Cover';
Trans.naming = 'LeafNutrient';
[~,SelectRadi] = StandNum(Input,List,Trans,SpeciesSegment,SelectRadi);

%% DRY Weight per leaf weighted. NUMERIC TRAIT, COVER ONLY
Trans.Traitmax = 1;
Trans.maxi =29;
Trans.tf = [1,6,NaN]; % Species list; trait list; strength nr
Trans.prop = Trans.maxi(1);
Trans.cover = Trans.maxi(1);
Traitmaxs(Trans.maxi) = Trans.Traitmax;
Trans.which = 'Cover';
Trans.naming = 'LeafNutrient';
[~,SelectRadi] = StandNum(Input,List,Trans,SpeciesSegment,SelectRadi);

%% Carbon per ROOT weighted. NUMERIC TRAIT, COVER ONLY
Trans.Traitmax = 1;
Trans.maxi = 30;
Trans.tf = [1,2,NaN]; % Species list; trait list; strength nr
Trans.prop = Trans.maxi(1);
Trans.cover = Trans.maxi(1);
Traitmaxs(Trans.maxi) = Trans.Traitmax;
Trans.which = 'Cover';
Trans.naming = 'RootNutrient';
[~,SelectRadi] = StandNum(Input,List,Trans,SpeciesSegment,SelectRadi);

%% Root DRY MASS Weighted. NUMERIC TRAIT, COVER ONLY
Trans.Traitmax = 1;
Trans.maxi =31;
Trans.tf = [1,6,NaN]; % Species list; trait list; strength nr
Trans.prop = Trans.maxi(1);
Trans.cover = Trans.maxi(1);
Traitmaxs(Trans.maxi) = Trans.Traitmax;
Trans.which = 'Cover';
Trans.naming = 'RootNutrient';
[~,SelectRadi] = StandNum(Input,List,Trans,SpeciesSegment,SelectRadi);

%% Rooting DEPTH weighted. NUMERIC TRAIT, COVER ONLY
Trans.Traitmax = 1;
Trans.maxi =32;
Trans.tf = [1,5,NaN]; % Species list; trait list; strength nr
Trans.prop = Trans.maxi(1);
Trans.cover = Trans.maxi(1);
Traitmaxs(Trans.maxi) = Trans.Traitmax;
Trans.which = 'Cover';
Trans.naming = 'RootNutrient';
[~,SelectRadi] = StandNum(Input,List,Trans,SpeciesSegment,SelectRadi);

%% LandScape diversity
Trans.Traitmax = 1;
Trans.maxi = 33;
Trans.prop = Trans.maxi(1);
Trans.cover = Trans.maxi(1);
Traitmaxs(Trans.maxi) = Trans.Traitmax;
Trans.which = 'Cover';
nrspecs = SelectRadi.Species.value(segmentnr,1);
SelectRadi.(genvarname([char(List.Outputs(Trans.cover))])).value(segmentnr,1)  =  Input.SegmentRecords.Shannon_LU(Input.SegmentRecords.SegmentID == segments_to_do);
SelectRadi.(genvarname([char(List.Outputs(Trans.cover))])).nrspecs(segmentnr,1) =  nrspecs;
clear RichnessTrait CoverTrait StrengthTrait TotalCoverSegment disp i x
clear Traitmax Selects Selects_Cover  nrspecs 
%% Life Form Diversity
Trans.Traitmax = 1;
Trans.maxi = 34;
Trans.prop = Trans.maxi(1);
Trans.cover = Trans.maxi(1);
Traitmaxs(Trans.maxi) = Trans.Traitmax;
RichnessTrait = zeros(nrspeciesSegment,Trans.Traitmax);
CoverTrait = zeros(nrspeciesSegment,Trans.Traitmax);
count = 0;
for x = 1:1:nrspeciesSegment
    PresList = find(Input.LifeFormSpecies(:,1) == SpeciesSegment(x,1));
    if isempty(PresList) ~= 1
        count_start = count +1;
        count = count + length(PresList);
        RichnessTrait(count_start:count) = (Input.LifeFormSpecies(PresList,2));
        CoverTrait(count_start:count) = SpeciesSegment(x,2);
    end
end
clear PresList x count
 %Remove species without No information
RemoveLife = find(RichnessTrait(:,1) == 10);
RichnessTrait(RemoveLife,:) = [];
CoverTrait(RemoveLife) = [];
clear RemoveLife
shannon = 0;
if isempty(RichnessTrait)~= 1
    to_take = find(isnan(RichnessTrait) == 0); % So only existing values are taken!!!
    RichnessTrait =  RichnessTrait(to_take);
    count = 0;
    for life = 1:1:9
       tet = find(RichnessTrait == life);
       if isempty(tet) ~= 1
           count = count + 1;
           Cover_tt(count) = sum(CoverTrait(tet));
       end
    end
    nr_categories = count;
    clear count
    Cover_tt = Cover_tt./sum(Cover_tt); % this is pi
    for i = 1:1:nr_categories
        shannon = shannon + (Cover_tt(i)*log(Cover_tt(i)));
    end
end
nrspecs = length(RichnessTrait);
SelectRadi.(genvarname(char(List.Outputs(Trans.cover)))).value(segmentnr,1)  =  -shannon;
SelectRadi.(genvarname(char(List.Outputs(Trans.cover)))).nrspecs(segmentnr,1) =  nrspecs;
clear RichnessTrait CoverTrait StrengthTrait TotalCoverSegment disp i x
clear Traitmax Selects Selects_Cover nrspecs 

%% Colour Diversity in Number of colours 
Trans.Traitmax = 1;
Trans.maxi = 35;
Trans.prop = Trans.maxi(1);
Trans.cover = Trans.maxi(1);
Traitmaxs(Trans.maxi) = Trans.Traitmax;
RichnessTrait = zeros(nrspeciesSegment,Trans.Traitmax);
CoverTrait = zeros(nrspeciesSegment,Trans.Traitmax);
count = 0;
for x = 1:1:nrspeciesSegment
    PresList = find(Input.ColourSpecies(:,1) == SpeciesSegment(x,1));
    if isempty(PresList) ~= 1
        count_start = count +1;
        count = count + length(PresList);
        RichnessTrait(count_start:count) = (Input.ColourSpecies(PresList,4));
        CoverTrait(count_start:count) = SpeciesSegment(x,2);
    end
end
clear PresList x count   
 %Remove species without No information
RemoveLife = find(RichnessTrait(:,1) == 9999);
RichnessTrait(RemoveLife,:) = [];
CoverTrait(RemoveLife) = [];
clear RemoveLife
shannon = 0;
if isempty(RichnessTrait)~= 1
    to_take = find(isnan(RichnessTrait) == 0); % So only existing values are taken!!!
    RichnessTrait =  RichnessTrait(to_take); %#ok<*FNDSB>
    count = 0;
    for colour = [1:28,9997,9998]
       tet = find(RichnessTrait == colour);
       if isempty(tet) ~= 1
           count = count + 1;
           Cover_tt(count) = sum(CoverTrait(tet)); %#ok<*SAGROW>
       end
    end
   nr_categories = count;
    clear count
    Cover_tt = Cover_tt./sum(Cover_tt); % this is pi
    for i = 1:1:nr_categories
        shannon = shannon + (Cover_tt(i)*log(Cover_tt(i)));
    end
end
nrspecs = length(RichnessTrait);
SelectRadi.(genvarname([char(List.Outputs(Trans.cover))])).value(segmentnr,1)  =  -shannon;
SelectRadi.(genvarname([char(List.Outputs(Trans.cover))])).nrspecs(segmentnr,1) =  nrspecs;
clear RichnessTrait CoverTrait StrengthTrait TotalCoverSegment disp i x
clear Traitmax Selects Selects_Cover nrspecs 
end
