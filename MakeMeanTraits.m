function MakeMeanTraits(Input,List)
% select relevant Species only and Order to traits;
% Segments are compared to their respective Focal sites as proportions.
    for landscape = 1:1:36 % per landscape to make 1:1 comparisons
        clc
        display('Running radiation statistics')
        display(landscape)
        Segments = Input.SegmentRecords((Input.SegmentRecords.Landscape(1:351) == landscape),1);
        Target = Input.SegmentRecords((Input.SegmentRecords.SegmentID(1:387) == (1000 + landscape)),1); %#ok<*NBRAK>
        Segments = [Target;Segments];
        segmentnr = 0;
        SelectRadi.whatever = [];
        %cut = 0
        for count = 1:1:length(Segments) % So per segments
            segments_to_do = Segments.SegmentID(count);
            segmentnr = segmentnr + 1;
            SpeciesSegmentFind = find(Input.GIReleveesSpecies(:,2) == segments_to_do);
            clear  SpeciesSegment
            clear SpeciesSegmentID
            clear SpeciesSegmentID2
            SpeciesSegment = Input.GIReleveesSpecies(SpeciesSegmentFind,[1,4,7]); % 1) Species,2) cover,3)linear/areal
            nrspeciesSegment = length(SpeciesSegment(:,1));
            Descriptives = Input.SegmentRecords((Input.SegmentRecords.SegmentID == segments_to_do),[1,3,2,14,10,11,15,16,17,20,21,28]);
            % 1.Segment ID, 2. Band,3. LandscapeID, 4. Country,
            % 5. linearGI, 6. ArealGI, 7 Road GI;
            % 8 Other linear GI 9 Segment yes/no 10 Shannon index 11 Total GI;
            % 12. Size
            Distances = Input.SegmentRecords((Input.SegmentRecords.SegmentID == segments_to_do),[4,5,8,9,18,19,22,23,24,25,26,27,42]);
            % 1 Mean Current Cattle; 2 Minumum current Cattle; 3 Mean wind;
            % 4 Min Wind 5 Min Road 6 Mean road 7 Min
            % Birds; 8 Mean birds; 9 min Human 10 mean
            % Human 11 min Mammal; 12 mean mammal;
            % 13 NO distance
            % make it one value per segment;
            SpeciesSegmentIDs(segmentnr,:) =[Descriptives, Distances]; %#ok<*SAGROW,*AGROW>
            
            %% Select only Focal site species for the respective Focal site
            % Check for species occurence in our list of species (so not identified
            % -9xx numbers are dropping out here).
            icount = 0;
            for x = 1:1:nrspeciesSegment
                PresList = find(Input.SpeciesList.AccSpecies(:) == SpeciesSegment(x,1));
                if isempty(PresList) ~= 1
                    icount =  icount +1;
                    PresenceSegment(icount,1) = SpeciesSegment(x,1); % AccSpecies nr
                    PresenceSegment(icount,2) = SpeciesSegment(x,2); % Abundance
                    PresenceSegment(icount,3) = SpeciesSegment(x,3); % Linear (1) or Area (2)
                end
            end
            clear PresList
            clear SpeciesSegment %!!!! Note name used again below, but that one cleaned
            
            % Select the species to work further with, depending on the
            % choice of type of species
            nrspeciesSegment = length( PresenceSegment(:,1));
            AllTargetSpecies = Input.SpeciesList.AccSpecies(Input.SpeciesList.TargetPresence ==1); % OPTION 2
            Species_Not_all_Target = Input.SpeciesList.AccSpecies(Input.SpeciesList.TargetPresence ==0); % OPTION 5
            AllSpecies = Input.SpeciesList.AccSpecies; % OPTION 3
            if  segments_to_do > 1000
                % so this is the actual determination of which are the target species, given the target is the first site in the series
                Species_this_Target =  PresenceSegment(:,1); % OPTION 1
                Species_this_Target(Species_this_Target == 0) = [];
                Species_Not_this_Target = AllSpecies;
                SpeciesSegment(:,1) = PresenceSegment(:,1); % AccSpecies nr
                SpeciesSegment(:,2) = PresenceSegment(:,2); % Abundance
                SpeciesSegment(:,3) = PresenceSegment(:,3);  % Linear (1) or Area (2)
                for x = nrspeciesSegment:-1:1
                    lst = find(Species_Not_this_Target == Species_this_Target(x));
                    Species_Not_this_Target(lst) = [];
                    clear lst
                end
            else
                ecount = 0;
                for x = 1:1:nrspeciesSegment
                    if Input.species_option == 1 % OPTION 1:
                        % determine whether the species occurs in the SINGLE TARGET as
                        % well, else don't include the species
                        PresList = find(Species_this_Target == PresenceSegment(x,1));
                    elseif Input.species_option == 2  % OPTION 2
                        % determine whether the species occurs in ANY TARGET
                        PresList = find(AllTargetSpecies == PresenceSegment(x,1));
                    elseif Input.species_option == 3  % OPTION 3:
                        %  Use ALL species
                        PresList = find(AllSpecies == PresenceSegment(x,1));
                    elseif  Input.species_option == 4  % OPTION 4:
                        % Species not in This argets
                        PresList = find(Species_Not_this_Target == PresenceSegment(x,1));
                    elseif  Input.species_option == 5  % OPTION 5:
                        % Incoming species: npt present in Targets
                        PresList = find(Species_Not_all_Target == PresenceSegment(x,1));
                    end
                    if isempty(PresList)~= 1
                        ecount = ecount + 1;
                        SpeciesSegment(ecount,1) = PresenceSegment(x,1);  %AccSpecies nr
                        SpeciesSegment(ecount,2) = PresenceSegment(x,2);  % Abundance
                        SpeciesSegment(ecount,3) = PresenceSegment(x,3); % Linear (1) or Area (2)
                    end  
                end
                if Input.Linear_option == 1 % Linear elements only
                    SpeciesSegment((SpeciesSegment(:,3) == 2),:) = [];
                elseif Input.Linear_option == 2 % Area elements only
                    SpeciesSegment((SpeciesSegment(:,3) == 1),:) = [];
                else
                    %all species are selected, linear and area
                end
            end
            clear PresenceSegment
            if exist('SpeciesSegment') == 1
                nrspeciesSegment = length(SpeciesSegment(:,1));
            else
                nrspeciesSegment = 0;
            end
            %%
            % Assign all traits based on the species within the segment; same
            % procedure as for All species differences. So this is per segment
            % and afterwards all analyses are based on proportional or numerical means
            % per traittype per segment
            [SelectRadi,maxtraits] =  CollateTraits(Input,List,SpeciesSegment,nrspeciesSegment,segmentnr,segments_to_do,SelectRadi);
            clear SpeciesSegment
        end
        clear Species_Not_this_Target
        % Make proportions and put in readable formats per traittype within
        % traits
        for trait = 1:List.traits_to_do;
            maxtrait = maxtraits(trait);
            for traittype = 1:maxtrait
                Traitmatrix = SelectRadi.(genvarname([char(List.Outputs(trait))])).value(:,traittype);
                Specmatrix = SelectRadi.(genvarname([char(List.Outputs(trait))])).nrspecs(:,traittype);
                Sizematrix =  SpeciesSegmentIDs.Size;
                if Traitmatrix(1) ~= 0
                    Prop_presence(:,traittype) = Traitmatrix./Traitmatrix(1); % THIS IS THE KEY LINE: normalising and making the proportions
                    Tot_presence(:,traittype) = SelectRadi.(genvarname([char(List.Outputs(trait))])).value(:,traittype);
                    nr_species(:,traittype) = SelectRadi.(genvarname([char(List.Outputs(trait))])).nrspecs(:,traittype);
                    if trait == 1
                        Spec_loss = Specmatrix./Specmatrix(1);
                        Size_prop = Sizematrix./Sizematrix(1);
                    end
                else
                    Prop_presence(:,traittype) = NaN;
                    Tot_presence(:,traittype) = NaN;
                    nr_species(:,traittype) = NaN;
                end
                clear traitmatrix
            end
            % Make it nice joined cell arrays
            if landscape == 1
                PropRadi.(genvarname([char(List.Outputs(trait))])) = Prop_presence;
                TotRadi.(genvarname([char(List.Outputs(trait))])) =  Tot_presence;
                NrSpec.(genvarname([char(List.Outputs(trait))])) =  nr_species;
                if trait  == 1
                    SpecLoss =  Spec_loss;
                    SizeProp = Size_prop;
                end
            else
                PropRadi.(genvarname([char(List.Outputs(trait))])) =  [PropRadi.(genvarname([char(List.Outputs(trait))]));Prop_presence];
                TotRadi.(genvarname([char(List.Outputs(trait))])) =  [TotRadi.(genvarname([char(List.Outputs(trait))]));Tot_presence];
                NrSpec.(genvarname([char(List.Outputs(trait))])) = [NrSpec.(genvarname([char(List.Outputs(trait))]));nr_species];
                if trait  == 1
                    SpecLoss =  [SpecLoss;Spec_loss];
                    SizeProp = [SizeProp;Size_prop];
                end
            end
            clear Prop_presence Tot_presence nr_species Spec_loss
        end
        if landscape == 1
            PropRadi.nrTargetspeciessegment = SelectRadi.nrTargetspeciessegment;
            PropRadi.SpeciesSegmentIDs = SpeciesSegmentIDs;
            
        else
            PropRadi.nrTargetspeciessegment =  [PropRadi.nrTargetspeciessegment;SelectRadi.nrTargetspeciessegment];
            PropRadi.SpeciesSegmentIDs = [ PropRadi.SpeciesSegmentIDs;SpeciesSegmentIDs];
        end
        clear SelectRadi Selection SpeciesSegmentIDs
    end % per landscape 
    %%
    save('Radiation_proportions', 'List','PropRadi','TotRadi','NrSpec','maxtraits', 'Input','SpecLoss','SizeProp','Species_Not_all_Target','AllTargetSpecies')
end