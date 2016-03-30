function IlluminationSpectraDataset_Analysis_Preprocess
% IlluminationSpectraDataset_Analysis_Preprocess
%
% Preprocesses the spectra and saves as CSV
%
% 4/1/2015    spitschan       Wrote it.

%% Figure out if we're in the right path
basePath = fileparts(mfilename('fullpath'));
dataRawPath = fullfile(basePath, '../dataraw');
dataPath = fullfile(basePath, '../data');

%% Now check what we have in the measurements folder. The hierarchy is typically:
% data/<location>/<date>
theLocations = dir(dataRawPath);
theLocations = theLocations(arrayfun(@(x) x.name(1), theLocations) ~= '.');
theLocations = theLocations(arrayfun(@(x) x.isdir, theLocations));

%% Find the available dates
for LocIndex = 1:length(theLocations)
    theDates{LocIndex} = dir(fullfile(dataRawPath, theLocations(LocIndex).name));
    theDates{LocIndex} = theDates{LocIndex}(arrayfun(@(x) x.name(1), theDates{LocIndex}) ~= '.');
    theDates{LocIndex} = theDates{LocIndex}([theDates{LocIndex}.isdir]);
end

%% Load the two calibration structures
cal{1} = IlluminationSpectraDataset_SpecCal_GetCalStruct('a');
cal{2} = IlluminationSpectraDataset_SpecCal_GetCalStruct('b');

%% Assign some empty variables
% Set the median filter span
medianFilterSpan = 8;
wlRange = [280 840];

for locIndex = [1 2];
    nTotalSpectra = 0;
    nSaturatedSpectra = 0;
    nNaNZeroSpectra = 0;
    nFinalSpectra = 0;
    
    %% Load the data files one by one
    d0 = 1;
    for d = 1:length(theDates{locIndex})
        fprintf('\n>> Loading from location %s [%g/%g], measurement date %s [%g/%g]', theLocations(locIndex).name, locIndex, length(theLocations), theDates{locIndex}(d).name, d, length(theDates{locIndex}));
        
        %% Reset variables
        c = 0;
        omniRawWls = [];
        omniRawSpectrum = [];
        omniIntegrationTime = [];
        omniBoardTempInCelsius = [];
        omniTimeStamp = [];
        omniWhichMeter = [];
        omniAbsSpectrum = [];
        luminanceVal = [];
        chromaticityXY = [];
        solEl = [];
        lunEl = [];
        lunFrac = [];
        
        %% Iterate over the files
        theFiles = dir(fullfile(dataRawPath, theLocations(locIndex).name, theDates{locIndex}(d).name, '*.mat'));
        fprintf('\n   * Loading file [total %g]: ', length(theFiles));
        for f = 1:length(theFiles)
            % Show counter
            if f > 1
                for j=0:log10(f-1)
                    fprintf('\b');
                end
            end
            fprintf('%d', f);
            
            % Construct the file name
            theFileName = fullfile(theLocations(locIndex).name, theDates{locIndex}(d).name, theFiles(f).name);
            
            % Load in the data
            tmp = load(fullfile(dataRawPath, theFileName));
            
            % Increment the counter
            c = c+1;
            
            % Pull out the interesting variables
            omniIntegrationTime(c) = tmp.meas.integrationTime;
            omniRawWls(:, c) = tmp.meas.wls;
            omniRawSpectrum(:, c) = mean(tmp.meas.spd, 2);
            
            % Correct for times when an intergation time of <1 ms is given.
            % This was done by mistake, and we correct it to be 1 ms, which is
            % the shortest available.
            if omniIntegrationTime(c) < 1000
                omniIntegrationTime(c) = 1000;
            end
            
            % Get some meta data
            omniBoardTempInCelsius(c) = mean(tmp.meas.boardTempInCelsius);
            omniTimeStamp(c) = datenum(tmp.meas.describe.invokeTime, 'yyyymmddTHHMMSS');
            origFileName{c} = theFileName;
            
            % Decide which meter we have
            omniWhichMeter(c) = IlluminationSpectraDataset_SpecCal_GetWhichMeterFromWls(tmp.meas.wls);
        end
        fprintf('. Done.');
        
        %% Get rid of the saturated spectra
        nTotalSpectra = nTotalSpectra + size(omniRawSpectrum, 2);
        fprintf('\n   * Removing saturated spectra... ');
        theSaturation = (omniRawSpectrum == 65535);
        [~, col] = find(theSaturation);
        indicesSaturated = unique(col);
        nSaturatedSpectra = nSaturatedSpectra + length(indicesSaturated);
        
        % Delete the spectra with saturated pixels
        omniRawWls(:, indicesSaturated) = [];
        omniRawSpectrum(:, indicesSaturated) = [];
        omniIntegrationTime(indicesSaturated) = [];
        omniBoardTempInCelsius(indicesSaturated) = [];
        omniTimeStamp(indicesSaturated) = [];
        omniWhichMeter(indicesSaturated) = [];
        fprintf('Done. ');
        if length(indicesSaturated) == 1
            fprintf('%g spectrum removed.', length(indicesSaturated));
        else
        fprintf('%g spectra removed.', length(indicesSaturated));
        end
        % Do dark substraction and absolute calibration
        nSpectra = size(omniRawSpectrum, 2);
        fprintf('\n   * Calibrating spectrum [total %g]: ', nSpectra);
        for i = 1:nSpectra
            % Show counter
            if i > 1
                for j=0:log10(i-1)
                    fprintf('\b');
                end
            end
            fprintf('%d', i);
            [omniAbsSpectrum(:, i), wls] = IlluminationSpectraDataset_SpecCal_RawToAbsolute(cal{omniWhichMeter(i)}, ...
                omniRawSpectrum(:, i), omniRawWls(:, i),...
                omniIntegrationTime(i), omniBoardTempInCelsius(i), ...
                medianFilterSpan);
        end
        fprintf('. Done.');
        
        % Remove stray values
        omniAbsSpectrum(omniAbsSpectrum < 0) = 0;
        fprintf('\n   * Removing NaN or zero spectra... ');
        
        % Delete the all zero or NaN pixes
        indicesNaNZero = find(all((isnan(omniAbsSpectrum)) | (omniAbsSpectrum == 0)));
        omniRawWls(:, indicesNaNZero) = [];
        omniRawSpectrum(:, indicesNaNZero) = [];
        omniIntegrationTime(indicesNaNZero) = [];
        omniBoardTempInCelsius(indicesNaNZero) = [];
        omniTimeStamp(indicesNaNZero) = [];
        omniWhichMeter(indicesNaNZero) = [];
        omniAbsSpectrum(:, indicesNaNZero) = [];
        fprintf('Done. ');
        if length(indicesNaNZero) == 1
            fprintf('%g spectrum removed.', length(indicesNaNZero));
        else
            fprintf('%g spectra removed.', length(indicesNaNZero));
        end
        nNaNZeroSpectra = nNaNZeroSpectra + length(indicesNaNZero);
        nFinalSpectra = nFinalSpectra + size(omniRawSpectrum, 2);
        
        % Get the number of spectra
        nSpectra = size(omniRawSpectrum, 2);
        
        fprintf('\n   * Getting ephemeris data... ');
        % Determine the solar and lunar altitudes, and fraction moon
        % illuminated.
        clear solAz solEl lunAz lunEl lunFrac lunPhase;
        [solAz solEl lunAz lunEl lunFrac lunPhase] = IlluminationSpectraDataset_Analysis_GetEphemerisDataTabulated(theLocations(LocIndex).name, omniTimeStamp);
        fprintf('Done.');
        
        %% Assign everything to a struct
        % Determine whether we have a sunset and and a sunrise. If we have
        % both (i.e. over night), we split it up based on the lowest point
        % of the sun.
        
        % Get the index of the lowest point
        [~, j] = min(solEl);
        if j+1 < nSpectra
            if  (solEl(j+1)) > solEl(j) && (solEl(j-1) > solEl(j)) && d0 ~= 3 && d0 ~= 7
                spectra(locIndex, d0).omniAbsSpectrum = omniAbsSpectrum(:, 1:j);
                spectra(locIndex, d0).omniRawSpectrum = omniRawSpectrum(:, 1:j);
                spectra(locIndex, d0).t = omniTimeStamp(1:j);
                spectra(locIndex, d0).lunFrac = lunFrac(1:j);
                spectra(locIndex, d0).lunPhase = lunPhase(1:j);
                spectra(locIndex, d0).lunEl = lunEl(1:j);
                spectra(locIndex, d0).solEl = solEl(1:j);
                spectra(locIndex, d0).lunAz = lunAz(1:j);
                spectra(locIndex, d0).solAz = solAz(1:j);
                spectra(locIndex, d0).integrationTime = omniIntegrationTime(1:j);
                spectra(locIndex, d0).d = repmat(d0, 1, length(spectra(locIndex, d0).t));
                d0 = d0+1;
                spectra(locIndex, d0).omniAbsSpectrum = omniAbsSpectrum(:, j+1:end);
                spectra(locIndex, d0).omniRawSpectrum = omniRawSpectrum(:, j+1:end);
                spectra(locIndex, d0).t = omniTimeStamp(j+1:end);
                spectra(locIndex, d0).lunFrac = lunFrac(j+1:end);
                spectra(locIndex, d0).lunPhase = lunPhase(j+1:end);
                spectra(locIndex, d0).lunEl = lunEl(j+1:end);
                spectra(locIndex, d0).solEl = solEl(j+1:end);
                spectra(locIndex, d0).lunAz = lunAz(j+1:end);
                spectra(locIndex, d0).solAz = solAz(j+1:end);
                spectra(locIndex, d0).integrationTime = omniIntegrationTime(j+1:end);
                spectra(locIndex, d0).d = repmat(d0, 1, length(spectra(locIndex, d0).t));
            else
                spectra(locIndex, d0).omniAbsSpectrum = omniAbsSpectrum;
                spectra(locIndex, d0).omniRawSpectrum = omniRawSpectrum;
                spectra(locIndex, d0).t = omniTimeStamp;
                spectra(locIndex, d0).lunFrac = lunFrac;
                spectra(locIndex, d0).lunPhase = lunPhase;
                spectra(locIndex, d0).lunEl = lunEl;
                spectra(locIndex, d0).solEl = solEl;
                spectra(locIndex, d0).lunAz = lunAz;
                spectra(locIndex, d0).solAz = solAz;
                spectra(locIndex, d0).integrationTime = omniIntegrationTime;
                spectra(locIndex, d0).d = repmat(d0, 1, length(spectra(locIndex, d0).t));
            end
        else
            spectra(locIndex, d0).omniAbsSpectrum = omniAbsSpectrum;
            spectra(locIndex, d0).omniRawSpectrum = omniRawSpectrum;
            spectra(locIndex, d0).t = omniTimeStamp;
            spectra(locIndex, d0).lunFrac = lunFrac;
            spectra(locIndex, d0).lunPhase = lunPhase;
            spectra(locIndex, d0).lunEl = lunEl;
            spectra(locIndex, d0).solEl = solEl;
            spectra(locIndex, d0).lunAz = lunAz;
            spectra(locIndex, d0).solAz = solAz;
            spectra(locIndex, d0).integrationTime = omniIntegrationTime;
            spectra(locIndex, d0).d = repmat(d0, 1, length(spectra(locIndex, d0).t));
            close all;
        end
        d0 = d0+1;
        
    end
    
    %% Collect all the data
    allSpectra_orig = [spectra(locIndex, :).omniAbsSpectrum];
    allSolEl_orig = [spectra(locIndex, :).solEl];
    allSolAz_orig = [spectra(locIndex, :).solAz];
    allLunEl_orig = [spectra(locIndex, :).lunEl];
    allLunAz_orig = [spectra(locIndex, :).lunAz];
    allLunFrac_orig = [spectra(locIndex, :).lunFrac];
    allLunPhase_orig = [spectra(locIndex, :).lunPhase];
    allTime_orig = [spectra(locIndex, :).t];
    allDates_orig = [spectra(locIndex, :).d];
    
    %% Remove NaN data
    allSolEl_orig(~any(~isnan(allSpectra_orig))) = [];
    allSolAz_orig(~any(~isnan(allSpectra_orig))) = [];
    allLunEl_orig(~any(~isnan(allSpectra_orig))) = [];
    allLunAz_orig(~any(~isnan(allSpectra_orig))) = [];
    allLunFrac_orig(~any(~isnan(allSpectra_orig))) = [];
    allLunPhase_orig(~any(~isnan(allSpectra_orig))) = [];
    allTime_orig(~any(~isnan(allSpectra_orig))) = [];
    allDates_orig(~any(~isnan(allSpectra_orig))) = [];
    allSpectra_orig(:, ~any(~isnan(allSpectra_orig))) = [];
    
    %% Remove all 0 data
    allSolEl_orig(all(omniRawSpectrum == 0)) = [];
    allSolAz_orig(all(omniRawSpectrum == 0)) = [];
    allLunEl_orig(all(omniRawSpectrum == 0)) = [];
    allLunAz_orig(all(omniRawSpectrum == 0)) = [];
    allLunFrac_orig(all(omniRawSpectrum == 0)) = [];
    allLunPhase_orig(all(omniRawSpectrum == 0)) = [];
    allTime_orig(all(omniRawSpectrum == 0)) = [];
    allDates_orig(all(omniRawSpectrum == 0)) = [];
    allSpectra_orig(:, all(omniRawSpectrum == 0)) = [];
    
    %% Save the spectra as CSV
    f = fopen(fullfile(dataPath, [theLocations(locIndex).name '_spectra.csv']), 'w');
    for i = 1:length(allTime_orig)
        fprintf(f, '%s,', datestr(allTime_orig(i)));
    end
    fprintf(f, '\n');
    fclose(f);
    fprintf('\n');
    M = [allDates_orig ; allSolEl_orig ; allLunEl_orig ; allLunFrac_orig ; allSpectra_orig(find(wls == wlRange(1)):find(wls == wlRange(2)), :)];
    dlmwrite(fullfile(dataPath, [theLocations(locIndex).name '_spectra.csv']), M, '-append', 'delimiter', ',', 'precision',10);
    
    f = fopen(fullfile(dataPath, [theLocations(locIndex).name '_dataquality.csv']), 'w');
    fprintf(f, 'N total spectra,N saturated spectra,N NaN or Zero spectra,N final spectra\n');
    fclose(f);
    dlmwrite(fullfile(dataPath, [theLocations(locIndex).name '_dataquality.csv']), [nTotalSpectra nSaturatedSpectra nNaNZeroSpectra nFinalSpectra], 'delimiter', ',');
end