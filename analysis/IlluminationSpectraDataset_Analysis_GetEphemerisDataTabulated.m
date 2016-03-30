function [solAz solEl lunAz lunEl lunFrac lunPhase] = IlluminationSpectraDataset_Analysis_GetEphemerisDataTabulated(location, tnum)
% [solAz solEl lunAz lunEl lunFrac lunPhase] = IlluminationSpectraDataset_Analysis_GetEphemerisDataTabulated(location, tnum)
%
% Returns ephemeris data from tables obtained from MICA.
%
% 3/26/2016		spitschan	Added to repository.

basePath = fileparts(mfilename('fullpath'));
ephemerisDataPath = fullfile(basePath, 'data/ephemeris');

A1 = table2cell(readtable(fullfile(ephemerisDataPath, [location '_072014.txt']), 'HeaderLines', 11, 'ReadVariableNames', false));
for i = 1:length(A1)
    tmp = strsplit(A1{i, :});
    timeTmp{i} = [tmp{1} ' ' tmp{2} ' ' tmp{3} ' ' tmp{4}];
    
    [Y,MO,D,H,MI,S] = datevec(datenum(timeTmp{i}, 'yyyy mmm dd HH:MM:SS'));
    H = H;
    timeStamp1{i} = datestr([Y,MO,D,H,MI,S], 30);
    
    phaseAngle(i) = str2double(tmp{5});
    fracIlluminated(i) = str2double(tmp{6});
end


A2 = table2cell(readtable(fullfile(ephemerisDataPath, [location '_MoonLocalZenith.txt']), 'HeaderLines', 12, 'ReadVariableNames', false));
for i = 1:length(A2)
    tmp = strsplit(A2{i, :});
    timeTmp{i} = [tmp{1} ' ' tmp{2} ' ' tmp{3} ' ' tmp{4}];
    
    [Y,MO,D,H,MI,S] = datevec(datenum(timeTmp{i}, 'yyyy mmm dd HH:MM:SS'));
    H = H;
    timeStamp2{i} = datestr([Y,MO,D,H,MI,S], 30);
    
    moonAltitude(i) = 90-str2double(tmp{5});
end

A3 = table2cell(readtable(fullfile(ephemerisDataPath, [location '_SunLocalZenith.txt']), 'HeaderLines', 12, 'ReadVariableNames', false));
for i = 1:length(A3)
    tmp = strsplit(A3{i, :});
    timeTmp{i} = [tmp{1} ' ' tmp{2} ' ' tmp{3} ' ' tmp{4}];
    
    [Y,MO,D,H,MI,S] = datevec(datenum(timeTmp{i}, 'yyyy mmm dd HH:MM:SS'));
    H = H;
    timeStamp3{i} = datestr([Y,MO,D,H,MI,S], 30);
    timeStampNum(i) = datenum(timeTmp{i}, 'yyyy mmm dd HH:MM:SS');
    
    sunAltitude(i) = 90-str2double(tmp{5});
end

timeStampNum = timeStampNum - 4/24;

% Check that all time stamps are the same

if (all(strcmp(timeStamp1, timeStamp2)) && all(strcmp(timeStamp1, timeStamp3)))
    timeStamp = timeStamp1;
else
    error('Time stamps are not the same')
end


% Get sol elevation with MICA
for i = 1:length(tnum)
    solEl(i) = interp1(timeStampNum, sunAltitude, tnum(i));
    lunEl(i) = interp1(timeStampNum, moonAltitude, tnum(i));
    lunFrac(i) = interp1(timeStampNum, fracIlluminated, tnum(i));
    lunPhase(i) = interp1(timeStampNum, phaseAngle, tnum(i));
end

solAz = NaN*ones(1, length(tnum));
lunAz = NaN*ones(1, length(tnum));