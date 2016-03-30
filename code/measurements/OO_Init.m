function [spectrometerObj, wrapper] = OO_Init()
% spectrometerObj = OO_Init()
%
% Initializes the USB2000+. Returns the spectrometer object.
%
% 6/26/14     ms      Written.

% Load the mdd driver
fprintf('> Loading OceanOptics_OmniDriver.mdd ... ');
spectrometerObj = icdevice('OceanOptics_OmniDriver.mdd');
fprintf('done.\n');

% Connect to the object
fprintf('> Connecting to object ... ');
tic;
connect(spectrometerObj);
fprintf('done: ');
toc

% Also open it up with the wrapper. This is a bit of a kluge, but we do
% that to get the board temperature measurements, which do not appear to be
% available with the IC toolbox method.
wrapper = com.oceanoptics.omnidriver.api.wrapper.Wrapper();
wrapper.openAllSpectrometers();


end