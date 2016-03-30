function tempInCelsius = OO_GetBoardTemperature(wrapper)
% tempInCelsius = OO_GetBoardTemperature(wrapper)
%
% Return board temperature in celsius.
%
% 6/28/14   ms    Written.

% Get and return the board temperature.
tmp = wrapper.getFeatureControllerBoardTemperature(0);
tempInCelsius = tmp.getBoardTemperatureCelsius;