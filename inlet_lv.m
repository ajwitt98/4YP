% Create a function to write the inlet flow rates into the combine data
function inletFlow = inlet_lv(RPM)
% close all 
% clear all
% clc

% Set testing RPM value
% RPM = 48;

% Set constants
m = 0.4880;
c = 0.-0.5834;

% Linear Equation 
inletFlow = m*RPM + c;
