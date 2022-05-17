% Create a function to write the inlet flow rates into the combine data
function inletFlow = inlet(RPM)
% close all 
% clear all
% clc

% Set constants
m = 0.87302;
c = 0.28107;

% Linear Equation 
inletFlow = m*RPM + c;
