% Create a function to combine flow rate and pressure data
function combinedData = combine
% close all 
% clear all
% clc
% Import Raw Data
flowData = flowrate_2;
pressureData = pressure;

% Initialise Values
combinedData = cell(length(flowData(:,1)),14);       % Combined Data Cell
initialtime = zeros(length(flowData(:,1)),1);        % Initial time vector
endtime = zeros(length(flowData(:,1)),1);            % End time vector

% Label Each of the Headers 
combinedData{1,1} = 'Run Number';
combinedData{1,2} = 'RPM';
combinedData{1,3} = 'Flow Time';
combinedData{1,4} = 'Flow Rate';
combinedData{1,5} = 'Inlet Pressure Time';
combinedData{1,6} = 'Inlet Pressure';
combinedData{1,7} = 'Outlet Pressure Time';
combinedData{1,8} = 'Outlet Pressure';
combinedData{1,9} = 'Pressure Difference Time';
combinedData{1,10} = 'Pressure Difference';
combinedData{1,11} = 'Inlet Flow Rate';
combinedData{1,12} = 'PPM';
combinedData{1,13} = 'DEF guarantee';
combinedData{1,14} = 'Mean Pressures';


%%
for i = 1:length(flowData(:,1))
    % Fill the immediately easy values
    combinedData{i+1,1} = flowData{i,1};
    combinedData{i+1,2} = flowData{i,2};
    combinedData{i+1,13} = 0;
    combinedData{i+1,14} = 0;

    % Work out the ppm from the name
    str = flowData{i,1};
    expression = '\_';
    splitStr = regexp(str,expression,'split');
    combinedData{i+1,12} = str2double(splitStr{2});

    % Use the rpm value to calculate the inlet flow rate
    RPM = flowData{i,2};
    combinedData{i+1,11} = inlet(RPM);

    % Find the difference between the start and end times 
    flow_start = flowData{i,3}(1);
    pressure_start = pressureData{i,3}(1);
    start_trimtime = abs(flow_start - pressure_start);
    flow_end = flowData{i,3}(end);
    pressure_end = pressureData{i,3}(end-5);
    end_trimtime = (pressure_end - flow_end);
   
    % Create numeric time vectors for comparison
    format long
    ptimenumeric = datenum(pressureData{i,3});
    ftimenumeric = datenum(flowData{i,3});
    
    % Find the start and end points of desired flow data
    f_initialtime = datenum(min(flowData{i,3}) + start_trimtime);
    f_endtime = datenum(max(flowData{i,3}) - end_trimtime);
    
    % Find the index where the smallest difference exists between the time vectors
    [~,f_time_start_index] = (min(abs(ftimenumeric - f_initialtime)));
    [~,f_time_end_index] = (min(abs(ftimenumeric - f_endtime)));

    % Find the times that these correspond to 
    f_start_time = seconds(flowData{i,3}(f_time_start_index) - flowData{i,3}(1));
    f_end_time = seconds(flowData{i,3}(f_time_end_index) - flowData{i,3}(1));

    % Find the corresponding pressure times
    [~,p_time_start_index] = (min(abs(ptimenumeric - f_initialtime)));
    [~,p_time_end_index] = (min(abs(ptimenumeric - f_endtime)));
    p_start_time = seconds(pressureData{i,3}(p_time_start_index) - pressureData{i,3}(1));
    p_end_time = seconds(pressureData{i,3}(p_time_end_index) - pressureData{i,3}(1));

    % Using indices write the data lines into the array
    combinedData{i+1,3} = flowData{i,4}(f_time_start_index:f_time_end_index) - f_start_time;
    combinedData{i+1,4} = flowData{i,5}(f_time_start_index:f_time_end_index);
    
    % Find the indices of each start and end pressure time
    [~,p1_time_start_index] = (min(abs(pressureData{i,7} - p_start_time)));
    p1_start_time = pressureData{i,7}(p1_time_start_index);
    [~,p2_time_start_index] = (min(abs(pressureData{i,8} - p_start_time)));
    p2_start_time = pressureData{i,8}(p2_time_start_index);
    [~,ppd_time_start_index] = (min(abs(pressureData{i,9} - p_start_time)));
    ppd_start_time = pressureData{i,9}(ppd_time_start_index);
    [~,p1_time_end_index] = (min(abs(pressureData{i,7} - p_end_time)));
    [~,p2_time_end_index] = (min(abs(pressureData{i,8} - p_end_time)));
    [~,ppd_time_end_index] = (min(abs(pressureData{i,9} - p_end_time)));

    % Using indices write these data lines in
    combinedData{i+1,5} = pressureData{i,7}(p1_time_start_index:p1_time_end_index) - p1_start_time;
    combinedData{i+1,6} = pressureData{i,4}(p1_time_start_index:p1_time_end_index);
    combinedData{i+1,7} = pressureData{i,8}(p2_time_start_index:p2_time_end_index)- p2_start_time;
    combinedData{i+1,8} = pressureData{i,5}(p2_time_start_index:p2_time_end_index); 
    combinedData{i+1,9} = pressureData{i,9}(ppd_time_start_index:ppd_time_end_index)- ppd_start_time;
    combinedData{i+1,10} = pressureData{i,6}(ppd_time_start_index:ppd_time_end_index); 

end


