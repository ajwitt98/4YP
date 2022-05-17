% Write a script to statistically handle and plot data
close all
clear all
clc
% Set partitioning parameters
n_smooth = 2;              % Number of gradient smoothing loops
interval = 30;              % Interval over which stability is required.
max_drift = 8;             % Set drift limit
gradient_cap = 0.020;       % Set gradient cap
grad_smoothing = 750;        % Set gradient vector smoothing value
band_smoothing = 1000;      % Set the smoothing of the banding function
p_upper = 300;              % Upper relevant pressure
p_lower = 60;               % Lower relevant pressure

% Import data sets
rawData = combine;
partitionData = partition(n_smooth,interval,max_drift,gradient_cap,grad_smoothing,band_smoothing,p_upper,p_lower);

%%
% Useful constants 
n_runs = length(partitionData(:,1))-1;
rpm_col = rawData(:,2);
rpm_col = cell2mat(rpm_col(2:end));
n_rpms = length(unique(rpm_col));

%%
% Evaluate effectiveness of partitioning constants
n_SS_sections = zeros(length(partitionData(:,1))-1,1);
avg_SS_run_time = zeros(length(partitionData(:,1))-1,1);
avg_inlet_pressure = zeros(length(partitionData(:,1))-1,1);
avg_pressure_difference = zeros(length(partitionData(:,1))-1,1);
SS_run_time = cell(length(partitionData(:,1))-1,2);
SS_inlet_pressure = cell(length(partitionData(:,1))-1,2);
SS_pressure_difference = cell(length(partitionData(:,1))-1,2);
SS_flow_rate = cell(length(partitionData(:,1))-1,2);

for i = 2:length(partitionData(:,1))

    % Set the rpm labels and find number of SS sections
    n_SS_sections(i-1) = length(partitionData{i,3});
    SS_run_time{i-1,1} = partitionData{i,2};
    SS_inlet_pressure{i-1,1} = partitionData{i,2};
    SS_pressure_difference{i-1,1} = partitionData{i,2};
    SS_flow_rate{i-1,1} = partitionData{i,2};
    
    % Calculate the values for each steady state section
    for j = 1:length(partitionData{i,3})
        SS_run_time{i-1,2}{j,1} = max(partitionData{i,3}{j} - min(partitionData{i,3}{j}));
        SS_inlet_pressure{i-1,2}{j,1} = mean(partitionData{i,6}{j});
        SS_pressure_difference{i-1,2}{j,1} = mean(partitionData{i,10}{j});
        SS_flow_rate{i-1,2}{j,1} = mean(partitionData{i,4}{j});
    end 
    
    % Calculate the average SS values for each bioreactor/rpm
    avg_SS_run_time(i-1) = mean(cell2mat(SS_run_time{i-1,2}));
    avg_inlet_pressure(i-1) = mean(cell2mat(SS_inlet_pressure{i-1,2}));
    avg_pressure_difference(i-1) = mean(cell2mat(SS_pressure_difference{i-1,2}));
    avg_flow_rate(i-1) = mean(cell2mat(SS_flow_rate{i-1,2}));

end

%%

% Create a vector of zeros to write the pressure differences into 
p = zeros(length(SS_pressure_difference(:,1)),1);
for i = 1:length(SS_pressure_difference(:,1))
    p(i) = length(SS_pressure_difference{i,2});
end
p_diff_string = zeros(sum(p),1);
rpm_value = zeros(sum(p));
p_inlet_string = zeros(sum(p),1);
for j = 1:length(SS_pressure_difference(:,1))
    for k = 1:length(SS_pressure_difference{j,2})        
        p_diff_string(j,k) = SS_pressure_difference{j,2}{k};
        p_inlet_string(j,k) = SS_inlet_pressure{j,2}{k};
        rpm_value(j,k) = SS_pressure_difference{j,1};
    end
end
 
rpm_value = nonzeros(rpm_value);
p_diff_string = nonzeros(p_diff_string);
p_inlet_string = nonzeros(p_inlet_string);

%%
% Plot Results
figure(1)
plot(rpm_value,p_diff_string,'x')

%%
figure(2)
v = violinplot(p_diff_string,rpm_value, ...
    'ShowData',false,'Width',0.18,'ViolinAlpha',0.6);
xlabel('Pump RPM')
ylabel('Pressure drop along lumen (mbar)')

%%
% Plot Flow Rate Curves 
clist = colormap(hsv(n_runs));
rsq = zeros(n_runs,1);
lobf = zeros(n_runs,2);


figure(3)
for i = 1:n_runs
    plot(cell2mat(SS_inlet_pressure{i,2}),cell2mat(SS_flow_rate{i,2}),'x')
    c = polyfit(cell2mat(SS_inlet_pressure{i,2}),cell2mat(SS_flow_rate{i,2}),1);
    lobf(i,1) = c(1);
    lobf(i,2) = c(2);
    % Display evaluated equation y = m*x + b
    disp(['Equation is y = ' num2str(c(1)) '*x + ' num2str(c(2))])
    % Evaluate fit equation using polyval
    y_est = polyval(c,cell2mat(SS_inlet_pressure{i,2}));
    % Add trend line to plot
    hold on
    plot(cell2mat(SS_inlet_pressure{i,2}),y_est,'color',clist(i,:),'LineWidth',2)
    hold off
    str = num2str(SS_inlet_pressure{i});
    legend([str 'rpm'])
    
    % Evaluate R^2 Value
    yresid = cell2mat(SS_flow_rate{i,2}) - y_est;
    SSresid = sum(yresid.^2);
    SStotal = (length(cell2mat(SS_flow_rate{i,2}))-1) * var(cell2mat(SS_flow_rate{i,2}));
    rsq(i) = 1 - SSresid/SStotal;

end
xlabel('Inlet Pressure (mbar)')
ylabel('Flow Rate (ml/min)')

%%
figure(4)
plot(p_inlet_string,p_diff_string,'x')


