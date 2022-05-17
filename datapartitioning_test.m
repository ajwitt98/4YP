%% Data Loading and String Creation
close all
clear all
clc
% Load in Data

[Time, Pressure1_bit, Pressure2_bit] = importfile("C:\Users\Alex Witt\Documents\Engineering work\Fourth Year\4YP\Data_w_Mass\" + ...
    "260_05_2", [1, Inf]);

% Required Constants
stamp_no = length(Time);
min_bit = 102.4;
max_bit = 921.6;
range_bit = max_bit - min_bit;
p_max = 500;

% Smoothing Constants
n_smooth = 1;                  %Number of gradient smoothing loops
interval = 30;                 %Interval over which stability is required
max_drift = 20;                 %Set drift limit
gradient_cap = 0.004;          %Set gradient cap
grad_smoothing = 4000;         %Set gradient vector smoothing value
turning_point = 20;            %Set turning point value
p_upper = 300;                 %Upper relevant pressure
p_lower = 60;                  %Lower relevant pressure
band_smoothing = 500;         %Set the smoothing of the banding function

% Create arrays
T = Time - datetime(Time(1));
outfmt = 'hh:mm:ss.SSSSS';
timetostring = duration(T,'Format',outfmt);
time_raw_1 = seconds(timetostring);
time_raw_2 = seconds(timetostring);
time_raw_PD = seconds(timetostring);
Pressure1_raw = zeros(stamp_no,1);
Pressure2_raw = zeros(stamp_no,1);
PressurePD_raw = zeros(stamp_no,1);
% Adjust Data for range of measurable values
for i = 1:stamp_no
    if Pressure1_bit(i) <= min_bit
        Pressure1_raw(i) = 0;
    elseif Pressure1_bit(i) >= max_bit
        Pressure1_raw(i) = p_max;
    else Pressure1_raw(i) = p_max*((Pressure1_bit(i) - min_bit)/range_bit);
    end

    if Pressure2_bit(i) <= min_bit
        Pressure2_raw(i) = 0;
    elseif Pressure2_bit(i) >= max_bit
        Pressure2_raw(i) = p_max;
    else Pressure2_raw(i) = p_max*((Pressure2_bit(i) - min_bit)/range_bit);
    end

    PressurePD_raw(i) = Pressure1_raw(i) - Pressure2_raw(i);

    if Pressure1_raw(i) == 0
        PressurePD_raw(i) = 0;
    elseif Pressure2_raw(i) == 0
        PressurePD_raw(i) = 0;
    end

end

% Take lengths of vectors
[a] = find(Pressure1_raw);
[b] = find(Pressure2_raw);
[c] = find(PressurePD_raw);
% Create New Arrays
Pressure1 = nonzeros(Pressure1_raw);
Pressure2 = nonzeros(Pressure2_raw);
PressurePD = nonzeros(PressurePD_raw);
% Take time data readings
time1 = time_raw_1(a);
time2 = time_raw_1(b);
timePD = time_raw_1(c);

% Take Rolling Averages
P_1_ravg = movmean(Pressure1,1000);
P_2_ravg = movmean(Pressure2,1000);
P_PD_ravg = movmean(PressurePD,1000);

% Remove irrelevant data points

% This is done on the smoothed data so as not to eliminate
% noise/contributions to legitimate data
P_1_banded = P_1_ravg;
% for i = 1:length(P_1_ravg)
%     if P_1_ravg(i) < p_lower
%         P_1_banded(i) = 0;
%     elseif P_1_ravg(i) > p_upper
%         P_1_banded(i) = 0;
%     end
% end
P_1_banded_ind = find(P_1_banded);
P_1_banded = P_1_banded(P_1_banded_ind);
time1_banded = time1(P_1_banded_ind);
  
% Find the time points and hence indices of where to partition SS data
P_1_mod = P_1_banded;       %Set initial modifiable pressure vector
time1_mod = time1_banded;   %Set initial modifiable time vector
   
for i = 1:n_smooth
    P_grad = zeros(length(P_1_mod),1);
    %Find absolute value of gradient of pressure
    for j = 2:length(P_1_mod)        
        P_grad(j) = Pressure1(j) - Pressure1(j-1);
    end

    P_grad = movmean(P_grad,1000);

    P_turning_point = zeros(length(P_grad),1);
    % Find turning points
    for j = 2:length(P_grad)
        P_turning_point(j) = P_grad(j) - P_grad(j-1);
    end

    P_turning_point = movmedian(P_turning_point,300);
    P_turning_point = movmean(P_turning_point,10000);
    
    %Smooth gradient data
    P_grad = movmean(P_grad,grad_smoothing);

    % Remove values above gradient cap, done outside initial loop so as not
    %to interfere with gradient calculation

    % Also remove data that is fails the turning point criteria

    for k = 1:length(P_grad)
        if abs(P_grad(k)) > gradient_cap
            P_1_mod(k) = 0;
        end
    end   
end

% Plot gradient and turning points
figure(2)
hold on
plot(time1_banded,P_grad)
yyaxis right
plot(time1_banded,P_turning_point)


%%
%Find the indices of the start and end points of steady state data
P_1_mod_ind = find(P_1_mod);
P_1_mod_ind_inv = find(~P_1_mod);
SS_start = zeros(length(P_1_mod),1);
SS_end = zeros(length(P_1_mod),1);
SS_start(1,1) = P_1_mod_ind(1,1);
% SS_end(1,1) = P_1_mod_ind_inv(1,1);
SS_end(end,1) = P_1_mod_ind(end,1);

for j = 2:length(P_1_mod_ind)
    if P_1_mod_ind(j) - P_1_mod_ind(j-1) ~=1
            SS_start(j) = P_1_mod_ind(j) ;
    end
end

for j = 2:length(P_1_mod_ind_inv)
    if P_1_mod_ind_inv(j) - P_1_mod_ind_inv(j-1) ~=1
            SS_end(j) = P_1_mod_ind_inv(j) - 1 ;
    end
end

if SS_start(1,1) == 1
    SS_end(1,1) = P_1_mod_ind_inv(1,1);
end


SS_start = nonzeros(SS_start);
SS_end = nonzeros(SS_end);

%Remove final SS_end value if needed 
 if SS_end(end) == SS_end(end - 1)            
     SS_end(end) = 0;
 end

SS_start = nonzeros(SS_start);
SS_end = nonzeros(SS_end);

%Look at the time duration of steady state data
for i = 1:length(SS_end)
    if time1_mod(SS_end(i)) - time1_mod(SS_start(i)) < interval
        SS_start(i) = 0;
        SS_end(i) = 0;
    end
end
  
%Break up the data into the steady state chunks to record in arrays
SS_start = nonzeros(SS_start);
SS_end = nonzeros(SS_end);
time1_SS_start = time1_mod(SS_start);
time1_SS_end = time1_mod(SS_end);
time1_SS_start_ind = zeros(length(SS_start),1);
time1_SS_end_ind = zeros(length(SS_start),1);
time2_SS_start = time1_SS_start;
time2_SS_end = time1_SS_end;
time2_SS_start_ind = time1_SS_start_ind;
time2_SS_end_ind = time1_SS_end_ind;
timePD_SS_start = time1_SS_start;
timePD_SS_end = time1_SS_end;
timePD_SS_start_ind = time1_SS_start_ind;
timePD_SS_end_ind = time1_SS_end_ind;
P_1_section = cell(length(SS_end),1);
T_1_section = cell(length(SS_end),1);
P_2_section = cell(length(SS_end),1);
T_2_section = cell(length(SS_end),1);
P_PD_section = cell(length(SS_end),1);
T_PD_section = cell(length(SS_end),1);
    
for i = 1:length(SS_start)
    [~,time1_SS_start_ind(i)] = (min(abs(time1 - time1_SS_start(i))));
    [~,time1_SS_end_ind(i)] = (min(abs(time1 - time1_SS_end(i))));
    P_1_section{i} = Pressure1(time1_SS_start_ind(i):time1_SS_end_ind(i));
    T_1_section{i} = time1(time1_SS_start_ind(i):time1_SS_end_ind(i));
    [~,time2_SS_start_ind(i)] = (min(abs(time2 - time2_SS_start(i))));
    [~,time2_SS_end_ind(i)] = (min(abs(time2 - time2_SS_end(i))));
    P_2_section{i} = Pressure2(time2_SS_start_ind(i):time2_SS_end_ind(i));
    T_2_section{i} = time2(time2_SS_start_ind(i):time2_SS_end_ind(i));
    [~,timePD_SS_start_ind(i)] = (min(abs(timePD - timePD_SS_start(i))));
    [~,timePD_SS_end_ind(i)] = (min(abs(timePD - timePD_SS_end(i))));
    P_PD_section{i} = PressurePD(timePD_SS_start_ind(i):timePD_SS_end_ind(i));
    T_PD_section{i} = timePD(timePD_SS_start_ind(i):timePD_SS_end_ind(i));
end

P_1_section_ravg = cell(length(SS_end),1);
for i = 1:length(SS_end)
    P_1_section_ravg{i} = movmean(P_1_section{i},band_smoothing);
end

%Look at the drift of steady state data
for i = 1:length(SS_start)
    if max(P_1_section_ravg{i}) - min(P_1_section_ravg{i}) > max_drift
        P_1_section{i} = [];
        T_1_section{i} = [];
        P_2_section{i} = [];
        T_2_section{i} = [];
        P_PD_section{i} = [];
        T_PD_section{i} = [];
    end
end

Empties = find(cellfun(@isempty,P_1_section));
P_1_section(Empties) = [];   
T_1_section(Empties) = [];
P_2_section(Empties) = [];
T_2_section(Empties) = [];
P_PD_section(Empties) = [];
T_PD_section(Empties) = [];

% Find the average pressure values of smoothed data
Average_Pressures = zeros(length(P_1_section),1);
Average_PD = zeros(length(P_1_section),1);
for i = 1:length(P_1_section)
    Average_Pressures(i,1) = mean(P_1_section{i});
    Average_PD(i,1) = mean(P_PD_section{i});
end


%Plotting Data to visualise smoothing 

figure(1)
hold on
for i = 1:length(P_1_section)
    plot(T_1_section{i},movmean(P_1_section{i},100))
end
plot(time1_mod,P_1_mod)
plot(time2,P_2_ravg)
plot(time1,P_1_ravg)
xlabel('Time (seconds)')
ylabel('Pressure (mbar)')

% yyaxis right
% plot(time1_banded,P_turning_point)
hold off

% figure(2)
% hold on
% plot(time1_banded,P_grad)
% yyaxis right
% plot(time1_banded,P_turning_point)
















