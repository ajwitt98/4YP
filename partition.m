% Function to partition data into steady state chunks
function partitionData = partition(interval,max_drift,gradient_cap,grad_smoothing,band_smoothing)

% close all 
% clear all
% clc

% Smoothing Constants
% interval = 30;                          %Interval over which stability is required
% max_drift = 100;                          %Set drift limit
% gradient_cap = 0.0015;                  %Set gradient cap
% grad_smoothing = 2050;                  %Set gradient vector smoothing value
% band_smoothing = 1000;                   %Set the smoothing of the banding function


% Import combined data 
Data = combine;

% Set initial values
partitionData = cell(length(Data(:,1)),16);
partitionData(1,1:14) = Data(1,1:14);
partitionData{1,15} = 'PD_diff';
partitionData{1,16} = 'PD_diff_ratio';
partitionData{1,17} = 'DEF PD';
partitionData{1,18} = 'DEF FR';


% Loop over every dataset
for x = 2:length(Data(:,1))

    % Set the given values within partitionData
    partitionData{x,1} = Data{x,1};
    partitionData{x,2} = Data{x,2};
    partitionData{x,11} = Data{x,11};
    partitionData{x,12} = Data{x,12};   

    % Definte Vectors for each file
    flow_rate = Data{x,4};
    flowtime = Data{x,3};
    Pressure1 = Data{x,6};
    Pressure2 = Data{x,8};
    PressurePD = Data{x,10};
    time1 = Data{x,5};
    time2 = Data{x,7};
    timePD = Data{x,9};


    % Take Rolling Averages
    P_1_ravg = movmean(Pressure1,1000);
    P_2_ravg = movmean(Pressure2,1000);
    P_PD_ravg = movmean(PressurePD,1000);

    % Set a modifiable banding pressure vector
    P_PD_mod = movmean(PressurePD,500);           %Set initial modifiable pressure vector
    P_1_banding = P_1_ravg;                     %Set initial pressure vector to write zeros to
    time1_mod = time1;                          %Set initial modifiable time vector
    time1_banding = time1;                      %Set initial modifiable time vector
    
       
    %Calculate the gradient and smooth the result
    P_grad = zeros(length(P_PD_mod),1);    
    for j = 2:length(P_PD_mod)        
        P_grad(j) = P_PD_mod(j) - P_PD_mod(j-1);
    end
    P_grad = movmean(P_grad,grad_smoothing);
    
    % % Find turning points
    % P_turning_point = zeros(length(P_grad),1);
    % for j = 2:length(P_grad)
    %     P_turning_point(j) = P_grad(j) - P_grad(j-1);
    % end
    % P_turning_point = movmean(P_turning_point,1000);    
    
    for k = 1:length(P_grad)
        if abs(P_grad(k)) > gradient_cap
            P_1_banding(k) = 0;
        end
    end
    time1_mod = time1_mod(find(P_PD_mod));
    P_PD_mod = nonzeros(P_PD_mod);
    
    %%
    %Find the indices of the start and end points of steady state data
    P_1_banding_ind = find(P_1_banding);
    P_1_banding_ind_inv = find(~P_1_banding);
    SS_start = zeros(length(P_1_banding),1);
    SS_end = zeros(length(P_1_banding),1);
    SS_start(1,1) = P_1_banding_ind(1,1);
    SS_end(end,1) = P_1_banding_ind(end,1);

    for j = 2:length(P_1_banding_ind)
        if P_1_banding_ind(j) - P_1_banding_ind(j-1) ~=1
                SS_start(j) = P_1_banding_ind(j) ;
        end
    end
    
    for j = 2:length(P_1_banding_ind_inv)
        if P_1_banding_ind_inv(j) - P_1_banding_ind_inv(j-1) ~=1
                SS_end(j) = P_1_banding_ind_inv(j) - 1 ;
        end
    end
    
    if SS_start(1,1) == 1
        SS_end(1,1) = P_1_banding_ind_inv(1,1);
    end


    SS_start = nonzeros(SS_start);
    SS_end = nonzeros(SS_end);
    
    % Remove final SS_end value if needed 
     if SS_end(end) == SS_end(end - 1)            
         SS_end(end) = 0;
     end

    SS_start = nonzeros(SS_start);
    SS_end = nonzeros(SS_end);

   
    %Look at the time duration of steady state data
    for i = 1:length(SS_end)
        if time1_banding(SS_end(i)) - time1_banding(SS_start(i)) < interval
            SS_start(i) = 0;
            SS_end(i) = 0;
        end
    end

    SS_start = nonzeros(SS_start);
    SS_end = nonzeros(SS_end);
      
    % Break up the data into the steady state chunks to record in arrays
    SS_start = nonzeros(SS_start);
    SS_end = nonzeros(SS_end);

    % Pressure values share same time values so use this to find the times where SS starts and ends  
    time1_SS_start = time1_banding(SS_start);
    time1_SS_end = time1_banding(SS_end);
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
    flowtime_SS_start = time1_SS_start;
    flowtime_SS_end = time1_SS_end;
    flowtime_SS_start_ind = time1_SS_start_ind;
    flowtime_SS_end_ind = time1_SS_end_ind;
    P_1_section = cell(length(SS_end),1);
    T_1_section = cell(length(SS_end),1);
    P_2_section = cell(length(SS_end),1);
    T_2_section = cell(length(SS_end),1);
    P_PD_section = cell(length(SS_end),1);
    T_PD_section = cell(length(SS_end),1);
    flow_rate_section = cell(length(SS_end),1);
    flowtime_section = cell(length(SS_end),1);

    % Find the corresponding indices where SS starts and ends

    % For pressure data        
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


    %%
    % For flow data
    for i = 1:length(SS_start)
        [~,flowtime_SS_start_ind(i)] = (min(abs(flowtime - flowtime_SS_start(i))));
        [~,flowtime_SS_end_ind(i)] = (min(abs(flowtime - flowtime_SS_end(i))));
        flow_rate_section{i} = flow_rate(flowtime_SS_start_ind(i):flowtime_SS_end_ind(i));
        flowtime_section{i} = flowtime(flowtime_SS_start_ind(i):flowtime_SS_end_ind(i));
    end

        P_PD_section_ravg = cell(length(SS_end),1);
    for i = 1:length(SS_end)
        P_PD_section_ravg{i} = movmean(P_PD_section{i},band_smoothing);
    end
    
    %Look at the drift of steady state data
    for i = 1:length(SS_start)
        if max(P_PD_section_ravg{i}) - min(P_PD_section_ravg{i}) > max_drift
            P_1_section{i} = [];
            T_1_section{i} = [];
            P_2_section{i} = [];
            T_2_section{i} = [];
            P_PD_section{i} = [];
            T_PD_section{i} = [];
        end
    end

    % Remove Empty Cells
    Empties = find(cellfun(@isempty,P_1_section));
    P_1_section(Empties) = [];   
    T_1_section(Empties) = [];
    P_2_section(Empties) = [];
    T_2_section(Empties) = [];
    P_PD_section(Empties) = [];
    T_PD_section(Empties) = [];
    flow_rate_section(Empties) = [];
    flowtime_section(Empties) = [];
    

    %%
    % Find Dead End Filtration Data
    DEF_mean_FR= zeros(length(P_1_section),1);
    DEF_median_FR = zeros(length(P_1_section),1);
    DEF_guarantee = zeros(length(P_1_section),1);
    DEF_PD = zeros(length(P_1_section),1);
    median_PD = zeros(length(P_1_section),1);
    mean_FR = zeros(length(P_1_section),1);
    median_FR = zeros(length(P_1_section),1);
    mode_FR = zeros(length(P_1_section),1);

    for j = 1:length(P_1_section)
        median_PD(j) = median(P_PD_section{j});
        mean_FR(j) = mean(flow_rate_section{j});
        median_FR(j) = median(flow_rate_section{j});
        mode_FR(j) = mode((flow_rate_section{j}));
    end
    partitionData{x,14} = median_PD;

    % Find the values that must be dead end filtration and calculate the
    % average pressure, then subtract this from each SS pressure difference
    for i = 1:length(flow_rate_section)        
        if  median_FR(i) <= 0 
            DEF_guarantee(i) = mean(P_PD_section{i});
        elseif mean_FR(i) <= 0
            DEF_guarantee(i) = mean(P_PD_section{i});
%         elseif mode_FR(i) == 0
%             DEF_guarantee(i) = mean(P_PD_section{i});
        end        
    end
    
    DEF_guarantee_avg = median(nonzeros(DEF_guarantee));
    DEF_guarantee_avg(isnan(DEF_guarantee_avg)) = 0;
    PD_diff = partitionData{x,14} - DEF_guarantee_avg;
    PD_diff_ratio = (100*PD_diff)./partitionData{x,14};
    PD_diff(isnan(PD_diff)) = 0;
    PD_diff_ratio(isnan(PD_diff_ratio)) = 0;
    partitionData{x,13} = DEF_guarantee_avg;
    partitionData{x,15} = PD_diff;
    partitionData{x,16} = PD_diff_ratio;

    % Eliminate the values that are decided to be too small in this case
    for i = 1:length(PD_diff_ratio)
        if median_PD(i) < 5
            P_1_section{i} = [];
            T_1_section{i} = [];
            P_2_section{i} = [];
            T_2_section{i} = [];
            P_PD_section{i} = [];
            T_PD_section{i} = [];
            flow_rate_section{i} = [];
            flowtime_section{i} = []; 
        end
        if PD_diff_ratio(i) < 5  
            DEF_mean_FR(i,1) = mean(flow_rate_section{i});
            DEF_median_FR(i,1) = median(flow_rate_section{i});            
            DEF_PD(i) = mean(P_PD_section{i});
            P_1_section{i} = [];
            T_1_section{i} = [];
            P_2_section{i} = [];
            T_2_section{i} = [];
            P_PD_section{i} = [];
            T_PD_section{i} = [];
            flow_rate_section{i} = [];
            flowtime_section{i} = [];            
        end
        if  median_FR(i) <= 0 
            DEF_mean_FR(i,1) = mean(flow_rate_section{i});
            DEF_median_FR(i,1) = median(flow_rate_section{i});            
            DEF_PD(i) = mean(P_PD_section{i});
            P_1_section{i} = [];
            T_1_section{i} = [];
            P_2_section{i} = [];
            T_2_section{i} = [];
            P_PD_section{i} = [];
            T_PD_section{i} = [];
            flow_rate_section{i} = [];
            flowtime_section{i} = [];            
        elseif mean_FR(i) <= 0
            DEF_mean_FR(i,1) = mean(flow_rate_section{i});
            DEF_median_FR(i,1) = median(flow_rate_section{i});            
            DEF_PD(i) = mean(P_PD_section{i});
            P_1_section{i} = [];
            T_1_section{i} = [];
            P_2_section{i} = [];
            T_2_section{i} = [];
            P_PD_section{i} = [];
            T_PD_section{i} = [];
            flow_rate_section{i} = [];
            flowtime_section{i} = [];            
        end
        partitionData{x,17} = nonzeros(DEF_PD);       
        partitionData{x,18} = nonzeros(DEF_mean_FR);
        
    end



    % Remove Empty Cells
    Empties = find(cellfun(@isempty,P_1_section));
    P_1_section(Empties) = [];   
    T_1_section(Empties) = [];
    P_2_section(Empties) = [];
    T_2_section(Empties) = [];
    P_PD_section(Empties) = [];
    T_PD_section(Empties) = [];
    flow_rate_section(Empties) = [];
    flowtime_section(Empties) = [];
         
    
    % Write partitioned data into cell array
    partitionData{x,3} = flowtime_section;
    partitionData{x,4} = flow_rate_section;
    partitionData{x,6} = P_1_section;
    partitionData{x,8} = P_2_section;
    partitionData{x,10} = P_PD_section;
    partitionData{x,5} = T_1_section;
    partitionData{x,7} = T_2_section;
    partitionData{x,9} = T_PD_section;   
    
     
end




