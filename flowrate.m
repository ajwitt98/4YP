% Create a function that outputs flow rate data
function flowData = flowrate
% close all
% clear all
% clc
% Import raw mass data
% Read all the files of the type required
fds = fileDatastore('*.csv', 'ReadFcn', @importmassfile);
fullFileNames = fds.Files;
num_files = length(fullFileNames);

% Create an array to read files into
flowData = cell(num_files,5);

% Loop over every file
for x = 1:num_files
    % Import each of the files
    fileName = fullFileNames{x};    
    [Time,Mass] = importmassfile(fileName, [1, Inf]);

    % Convert time into correct format
    T = Time - datetime(Time(1));
    outfmt = 'hh:mm:ss.SSSSS';
    timetostring = duration(T,'Format',outfmt);
    time_raw_1 = seconds(timetostring);

    % Break time from seconds into fractions of a second
    time_change_loc = zeros(length(time_raw_1),1);
    time_change_loc(end) = length(time_raw_1);
    
    for i = 2:length(time_raw_1)
        if time_raw_1(i) - time_raw_1(i-1) == 1
            time_change_loc(i) = i - 1;            
        end
    end

    time_change_loc = nonzeros(time_change_loc);
    time_change_size = zeros(length(time_change_loc),1);
    time_change_size(1) = time_change_loc(1);
    
    for i = 2:length(time_change_loc)
        time_change_size(i) = (time_change_loc(i) - time_change_loc(i-1));
    end
    
    time_just_seconds = [0:1:max(time_raw_1)].';
    increment = 1./time_change_size;
    
    time_cell = cell(length(time_just_seconds),1);
    
    for i = 1:length(time_change_loc)-1
        time_cell{i} = [time_just_seconds(i):increment(i):time_just_seconds(i+1)-increment(i)].';
    end
    time_cell{end} = max(time_just_seconds);
    flowtime = cell2mat(time_cell);
    

    % Actually calculate flow rate with given timesteps 
    flow_rate = zeros(length(flowtime),1);
    flow_time_increment = zeros(length(flowtime),1);
    for i = 2:length(flowtime)
        flow_time_increment(i) = flowtime(i) - flowtime(i - 1);
        flow_rate(i) = 60*(Mass(i) - Mass(i - 1))/flow_time_increment(i);
    end

    singlefullfilename = fullFileNames{x};
    rpmidentifier = regexp(singlefullfilename,'\d\d\d','match','once');
    runidentifier = regexp(singlefullfilename,'\d\d\d\_\w*','match','once');
    rpm = str2double(rpmidentifier);
    flowData{x,1} = runidentifier;
    flowData{x,2} = rpm/10;
    flowData{x,3} = Time;
    flowData{x,4} = flowtime;
    flowData{x,5} = flow_rate;

end


