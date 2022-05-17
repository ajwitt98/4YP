% Create a function that outputs flow rate data
function flowData = flowrate_2
% close all
% clear all
% clc
% Import raw mass data
% Read all the files of the type required
fds = fileDatastore('*.csv', 'ReadFcn', @importmassfile);
fullFileNames = fds.Files;
num_files = length(fullFileNames);

% Create an array to read files into
flowData = cell(num_files,6);

%%

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
    unique_time = unique(time_raw_1);
    flow_vect = cell(length(unique_time),3);
    raw_flow = zeros(length(Mass)-1,1);
    flow_time = time_raw_1(2:end);

    for i = 2:length(Mass)
        raw_flow(i-1) = Mass(i) - Mass(i-1);
    end

    for i = 1:length(unique_time)
        for j = 1:length(flow_time)
            if flow_time(j) == unique_time(i)
                flow_vect{i,1}(j,1) = raw_flow(j);
            end
        end
        flow_vect{i,2} = sum(cell2mat(flow_vect(i,1)));
        flow_vect{i,3} = flow_vect{i,2}*60;
    end

    flow_rate = cell2mat(flow_vect(:,3));
    flow_time = unique_time;


    singlefullfilename = fullFileNames{x};
    rpmidentifier = regexp(singlefullfilename,'\d\d\d','match','once');
    runidentifier = regexp(singlefullfilename,'\d\d\d\_\w*','match','once');
    rpm = str2double(rpmidentifier);
    flowData{x,1} = runidentifier;
    flowData{x,2} = rpm/10;
    flowData{x,3} = unique(Time(2:end));
    flowData{x,4} = flow_time;
    flowData{x,5} = flow_rate;
    flowData{x,6} = mean(flow_rate);
end