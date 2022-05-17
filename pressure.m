% Create a function that outputs pressure data
function pressureData = pressure
% close all
% clear all
% clc
%Read all the files of the type required
fds = fileDatastore('*.txt', 'ReadFcn', @importpressurefile);
fullFileNames = fds.Files;
num_files = length(fullFileNames);
%Create an array to read files into
pressureData = cell(num_files,9);    
% Data handling section
for x = 1:num_files
    %Import the data from file 
    fileName = fullFileNames{x};    
    [Time, Pressure1_bit, Pressure2_bit] = importpressurefile(fileName, [1, Inf]);
    
    %Define constants
    stamp_no = length(Time);
    min_bit = 102.4;
    max_bit = 921.6;
    range_bit = max_bit - min_bit;
    p_max = 500;

    %Adjust time into usable format
    T = Time - datetime(Time(1));
    outfmt = 'hh:mm:ss.SSSSS';
    timetostring = duration(T,'Format',outfmt);
    time_raw_1 = seconds(timetostring);
        
    %Initialise Arrays
    Pressure1_raw = zeros(stamp_no,1);
    Pressure2_raw = zeros(stamp_no,1);
    PressurePD_raw = zeros(stamp_no,1);
    
    % Adjust Data for range of measurable values and remove zeros
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
    
    singlefullfilename = fullFileNames{x};
    rpmidentifier = regexp(singlefullfilename,'\d\d\d','match','once');
    runidentifier = regexp(singlefullfilename,'\d\d\d\_\w*','match','once');
    rpm = str2double(rpmidentifier);
    pressureData{x,1} = runidentifier;
    pressureData{x,2} = rpm/10;
    pressureData{x,3} = Time;
    pressureData{x,4} = Pressure1;
    pressureData{x,5} = Pressure2;
    pressureData{x,6} = PressurePD;
    pressureData{x,7} = rmmissing(time1);
    pressureData{x,8} = rmmissing(time2);
    pressureData{x,9} = rmmissing(timePD);   

end