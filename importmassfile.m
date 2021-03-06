function [Time1, Data] = importmassfile(filename, dataLines)
%IMPORTFILE Import data from a text file
%  [TIME1, DATA] = IMPORTMASSFILE(FILENAME) reads data from text file
%  FILENAME for the default selection.  Returns the data as column
%  vectors.
%
%  [TIME1, DATA] = IMPORTMASSFILE(FILE, DATALINES) reads data for the
%  specified row interval(s) of text file FILENAME. Specify DATALINES as
%  a positive scalar integer or a N-by-2 array of positive scalar
%  integers for dis-contiguous row intervals.
%
%  Example:
%  [Time1, Data] = importmassfile("C:\Users\Alex Witt\Documents\Engineering work\Fourth Year\4YP\Data_w_Mass\Mass Data\massbalancetest.csv", [2, Inf]);
%
%  See also READTABLE.
%
% Auto-generated by MATLAB on 12-Apr-2022 13:31:19

%% Input handling

% If dataLines is not specified, define defaults
if nargin < 2
    dataLines = [2, Inf];
end

%% Set up the Import Options and import the data
opts = delimitedTextImportOptions("NumVariables", 11);

% Specify range and delimiter
opts.DataLines = dataLines;
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["Var1", "Var2", "Time1", "Var4", "Data", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11"];
opts.SelectedVariableNames = ["Time1", "Data"];
opts.VariableTypes = ["string", "string", "datetime", "string", "double", "string", "string", "string", "string", "string", "string"];

% Specify file level properties
opts.ImportErrorRule = "omitrow";
opts.MissingRule = "omitrow";
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, ["Var1", "Var2", "Var4", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11"], "WhitespaceRule", "preserve");
opts = setvaropts(opts, ["Var1", "Var2", "Var4", "Var6", "Var7", "Var8", "Var9", "Var10", "Var11"], "EmptyFieldRule", "auto");
opts = setvaropts(opts, "Time1", "InputFormat", "HH:mm:ss");

% Import the data
tbl = readtable(filename, opts);

%% Convert to output type
Time1 = tbl.Time1;
Data = tbl.Data;
end