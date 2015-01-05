close all
clear all
clc

folder = uigetdir; % check the help for uigetdir to see how to specify a starting path, which makes your life easier

% get the names of all files. dirListing is a struct array.
dirListing = dir(folder);

% loop through the files and open. Note that dir also lists the directories, so you have to check for them.

x = [];
y = [];
cat = [];

for d = 4:length(dirListing)-1
    if dirListing(1).isdir
        
        fileName = fullfile(folder,dirListing(d).name) % use full path because the folder may not be the active path

        cat = [cat; str2num(fileName(79:81))];
        x = [x; str2num(fileName(83:85))];
        y = [y; str2num(fileName(87:89))];
        
    end
end