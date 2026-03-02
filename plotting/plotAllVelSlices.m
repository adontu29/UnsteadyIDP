function plotAllVelSlices()
% plotAllVelSlices  Call plotVelSlices for all .mat files in the folder.

    files = dir('*.mat');   % list all .mat files in current directory

    for k = 1:numel(files)
        fname = files(k).name;
        fprintf('Plotting %s\n', fname);
        plotVelSlices(fname);   % your existing function
    end

end