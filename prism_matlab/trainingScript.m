% Joydeep Munshi - 2021 Jan
% script for generating training data for the given input atoms supercells

%f = waitbar(0, 'Starting');
clear
clc
cellDim = [50 50 500];
myDir = './supercell';
myFiles = dir(fullfile(myDir,'*.mat'));

for k = 1:length(myFiles)
    %waitbar(k/length(myFiles), f, sprintf('Calculating: %d %%', floor(k/(length(myFiles))*100)));
    baseFileName = myFiles(k).name;
    fullFileName = fullfile('supercell', baseFileName);
    load('Atoms.mat')
    fprintf('Reading: %s \n', baseFileName);
    supercell = load(fullFileName);
    
    tic
    if supercell.label == "unrotated"
        coord = supercell.supercell_unrotated;
    elseif supercell.label == "zoneaxis"
        coord = supercell.supercell_zoneaxis;
    else
        coord = supercell.supercell_rotated;
    end
    
    emdSTEM = PRISM01_potential(coord,cellDim);
    counter = 0;   

    for i = [0:0.2:1.8] / 1000
         [emdSTEM,EWamp] = PRISM02_multislice(emdSTEM, i);
         sTraining = processPotentialML(emdSTEM);
         
         outFile = strsplit(baseFileName,'.');
         counter = counter + 1;
         tiltFile = strcat('_tilt',int2str(counter));
         outFileName = strcat(outFile(1),strcat(tiltFile, '.mat'));
        
         outFileName = string(fullfile('training', outFileName(1)));
         disp(outFileName)
         probe = sTraining.dataProbe;
         meas = sTraining.dataMeas;
         pot = sTraining.dataPot;
         qz = sTraining.dataQz;    
         save(outFileName,'probe','meas','pot','qz');
    end

    clearvars -except myFiles cellDim
    close all
    toc
end

close(f)

