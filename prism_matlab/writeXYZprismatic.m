function [] = writeXYZprismatic(fileName,comment,...
    cellDim,atoms,occArray,uArray)

% Colin Ophus - 2020 Nov

% Write .xyz input file for Prismatic
ID_array = atoms(:,4);
xyz_array = atoms(:,1:3);
if length(occArray) == 1
    occArray = occArray*ones(size(xyz_array,1),1);
end
if length(uArray) == 1
    uArray = uArray*ones(size(xyz_array,1),1);
end

% Initialize file
fid = fopen(fileName,'wt');

% Write comment line (1st)
fprintf(fid,[comment '\n']);

% Write cell dimensions
fprintf(fid,'    %f %f %f \n',cellDim(1:3));

% Write atomic data
dataAll = [ID_array xyz_array occArray uArray];
fprintf(fid,'%d  %f  %f  %f  %d  %f \n',dataAll');

% Write -1 at end of file, for computem compatibility
fprintf(fid,'-1');
fprintf(fid,'\n');

% Close file
fclose(fid);

end