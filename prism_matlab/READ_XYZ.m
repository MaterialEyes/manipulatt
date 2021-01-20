function [ atoms ] = READ_XYZ( filename )
%% Joydeep Munshi - 2020 Dec
% read xyz (pymatgen format) files and convert to atoms array for prismatic

    %[filename, pathname] = uigetfile({'*.xyz'},'File Selector');
    fid = fopen(filename);
    tline = fgetl(fid);
    line_NUM=1;
        while ischar(tline)
            tline = fgetl(fid);
            line_NUM = line_NUM +1;
        end
    TOTAL_LINES = line_NUM-1;
    fclose(fid); 
    
    atoms = zeros(line_NUM-4,4);
    fid = fopen(filename);
    tline = fgetl(fid);
    %display(tline)
    line_NUM=1;
    COUNTER = 1;
    
        while ischar(tline)
            if(line_NUM > 2 && line_NUM < TOTAL_LINES )
                C = strsplit(tline," ");
                atoms(COUNTER,1) = str2double(C{2});
                atoms(COUNTER,2) = str2double(C{3});
                atoms(COUNTER,3) = str2double(C{4});
                atoms(COUNTER,4) = str2double(C{1});
                COUNTER = COUNTER + 1;
            end
            tline = fgetl(fid);
            line_NUM = line_NUM +1;
        end
        
    fclose(fid);  
    %msgbox('File read into COORDINATES');
end
