clear; clc; close all;

for camid = 1:8
    % 1-8 patients
    num = sprintf('%02d', camid);

    % Read the file
    filename = ['CLData/AP04-CAM-' num '-XT-CLOSED.txt'];  
    disp(filename)
    fid = fopen(filename, 'r');
    fileContent = textscan(fid, '%s', 'Delimiter', '\n');
    fclose(fid);
    lines = fileContent{1};
    
    % Struture setup
    data = struct();
    data.ID = '';
    data.Weight_KG = [];
    data.Basal_UDay = [];
    data.Basal_rate_per30min_Uh = zeros(48, 1);
    
    % Basic information
    partid = strsplit(lines{1}, ':');  % ID
    data.ID = strtrim(partid{2});
    partweight = strsplit(lines{2}, ':');  % Weight (KG)
    data.Weight_KG = sscanf(partweight{2}, '%f');
    partbasal = strsplit(lines{3}, ':');   % Basal (U/day)
    data.Basal_UDay = sscanf(partbasal{2}, '%f');
    
    % Basal rate in 30min steps (U/h)
    idxbr = find(contains(lines, 'Basal rate in 30min steps'));
    partbr = lines(idxbr+1:idxbr+4);
    partbrnum = strjoin(partbr, ' ');
    partbrnum = sscanf(partbrnum, '%f');
    data.Basal_rate_per30min_Uh = reshape(partbrnum, [], 1);
    
    % Enteral_bolus (meal)
    idxeb = find(contains(lines, 'Enteral_bolus'));
    start_idxeb = idxeb + 3;
    cell_eb = [];
    for i = start_idxeb:length(lines)
        if strlength(strtrim(lines(i))) == 0
            break;  % space
        end
        cell_eb = [cell_eb; lines(i)];
    end
    neb = length(cell_eb);
    data.Enteral_bolus = cell(neb, 3);
    for i = 1:neb   % Date   Time   CHO(g)
        tokens = regexp(cell_eb(i), '(\d{2}/\d{2}/\d{4})\s+(\d{2}:\d{2})\s+([\d.]+)', 'tokens');
        if ~isempty(tokens)
            data.Enteral_bolus{i,1} = tokens{1}{1}{1};             
            data.Enteral_bolus{i,2} = tokens{1}{1}{2};              
            data.Enteral_bolus{i,3} = str2double(tokens{1}{1}{3});  
        end
    end
    
    % Insulin_bolus
    idxib = find(contains(lines, 'Insulin_bolus'));
    start_idxib = idxib + 3;
    cell_ib = [];
    for i = start_idxib:length(lines)
        if strlength(strtrim(lines(i))) == 0
            break;  % space
        end
        cell_ib = [cell_ib; lines(i)];
    end
    nib = length(cell_ib);
    data.Insulin_bolus = cell(nib, 5);
    for i = 1:nib   % Date   Time   Bolus(U)   Dutration(min)   Insulin(S|R|N)
        tokensib = regexp(cell_ib(i), '(\d{2}/\d{2}/\d{4})\s+(\d{2}:\d{2})\s+([\d.]+)\s+(\d+)\s+([A-Za-z])', 'tokens');
        if ~isempty(tokensib)
            data.Insulin_bolus{i,1} = tokensib{1}{1}{1};             
            data.Insulin_bolus{i,2} = tokensib{1}{1}{2};              
            data.Insulin_bolus{i,3} = str2double(tokensib{1}{1}{3});  
            data.Insulin_bolus{i,4} = str2double(tokensib{1}{1}{4}); 
            data.Insulin_bolus{i,5} = tokensib{1}{1}{5};
        end
    end
    
    % Insulin_infusion
    idxii = find(contains(lines, 'Insulin_infusion'));
    start_idxii = idxii(1) + 3;
    cell_ii = [];
    for i = start_idxii:length(lines)
        if strlength(strtrim(lines(i))) == 0
            break;  % space
        end
        cell_ii = [cell_ii; lines(i)];
    end
    nii = length(cell_ii);
    data.Insulin_infusion = cell(nii, 4);
    for i = 1:nii   % Date   Time   Rate(U/h)   Insulin(S|R)
        tokensii = regexp(cell_ii(i), '(\d{2}/\d{2}/\d{4})\s+(\d{2}:\d{2})\s+([\d.]+)\s+([A-Za-z])', 'tokens');
        if ~isempty(tokensii)
            data.Insulin_infusion{i,1} = tokensii{1}{1}{1};             
            data.Insulin_infusion{i,2} = tokensii{1}{1}{2};              
            data.Insulin_infusion{i,3} = str2double(tokensii{1}{1}{3});  
            data.Insulin_infusion{i,4} = tokensii{1}{1}{4};
        end
    end
    
    % Glucose_concentration
    idxgc = find(contains(lines, 'Glucose_concentration'));
    start_idxgc = idxgc(1) + 3;
    cell_gc = [];
    for i = start_idxgc:length(lines)
        if strlength(strtrim(lines(i))) == 0
            break;  % space
        end
        cell_gc = [cell_gc; lines(i)];
    end
    ngc = length(cell_gc);
    data.Glucose_concentration = cell(ngc, 3);
    for i = 1:ngc   % Date   Time   conc(mmol/L)
        tokensgc = regexp(cell_gc(i), '(\d{2}/\d{2}/\d{4})\s+(\d{2}:\d{2})\s+([\d.]+)', 'tokens');
        if ~isempty(tokensgc)
            data.Glucose_concentration{i,1} = tokensgc{1}{1}{1};             
            data.Glucose_concentration{i,2} = tokensgc{1}{1}{2};              
            data.Glucose_concentration{i,3} = str2double(tokensgc{1}{1}{3});  
        end
    end
    
    % Reference_glucose_concentration
    idxrgc = find(contains(lines, 'Reference_glucose_concentration'));
    start_idxrgc = idxrgc(1) + 3;
    cell_rgc = [];
    for i = start_idxrgc:length(lines)
        if strlength(strtrim(lines(i))) == 0
            break;  % space
        end
        cell_rgc = [cell_rgc; lines(i)];
    end
    nrgc = length(cell_rgc);
    data.Reference_glucose_concentration = cell(nrgc, 3);
    for i = 1:nrgc   % Date   Time   conc(mmol/L)
        tokensrgc = regexp(cell_rgc(i), '(\d{2}/\d{2}/\d{4})\s+(\d{2}:\d{2})\s+([\d.]+)', 'tokens');
        if ~isempty(tokensrgc)
            data.Reference_glucose_concentration{i,1} = tokensrgc{1}{1}{1};             
            data.Reference_glucose_concentration{i,2} = tokensrgc{1}{1}{2};              
            data.Reference_glucose_concentration{i,3} = str2double(tokensrgc{1}{1}{3});  
        end
    end
    
    % Save the .mat file
    save(['cam' num '.mat'], 'data');  % 01-08
end
