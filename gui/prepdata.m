function prepdata(filename)
% Remove duplicate rows
    dat = readtable(filename);
    
    varlist = [dat.T_oven dat.T_Cu1 dat.T_Cu2 dat.T_Inco1 dat.T_H25 dat.T_Inco2 dat.T_Cu3];
    filter = zeros(size(varlist));
    
    % pick all cells where entry same as one below
    for i = 1:size(varlist, 2) % columns in varlist
        for j = 1:size(varlist, 1)-1 % rows -1
            filter(j, i) = varlist(j, i) == varlist(j+1, i);
        end
    end
    
    % pick all rows where all entries not same as row below
    filter = not(all(filter, 2));
    
    dat_filt = dat(filter, :);
    
    writetable(dat_filt, filename)
end