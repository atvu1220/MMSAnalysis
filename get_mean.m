function [average] = get_mean(timedata,plasmadata,date_start,date_end)
    %Returns average of parameter within time intervals
    %crop data
    [~,plasmadata,~,~] = crop(timedata,plasmadata,date_start,date_end);
    
    average = mean(plasmadata(:,:));
end

