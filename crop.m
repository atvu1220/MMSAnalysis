function [timedata,plasmadata,start_index,end_index] = crop(timedata,plasmadata,event_start,event_end)
    
    %crops data to within specified time range
    if isduration(event_start) == 1
        formatOut='yyyy-mm-dd HH:MM:SS.FFF';
        event_start = datestr(event_start,formatOut);
        event_end = datestr(event_end,formatOut);
    end
    
    
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    
    tstart = datenum(event_start,formatIn);
    tend = datenum(event_end,formatIn);
    
    % %Find the start and end limits of the event in the data
    start_index = find(timedata >= tstart, 1);
    end_index = find(timedata >= tend, 1);
    
    %Convert the datetime to date String, and then crop to our event timeframe
    timedata = timedata(start_index:end_index,1);
    plasmadata = plasmadata(start_index:end_index,:);
end

