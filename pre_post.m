function [pre,post] = pre_post(timedata,plasmadata,date_start,date_end,datatype)
    %inputs the raw data, time and values and then crops the data to the
    %event time interval. then grabs the mean beginning and end values for
    %the values of plasma parameters before and after the event
    
    %crop data
    [~,plasmadata,start_index,end_index] = crop(timedata,plasmadata,date_start,date_end);
    %Number of data points to average is 1.25% of the total number of data
    %points.
    %     n = round(length(plasmadata)*(0.0125));
    if strcmp(datatype,'B')
        n=5*128;
    elseif strcmp(datatype,'i')
        n = round(5.0/0.15);
    elseif strcmp(datatype,'e')
        n = round(5.0/0.03);
    end
    
    
    pre = mean(plasmadata(1:1+n,:));
    post = mean(plasmadata(end-n:end,:));
    
end

