function [pre,post] = calculate_prepostAverages(timedata,parameterdata,date_start,date_end,datatype,durationForAverage)
    %inputs the raw data, time and values and then crops the data to the
    %event time interval. then grabs the mean beginning and end values for
    %the values of plasma parameters before and after the event
    
    %crop data
    [~,parameterdata,start_index,end_index] = crop(timedata,parameterdata,date_start,date_end);

    if strcmp(datatype,'B')
        n=durationForAverage*128;
    elseif strcmp(datatype,'i')
        n = round(durationForAverage/0.15);
    elseif strcmp(datatype,'e')
        n = round(durationForAverage/0.03);
    end
    
    pre = mean(parameterdata(1:1+n,:));
    post = mean(parameterdata(end-n:end,:));
    
end

%Number of data points to average is 1.25% of the total number of data
%points.
%     n = round(length(plasmadata)*(0.0125));