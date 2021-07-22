function [Sigma,CV] = calculate_coreDensitySTD(left_InnerEdge,right_InnerEdge,fpi_timedata,fpi_ndata)
    
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    
    %Core Time limits
    coreStart = datenum(left_InnerEdge,formatIn);
    coreEnd = datenum(right_InnerEdge,formatIn);
    
    coreStart_index = find(fpi_timedata > coreStart, 1);
    coreEnd_index = find(fpi_timedata > coreEnd, 1);
    
    fpi_coreTimedata = fpi_timedata(coreStart_index:coreEnd_index,1);
    fpi_coreNdata = fpi_ndata(coreStart_index:coreEnd_index,1);
    
    
    p = polyfit(fpi_coreTimedata,fpi_coreNdata,1); %parameters for linear fit
    trendDensity = polyval(p,fpi_coreTimedata); %Trend Line
    Sigma = std(fpi_coreNdata - trendDensity); %1 STD
    
    average = mean(trendDensity);
    
    CV = Sigma/average; %Coefficient of Variance, STD/AVG(TrendLine)
end