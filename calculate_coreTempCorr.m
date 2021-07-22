function [tempCorr] = calculate_coreTempCorr(left_InnerEdge,right_InnerEdge,fpi_timedata,fpi_temppara,fpi_tempperp)
    
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    
    %Core Time limits
    coreStart = datenum(left_InnerEdge,formatIn);
    coreEnd = datenum(right_InnerEdge,formatIn);
    
    coreStart_index = find(fpi_timedata > coreStart, 1);
    coreEnd_index = find(fpi_timedata > coreEnd, 1);
    
    fpi_coreTimedata = fpi_timedata(coreStart_index:coreEnd_index,1);
    fpi_coretempparadata = fpi_temppara(coreStart_index:coreEnd_index,1);
    fpi_coretempperpdata = fpi_tempperp(coreStart_index:coreEnd_index,1);
    
    
    tempCorr = corrcoef(fpi_coretempparadata,fpi_coretempperpdata);
    tempCorr = tempCorr(1,2);
    
%     p = polyfit(fpi_coreTimedata,fpi_coreNdata,1); %parameters for linear fit
%     trendDensity = polyval(p,fpi_coreTimedata); %Trend Line
%     Sigma = std(fpi_coreNdata - trendDensity); %1 STD
end