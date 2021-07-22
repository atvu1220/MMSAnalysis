function [core_BmagCorr,core_vxCorr,core_vmagCorr,core_tparaCorr,core_tperpCorr] = calculate_coreCorr(left_InnerEdge,right_InnerEdge,fgm_timedata,fgm_bdata,fpi_timedata,fpi_ndata,fpi_vdata,fpi_temppara,fpi_tempperp)
    
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    %Downsample FGM to FPi
    [~,fgm_bdata_interp] = interpxyz(fgm_timedata,fgm_bdata(:,1:3),fpi_timedata);

    %Core Time limits
    coreStart = datenum(left_InnerEdge,formatIn);
    coreEnd = datenum(right_InnerEdge,formatIn);
    
    coreStart_index = find(fpi_timedata > coreStart, 1);
    coreEnd_index = find(fpi_timedata > coreEnd, 1);
    
    
    fgm_coreBdata = fgm_bdata_interp(coreStart_index:coreEnd_index,1:3);
    fgm_coreBdata_mag = vecnorm(fgm_coreBdata,2,2);
    
    fpi_coreTimedata = fpi_timedata(coreStart_index:coreEnd_index,1);
    fpi_corendata = fpi_ndata(coreStart_index:coreEnd_index,1);
    fpi_corevdata = fpi_vdata(coreStart_index:coreEnd_index,1:3);
    fpi_corevdata_mag = vecnorm(fpi_vdata,2,2);
    fpi_coretempparadata = fpi_temppara(coreStart_index:coreEnd_index,1);
    fpi_coretempperpdata = fpi_tempperp(coreStart_index:coreEnd_index,1);
    
    
    core_BmagCorr = corrcoef(fpi_corendata,fgm_coreBdata_mag);
    core_BmagCorr = core_BmagCorr(1,2);
    
    core_vxCorr = corrcoef(fpi_corendata,fpi_corevdata(:,1));
    core_vxCorr = core_vxCorr(1,1);
    core_vmagCorr = corrcoef(fpi_corendata,fpi_corevdata_mag);
    core_vmagCorr = core_vmagCorr(1,2);
    core_tparaCorr = corrcoef(fpi_corendata,fpi_coretempparadata);
    core_tparaCorr = core_tparaCorr(1,2);
    core_tperpCorr = corrcoef(fpi_corendata,fpi_coretempperpdata);
    core_tperpCorr = core_tperpCorr(1,2);
   
%     tempCorr = corrcoef(fpi_coretempparadata,fpi_coretempperpdata);
%     tempCorr = tempCorr(1,2);
    
%     p = polyfit(fpi_coreTimedata,fpi_coreNdata,1); %parameters for linear fit
%     trendDensity = polyval(p,fpi_coreTimedata); %Trend Line
%     Sigma = std(fpi_coreNdata - trendDensity); %1 STD
end