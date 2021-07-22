function [vCorr,vxCorr,tempparaCorr,tempperpCorr,tempCorr,BmagCorr] = calculate_correlation(leftTime,rightTime,fpi_timedata,fpi_ndata,fpi_vdata,fpi_tparadata,fpi_tperpdata,fgm_bdata)
    %calculates the correlation coefficient [-1,1] for two fpi parameters.
    %density with velocity mag, x-velocity, temp para, temp perp, temp, B
    
    
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    tstart = datenum(leftTime,formatIn);
    tend = datenum(rightTime,formatIn);
    
    %Find the start and end limits of the event in the data
    start_index = find(fpi_timedata > tstart, 1);
    end_index = find(fpi_timedata > tend, 1);
    
    %crop data
    fpi_timedata = fpi_timedata(start_index:end_index,1);
    fpi_ndata = fpi_ndata(start_index:end_index,:);
    fpi_vdata = fpi_vdata(start_index:end_index,:);
    fpi_vxdata = fpi_vdata(:,1);
    fpi_vdata = vecnorm(fpi_vdata,2,2);
    fpi_tparadata = fpi_tparadata(start_index:end_index,:);
    fpi_tperpdata = fpi_tperpdata(start_index:end_index,:);
    fpi_tdata = (fpi_tparadata.^2 + fpi_tperpdata.^2).^(1/2);
    
    fgm_bdata = fgm_bdata(start_index:end_index,:);
    fgm_bmagdata = vecnorm(fgm_bdata,2,2);
    
    %Calculate correlation coefficients
    vCorr = corrcoef(fpi_ndata,fpi_vdata);
    vxCorr = corrcoef(fpi_ndata,fpi_vxdata);
    tempparaCorr = corrcoef(fpi_ndata,fpi_tparadata);
    tempperpCorr = corrcoef(fpi_ndata,fpi_tperpdata);
    tempCorr = corrcoef(fpi_ndata,fpi_tdata);
    BmagCorr = corrcoef(fpi_ndata,fgm_bmagdata);
    
    vCorr = vCorr(1,2);
    vxCorr = vxCorr(1,2);
    tempparaCorr = tempparaCorr(1,2);
    tempperpCorr = tempperpCorr(1,2);
    tempCorr = tempCorr(1,2);
    BmagCorr = BmagCorr(1,2);