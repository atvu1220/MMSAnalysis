function [edge_nBCorr,edge_nvCorr,edge_vBCorr,edge_nTCorr,...
        edge_meanB,edge_maxB,deltaB, maxBoverAmbientB,...
        edge_meanN,edge_maxN,deltaN, maxNoverAmbientN,...
        edge_meanV,edge_maxV,deltaV, maxVoverAmbientV] = calculate_edgesAnalysis(event_start,event_end,leftTime,rightTime,fpi_timedata,fpi_ndata,fpi_vdata,fpi_tparadata,fpi_tperpdata,fgm_bdata,durationForAverage,side)
    %analyzes parametrs within event edges. and outside, 5s beyound event edges average.
    outsidePoints = floor(5/0.15);
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    
    edgeStart = datenum(leftTime,formatIn);
    edgeEnd = datenum(rightTime,formatIn);
    eventStart = datenum(event_start,formatIn);
    eventEnd = datenum(event_end,formatIn);
    
    
    %Find the start and end limits of the event in the data
    edgeStart_index = find(fpi_timedata > edgeStart, 1);
    edgeEnd_index = find(fpi_timedata > edgeEnd, 1);
    
    
    eventStart_index = find(fpi_timedata > eventStart, 1);
    eventEnd_index = find(fpi_timedata > eventEnd, 1);
    
    %crop data
    edge_ndata = fpi_ndata(edgeStart_index:edgeEnd_index,:);
    
    fpi_vdata = vecnorm(fpi_vdata,2,2);
    edge_vdata = fpi_vdata(edgeStart_index:edgeEnd_index);
    
    
    edge_tparadata = fpi_tparadata(edgeStart_index:edgeEnd_index,:);
    edge_tperpdata = fpi_tperpdata(edgeStart_index:edgeEnd_index,:);
    edge_tdata = (edge_tparadata.^2 + edge_tperpdata.^2).^(1/2);
    
    fgm_bdata = vecnorm(fgm_bdata,2,2);
    edge_bdata = fgm_bdata(edgeStart_index:edgeEnd_index,:);
    
    
    %Calculations
    %Correlations
    edge_nBCorr = corrcoef(edge_ndata,edge_bdata);
    edge_nBCorr = edge_nBCorr(1,2);
    
    edge_nvCorr = corrcoef(edge_ndata,edge_vdata);
    edge_nvCorr = edge_nvCorr(1,2);
    
    edge_vBCorr = corrcoef(edge_vdata,edge_bdata);
    edge_vBCorr = edge_vBCorr(1,2);
    
    edge_nTCorr = corrcoef(edge_ndata,edge_tdata);
    edge_nTCorr = edge_nTCorr(1,2);
    
    %Max Means inside edge and also outside mean.
    edge_meanB = mean(edge_bdata);
    edge_maxB = max(edge_bdata);
    
    edge_meanN = mean(edge_ndata);
    edge_maxN = max(edge_ndata);
    
    edge_meanV = mean(edge_vdata);
    edge_maxV = max(edge_vdata);
    
    
    %Find outside parameters, consistent with other code for calculate_criteria
    
    if strcmp(side,'leading')
        outside_meanB = mean(fgm_bdata(eventStart_index:edgeStart_index));
        outside_meanN = mean(fpi_ndata(eventStart_index:edgeStart_index));
        outside_meanV = mean(fpi_vdata(eventStart_index:edgeStart_index));
    elseif  strcmp(side,'trailing')
        outside_meanB = mean(fgm_bdata(edgeEnd_index:eventEnd_index));
        outside_meanN = mean(fpi_ndata(edgeEnd_index:eventEnd_index));
        outside_meanV = mean(fpi_vdata(edgeEnd_index:eventEnd_index));
    end
    
    [outside_meanB_left,outside_meanB_right] = calculate_prepostAverages(fpi_timedata,fgm_bdata,event_start,event_end,'i',durationForAverage);
    [outside_meanN_left,outside_meanN_right] = calculate_prepostAverages(fpi_timedata,fpi_ndata,event_start,event_end,'i',durationForAverage);
    [outside_meanV_left,outside_meanV_right] = calculate_prepostAverages(fpi_timedata,fpi_vdata,event_start,event_end,'i',durationForAverage);
%     
    % % %     %Make sure it isn't out of bounds
    % % %     if edgeStart_index-outsidePoints < 1
    % % %         outside_meanB_left = mean(fgm_bdata(1:edgeStart_index,:));
    % % %         outside_meanN_left = mean(fpi_ndata(1:edgeStart_index,:));
    % % %         outside_meanV_left = mean(fpi_vdata(1:edgeStart_index,:));
    % % %     else
    % % %         outside_meanB_left = mean(fgm_bdata(edgeStart_index-outsidePoints:edgeStart_index,:));
    % % %         outside_meanN_left = mean(fpi_ndata(edgeStart_index-outsidePoints:edgeStart_index,:));
    % % %         outside_meanV_left = mean(fpi_vdata(edgeStart_index-outsidePoints:edgeStart_index,:));
    % % %     end
    % % %
    % % %
    % % %     if edgeEnd_index+outsidePoints > length(fgm_bdata)
    % % %         outside_meanB_right = mean(fgm_bdata(edgeEnd_index:end,:));
    % % %         outside_meanN_right = mean(fpi_ndata(edgeEnd_index:end,:));
    % % %         outside_meanV_right = mean(fpi_vdata(edgeEnd_index:edgeEnd_index,:));
    % % %     else
    % % %         outside_meanB_right = mean(fgm_bdata(edgeEnd_index:edgeEnd_index+outsidePoints,:));
    % % %         outside_meanN_right = mean(fpi_ndata(edgeEnd_index:edgeEnd_index+outsidePoints,:));
    % % %         outside_meanV_right = mean(fpi_vdata(edgeEnd_index:edgeEnd_index+outsidePoints,:));
    % % %     end
    
    deltaB = abs(outside_meanB_left - outside_meanB_right);
%     meandeltaB = mean([outside_meanB_left,outside_meanB_right]);
    %maxBovermeandeltaB = edge_maxB/meandeltaB;
    
    
    deltaN = abs(outside_meanN_left - outside_meanN_right);
%     meandeltaN = mean([outside_meanN_left,outside_meanN_right]);
    %maxNovermeandeltaN = edge_maxN/meandeltaN;
    
    
    deltaV = abs(outside_meanV_left - outside_meanV_right);
%     meandeltaV = mean([outside_meanV_left,outside_meanV_right]);
    %maxVovermeandeltaV = edge_maxV/meandeltaV;
    maxBoverAmbientB = edge_maxB/outside_meanB;
    edge_maxN
    outside_meanN
    maxNoverAmbientN = edge_maxN/outside_meanN
    maxVoverAmbientV = edge_maxV/outside_meanV;
        
%         
%     if strcmp(side,'leading')
%         maxBoverAmbientB = edge_maxB/outside_meanB_left;
%         maxNoverAmbientN = edge_maxN/outside_meanN_left;
%         maxVoverAmbientV = edge_maxV/outside_meanV_left;
%     elseif strcmp(side,'trailing')
%         maxBoverAmbientB = edge_maxB/outside_meanB_right;
%         maxNoverAmbientN = edge_maxN/outside_meanN_right;
%         maxVoverAmbientV = edge_maxV/outside_meanV_right;
%     end