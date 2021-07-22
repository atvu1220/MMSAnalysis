function [substructurePresent,fpi_coreTimedata,fpi_coreNdata,trendDensity,Sigma,substructureStarts,substructureEnds,...
        maxNcoreratio,CoreMaxDeflectionAngle,maxVcoreratio,maxDPcoreratio,maxTcoreratio,nVcorr,nBcorr,...
        maxNSWratio,SWMaxDeflectionAngle,maxVSWratio,maxDPSWratio,maxTSWratio,...
        maxBcoreratio,maxBSWratio,...
        SSSize,SSSize_in_ElectronScales,SS_speed,SS_duration,durationRatio,core_ion_gyroradius,SSstart_indexFraction,SSSizeBulkFlow,SSSizeBulkVflow_ionGyro,...
        CoreSWDeflectionAngle, maxTecoreratio,maxTeSWRatio,...
        SS_i_nTCorr,SS_e_nTCorr] = ...
        calculate_densityStats(fpi_timedata,...
        fpi_ndata,fpi_vdata,fpi_i_tparadata,fpi_i_tperpdata,...
        mms1_fgm_timedata_srvy,mms1_fgm_bdata_srvy,fpi_e_timedata,fpi_e_ndata,fpi_e_tperpdata,fpi_e_tparadata,...
        event_start,event_end,left_OuterEdge,right_OuterEdge,left_InnerEdge,right_InnerEdge,ssSTDs,minSSDataPoints,...
        mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
        mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
        mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
        mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
        mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
        mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
        mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
        mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
        downstreamV,upstreamV,downstreamN,upstreamN)
    
    
    %Core Time limits
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    coreStart = datenum(left_InnerEdge,formatIn);
    coreEnd = datenum(right_InnerEdge,formatIn);
    
    %Interpolate B data to FPI time cadence
    [~,fgm_Bdata] = interpxyz(mms1_fgm_timedata_srvy,mms1_fgm_bdata_srvy,fpi_timedata);
    %     [~,e_tperpdata] = interpxyz(fpi_e_timedata,fpi_e_tperpdata,fpi_timedata);
    
    coreStart_index = find(fpi_timedata > coreStart, 1);
    coreEnd_index = find(fpi_timedata > coreEnd, 1);
    
    fpi_coreTimedata = fpi_timedata(coreStart_index:coreEnd_index,1);
    fpi_coreNdata = fpi_ndata(coreStart_index:coreEnd_index,1);
    fpi_coreVdata = fpi_vdata(coreStart_index:coreEnd_index,:);
    fpi_coreTparadata = fpi_i_tparadata(coreStart_index:coreEnd_index,1);
    fpi_coreTperpdata = fpi_i_tperpdata(coreStart_index:coreEnd_index,1);
    fpi_coreTempdata = 1/3 .* (fpi_coreTparadata + 2.*fpi_coreTperpdata);
    fgm_coreBdata = fgm_Bdata(coreStart_index:coreEnd_index,:);
    
    p = polyfit(fpi_coreTimedata,fpi_coreNdata,1); %parameters for linear fit
    trendDensity = polyval(p,fpi_coreTimedata); %Trend Line
    Sigma = std(fpi_coreNdata - trendDensity); %1 STD
    
    %[downstreamN,upstreamN]=pre_post(fpi_timedata,fpi_ndata,event_start,event_end,'i');
    %[downstreamDP,upstreamDP]=pre_post(fpi_timedata,calculate_dynamic_pressure(fpi_ndata,vecnorm(fpi_vdata,2,2),'i'),event_start,event_end,'i');
    %[downstreamV,upstreamV]=pre_post(fpi_timedata,fpi_vdata,event_start,event_end,'i');
    
    %April 2020, used results from SPEDAS true solar wind velocity and density for each event
    downstreamDP = calculate_dynamic_pressure(downstreamN,downstreamV,'i');
    upstreamDP = calculate_dynamic_pressure(upstreamN,upstreamV,'i');
    
    [downstreamT,upstreamT]=pre_post(fpi_timedata, 1/3 .* (fpi_i_tparadata + 2.*fpi_i_tperpdata),event_start,event_end,'i');
    [downstreamB,upstreamB]=pre_post(fpi_timedata,fgm_Bdata,event_start,event_end,'i');
    [downstreamTe,upstreamTe]=pre_post(fpi_e_timedata,1/3 .* (fpi_e_tparadata + 2.*fpi_e_tperpdata),event_start,event_end,'e');
    
    
    
    
    %if density dips below requirement for 1 point, take those three points out of the transition matrix,
    %so that we can have a longer substructure.
    transitions = diff([0;((fpi_coreNdata >= ssSTDs*Sigma+trendDensity))~=0;0]);
    
    
    substructureStarts=find(transitions==-1);
    substructureEnds = find(transitions==1);
    substructureLengths = substructureStarts-substructureEnds;
    substructureStarts = substructureStarts(abs(substructureLengths) >= minSSDataPoints);
    substructureEnds = substructureEnds(abs(substructureLengths) >= minSSDataPoints);
    %substructureLengths(abs(substructureLengths) >= minSSDataPoints);
    substructureStarts_index = [];
    substructureEnds_index = [];
    
    if ~isempty(substructureStarts)
        substructurePresent = 1;
        
        % FIRST CHECK OF MAX SS -If more than one SS, choose the one with the largest density
        if length(substructureStarts) > 1
            maxNofEachSS = [];
            
            for i=1:length(substructureStarts)
                if ~(substructureStarts(i) >= length(fpi_coreTimedata)) && ~(substructureEnds(i) > length(fpi_coreTimedata))
                    
                    maxNofEachSS = [maxNofEachSS,max(...
                        fpi_coreNdata(min([substructureStarts(i),substructureEnds(i)])...
                        :max([substructureStarts(i),substructureEnds(i)])))];
                else
                    %This "substructure" is not within bounds
                end
            end
            
            [~,indexofMaxNofSS] = max(maxNofEachSS);
            substructureStarts_index = min([substructureStarts(indexofMaxNofSS),substructureEnds(indexofMaxNofSS)]);
            substructureEnds_index = max([substructureStarts(indexofMaxNofSS),substructureEnds(indexofMaxNofSS)]);
            
        else
            indexofMaxNofSS = 1;
            substructureStarts_index = min([substructureStarts(indexofMaxNofSS),substructureEnds(indexofMaxNofSS)]);
            substructureEnds_index = max([substructureStarts(indexofMaxNofSS),substructureEnds(indexofMaxNofSS)]);
            
        end
        
        
        
        
        %         for i=1:20
        %             if substructureStarts_index-i >= 1
        %                 if fpi_coreNdata(substructureStarts_index-i) <= trendDensity(substructureStarts_index-i)
        %                     substructureStarts_index = substructureStarts_index-i;
        %                     break;
        %                 else
        %                     continue;
        %                 end
        %             else
        %                 break
        %             end
        %         end
        %
        %         for i=0:20
        %             if substructureEnds_index+i <= length(fpi_coreNdata)
        %                 if fpi_coreNdata(substructureEnds_index+i) <= trendDensity(substructureEnds_index+i)
        %                     substructureEnds_index = substructureEnds_index+i;
        %                     break;
        %                 else
        %                     continue;
        %                 end
        %             else
        %                 break;
        %             end
        %
        %         end
        
        %SECOND CHECK of MAX SS - If indices of "substructure" are greater than core length, then it was incorrectedly labeled as
        %a SS
        if (substructureStarts_index >= length(fpi_coreTimedata)) || (substructureEnds_index > length(fpi_coreTimedata))
            substructurePresent = 0;
            maxNcoreratio = {[]};
            maxNSWratio = {[]};
            CoreMaxDeflectionAngle = {[]};
            SWMaxDeflectionAngle = {[]};
            maxVcoreratio = {[]};
            maxVSWratio = {[]};
            maxDPcoreratio = {[]};
            maxDPSWratio = {[]};
            maxTcoreratio = {[]};
            maxTSWratio = {[]};
            nVcorr = {[]};
            nBcorr = {[]};
            maxBcoreratio=  {[]};
            maxBSWratio =  {[]};
            SSSize =  {[]};
            SSSize_in_ElectronScales =  {[]};
            SS_speed = {[]};
            SS_duration = {[]};
            durationRatio = {[]};
            core_ion_gyroradius = {[]};
            SSstart_indexFraction = {[]};
            SSSizeBulkFlow = {[]};
            SSSizeBulkVflow_ionGyro = {[]};
            meanSWV=mean([downstreamV;upstreamV],1);
            CoreSWDeflectionAngle = angle(mean(fpi_coreVdata),meanSWV);
            maxTecoreratio = {[]};
            maxTeSWRatio = {[]};
            SS_i_nTCorr= {[]};
            SS_e_nTCorr =  {[]};
            
            
        else %If it is still within the bounds of the core, proceed
            
            
            %Calculate Dynamic Pressure
            [dynamic_pressure] = calculate_dynamic_pressure(fpi_coreNdata,fpi_coreVdata,'i');
            
            %Density Enhancement
            [~,~,maxNcoreratio,maxSSN,~] = calculate_ratioCore(fpi_coreNdata,substructureStarts_index,substructureEnds_index);
            maxNSWratio = maxSSN/mean([downstreamN,upstreamN]);
            
            %Magnetic Field Strength
            [~,~,maxBcoreratio,maxSSB,~] = calculate_ratioCore(fgm_coreBdata(:,4),substructureStarts_index,substructureEnds_index);
            maxBSWratio = maxSSB/(mean([downstreamB(4);upstreamB(4)]));
            
            %Velocity Deflection
            CoreMaxDeflectionAngle=calculate_deflectionAngle(fpi_coreVdata,substructureStarts_index,substructureEnds_index);
            meanSWV=mean([downstreamV;upstreamV],1);
            SWMaxDeflectionAngle=calculate_SWdeflectionAngle(meanSWV,fpi_coreVdata,substructureStarts_index,substructureEnds_index);
            
            %Velocity Magnitude
            [~,~,maxVcoreratio,maxSSV,~] = calculate_ratioCore(abs(fpi_coreVdata(:,1)),substructureStarts_index,substructureEnds_index);
            maxVSWratio = maxSSV/abs(meanSWV(1));
            
            %Dynamic Pressure
            [~,~,maxDPcoreratio,maxSSDP,~] = calculate_ratioCore(dynamic_pressure,substructureStarts_index,substructureEnds_index);
            maxDPSWratio=maxSSDP/mean([downstreamDP,upstreamDP]);
            
            %Ion Temperature
            [~,~,maxTcoreratio,maxSST,~] = calculate_ratioCore(fpi_coreTempdata,substructureStarts_index,substructureEnds_index);
            maxTSWratio=maxSST/mean([downstreamT,upstreamT]);
            
            %Convert ion SS indices to electron SS indices
            substructureStarts_index_e = find(fpi_coreTimedata(substructureStarts_index) >= fpi_e_timedata,1,'last');
            substructureEnds_index_e = find(fpi_coreTimedata(substructureEnds_index) <= fpi_e_timedata,1,'first');
            
            [~,~,maxTecoreratio,maxSSTe,~] = calculate_ratioCore(1/3.*(2.*fpi_e_tperpdata + fpi_e_tparadata),substructureStarts_index_e,substructureEnds_index_e);
            
            meanSWTe = mean([downstreamTe;upstreamTe]);
            maxTeSWRatio = maxSSTe/meanSWTe;
            
            
            %NVB correlation
            nVcorr = calculate_CorrelationCore(fpi_coreNdata,abs(vecnorm(fpi_coreVdata,2,2)),substructureStarts_index,substructureEnds_index);
            nBcorr = calculate_CorrelationCore(fpi_coreNdata,fgm_coreBdata(:,4),substructureStarts_index,substructureEnds_index);
            
            
            %Core,SW Velocity Deflection
            CoreSWDeflectionAngle = angle(mean([fpi_coreVdata(1:substructureStarts_index-1,:);fpi_coreVdata(substructureEnds_index+1:end,:)]),meanSWV);
            
            %%
            %Calculate 2-D normal of substructure using the Timing method
            
            %B and Position Data Preparation
            %%%%%%%%%%Interpolate position and bdata with bdata from mms1 %%%%%%%%%%%%%
            %Interpolate the Data
            
            [~,mms1_mec_rdata_interp] = interpxyz(mms1_mec_timedata_raw,mms1_mec_rdata_raw,mms1_fgm_timedata_raw);
            
            [~,mms2_fgm_bdata_interp] = interpxyz(mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,mms1_fgm_timedata_raw);
            [~,mms2_mec_rdata_interp] = interpxyz(mms2_mec_timedata_raw,mms2_mec_rdata_raw,mms1_fgm_timedata_raw);
            
            [~,mms3_fgm_bdata_interp] = interpxyz(mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,mms1_fgm_timedata_raw);
            [~,mms3_mec_rdata_interp] = interpxyz(mms3_mec_timedata_raw,mms3_mec_rdata_raw,mms1_fgm_timedata_raw);
            
            [~,mms4_fgm_bdata_interp] = interpxyz(mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,mms1_fgm_timedata_raw);
            [~,mms4_mec_rdata_interp] = interpxyz(mms4_mec_timedata_raw,mms4_mec_rdata_raw,mms1_fgm_timedata_raw);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            
            %%%%%%%%%%%%%%%%%Data Cropping and Matrix Formation%%%%%%%%%%%%%%%%%%%%%%%%
            %Find the start and end limits of the event in the data, index of, in
            %datenum format
            
            % %Crop mec data to specific time period
            [~,mms1_mec_rdata,~,~] = crop(mms1_fgm_timedata_raw,mms1_mec_rdata_interp,left_InnerEdge,right_InnerEdge);
            [~,mms2_mec_rdata,~,~] = crop(mms1_fgm_timedata_raw,mms2_mec_rdata_interp,left_InnerEdge,right_InnerEdge);
            [~,mms3_mec_rdata,~,~] = crop(mms1_fgm_timedata_raw,mms3_mec_rdata_interp,left_InnerEdge,right_InnerEdge);
            [~,mms4_mec_rdata,~,~] = crop(mms1_fgm_timedata_raw,mms4_mec_rdata_interp,left_InnerEdge,right_InnerEdge);
            
            %Crop fgm data to specific time period
            [mms1_fgm_timedata,mms1_fgm_bdata,~,~] = crop(mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,left_InnerEdge,right_InnerEdge);
            [mms2_fgm_timedata,mms2_fgm_bdata,~,~] = crop(mms1_fgm_timedata_raw,mms2_fgm_bdata_interp,left_InnerEdge,right_InnerEdge);
            [mms3_fgm_timedata,mms3_fgm_bdata,~,~] = crop(mms1_fgm_timedata_raw,mms3_fgm_bdata_interp,left_InnerEdge,right_InnerEdge);
            [mms4_fgm_timedata,mms4_fgm_bdata,~,~] = crop(mms1_fgm_timedata_raw,mms4_fgm_bdata_interp,left_InnerEdge,right_InnerEdge);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Find the time of max density
            maxN_FPItime = fpi_coreTimedata(fpi_coreNdata == max(fpi_coreNdata(substructureStarts_index:substructureEnds_index)));
            SSStart_FPItime = fpi_coreTimedata(substructureStarts_index,1);
            SSEnd_FPItime = fpi_coreTimedata(substructureEnds_index,1);
            
            SSStart_FGMtime = find(mms1_fgm_timedata >= SSStart_FPItime,1);
            SSEnd_FGMtime = find(mms1_fgm_timedata >= SSEnd_FPItime,1);
            
            maxN_FGMtime=find(mms1_fgm_timedata >= maxN_FPItime,1);
            
            
            t_0= maxN_FGMtime;
            data_points = (SSEnd_FGMtime-SSStart_FGMtime);
            data_points = 4*round(data_points/4);
            [n_boundary,v_boundary,ccmean,cc12, cc13, cc14, cc12max, cc13max, cc14max, time_boundary,unique_timings] = timing_method(t_0,data_points,...
                mms1_fgm_timedata,mms1_fgm_bdata(:,1:3),...
                mms2_fgm_bdata(:,1:3),...
                mms3_fgm_bdata(:,1:3),...
                mms4_fgm_bdata(:,1:3),...
                mms1_mec_rdata(:,1:3),...
                mms2_mec_rdata(:,1:3),...
                mms3_mec_rdata(:,1:3),...
                mms4_mec_rdata(:,1:3));
            
            [~,maxCorrIndex] = max(ccmean);
            n_boundary = n_boundary(:,maxCorrIndex);
            v_boundary = v_boundary(:,maxCorrIndex);
            SS_speed = v_boundary;
            SS_duration = (86400*(mms1_fgm_timedata(SSEnd_FGMtime)-mms1_fgm_timedata(SSStart_FGMtime)))
            %%SS Size
            SSSize = v_boundary*SS_duration; %in km
            
            %Constants
            mion = 1.67e-27;%kg
            me = 9.11e-31;%kg
            epsilon_0 = 8.85e-12;
            q=1.60e-19;%C
            speedofLight = 2.999e8; %m/s
            
            
            
            
            %              SSTperp = mean(fpi_e_tperpdata(find(fpi_e_timedata>...
            %                  fpi_coreTimedata(substructureStarts_index),1):find(fpi_e_timedata>fpi_coreTimedata(substructureEnds_index),1)));
            %             vperp = (2*q*SSTperp/(9.11e-31))^(1/2); %in m/s
            %
            %             electron_gyroradii = (9.11e-31)*vperp / (q*SSB)/1000; %in km
            
            % Calculate Correlations between n and Temp
            SS_n_e_values = fpi_e_ndata(find(fpi_e_timedata>...
                fpi_coreTimedata(substructureStarts_index),1):find(fpi_e_timedata>fpi_coreTimedata(substructureEnds_index)));
            SS_temp_e_values = (1/3).*(fpi_e_tparadata(find(fpi_e_timedata>...
                fpi_coreTimedata(substructureStarts_index),1):find(fpi_e_timedata>fpi_coreTimedata(substructureEnds_index))) + ...
                2.*fpi_e_tperpdata(find(fpi_e_timedata>...
                fpi_coreTimedata(substructureStarts_index),1):find(fpi_e_timedata>fpi_coreTimedata(substructureEnds_index))));
            SS_e_nTCorr = corrcoef(SS_n_e_values,SS_temp_e_values);
            SS_e_nTCorr = SS_e_nTCorr(1,2);
            
            SS_n_i_values = fpi_coreNdata(substructureStarts_index:substructureEnds_index);
            SS_temp_i_values = fpi_coreTempdata(substructureStarts_index:substructureEnds_index);
            SS_i_nTCorr = corrcoef(SS_n_i_values,SS_temp_i_values);
            SS_i_nTCorr = SS_i_nTCorr(1,2);
            
            %%Electron Scale Calculation
            SS_n_e = 10^6 * mean(fpi_e_ndata(find(fpi_e_timedata>...
                fpi_coreTimedata(substructureStarts_index),1):find(fpi_e_timedata>fpi_coreTimedata(substructureEnds_index),1)));
            omega_pe = ( SS_n_e * q^2 / (me*epsilon_0) ) ^(1/2)
            electron_skinDepth =  speedofLight/omega_pe / 1000 %in km
            
            SSSize_in_ElectronScales = SSSize/electron_skinDepth
            
            durationRatio = SS_duration / (86400*(fpi_coreTimedata(end)-fpi_coreTimedata(1)))
            
            %Event Core Size
            CoreB = mean(fgm_coreBdata(:,4))*10^-9;
            CoreTperp = mean(fpi_coreTperpdata);
            core_vperp = (2*q*CoreTperp/mion)^(1/2);
            core_ion_gyroradius = mion*core_vperp / (q*CoreB)/1000
            SSstart_indexFraction = floor(mean([substructureStarts_index,substructureEnds_index])) / length(fpi_coreNdata)
            
            %SS Size using Ion Bulk Flow
            SSSizeBulkFlow = SS_duration*mean(vecnorm(fpi_coreVdata(substructureStarts_index:substructureEnds_index,:),2,2)) %Size in km
            
            SSB = mean(mms1_fgm_bdata(SSStart_FGMtime:SSEnd_FGMtime,4))*10^-9;
            SSTperp = mean(fpi_coreTperpdata(substructureStarts_index:substructureEnds_index)); %Ion
            vperp = (2*q*SSTperp/(mion))^(1/2);
            ion_gyroradius = (mion)*vperp / (q*SSB)/1000; %in km
            SSSizeBulkVflow_ionGyro= SSSizeBulkFlow/ion_gyroradius
            %%
        end
        
    else
        substructurePresent = 0;
        maxNcoreratio = {[]};
        maxNSWratio = {[]};
        CoreMaxDeflectionAngle = {[]};
        SWMaxDeflectionAngle = {[]};
        maxVcoreratio = {[]};
        maxVSWratio = {[]};
        maxDPcoreratio = {[]};
        maxDPSWratio = {[]};
        maxTcoreratio = {[]};
        maxTSWratio = {[]};
        nVcorr = {[]};
        nBcorr = {[]};
        maxBcoreratio=  {[]};
        maxBSWratio =  {[]};
        SSSize =  {[]};
        SSSize_in_ElectronScales =  {[]};
        SS_speed = {[]};
        SS_duration = {[]};
        durationRatio = {[]};
        core_ion_gyroradius = {[]};
        SSstart_indexFraction = {[]};
        SSSizeBulkFlow = {[]};
        SSSizeBulkVflow_ionGyro = {[]};
        meanSWV=mean([downstreamV;upstreamV],1);
        CoreSWDeflectionAngle = angle(mean(fpi_coreVdata),meanSWV);
        maxTecoreratio = {[]};
        maxTeSWRatio = {[]};
        SS_i_nTCorr= {[]};
        SS_e_nTCorr =  {[]};
    end
    
    
    
    
end
function [maxDeflectionAngle] = calculate_SWdeflectionAngle(SWMeanVelocity,CoreVelocityVector,substructureStarts_index,substructureEnds_index)
    SSpoints = CoreVelocityVector(substructureStarts_index:substructureEnds_index,:);
    SSdeflectionAngles = [];
    for i=1:length(SSpoints)
        SSdeflectionAngles = [SSdeflectionAngles,angle(SSpoints(i,:),SWMeanVelocity)];
    end
    maxDeflectionAngle = max(SSdeflectionAngles);
    
end
function [maxDeflectionAngle] = calculate_deflectionAngle(CoreVelocityVector,substructureStarts_index,substructureEnds_index)
    coreVMinusSS = mean([CoreVelocityVector(1:substructureStarts_index-1,:);CoreVelocityVector(substructureEnds_index+1:end,:)]);
    SSpoints = CoreVelocityVector(substructureStarts_index:substructureEnds_index,:);
    SSdeflectionAngles = [];
    for i=1:length(SSpoints)
        SSdeflectionAngles = [SSdeflectionAngles,angle(SSpoints(i,:),coreVMinusSS)];
    end
    maxDeflectionAngle = max(SSdeflectionAngles);
    
end
function [corr] = calculate_CorrelationCore(parameterCoredata1,parameterCoredata2,substructureStarts_index,substructureEnds_index)
    corrMatrix = corrcoef(parameterCoredata1(substructureStarts_index:substructureEnds_index),parameterCoredata2(substructureStarts_index:substructureEnds_index));
    corr = corrMatrix(1,2);
    
end
function [meanRatio,minRatio,maxRatio,maxSS,meanSS] = calculate_ratioCore(parameterCoredata,substructureStarts_index,substructureEnds_index)
    coreParameterMinusSS = [parameterCoredata(1:substructureStarts_index-1);parameterCoredata(substructureEnds_index+1:end)];
    meancoreParameterMinusSS = mean(coreParameterMinusSS);
    maxSS = max(parameterCoredata(substructureStarts_index:substructureEnds_index));
    minSS = min(parameterCoredata(substructureStarts_index:substructureEnds_index));
    meanSS = mean(parameterCoredata(substructureStarts_index:substructureEnds_index));
    maxRatio = maxSS/meancoreParameterMinusSS;
    minRatio = minSS/meancoreParameterMinusSS;
    meanRatio = meanSS/meancoreParameterMinusSS;
end
function [] = calculate_SStiming()
    %Take input B at srvy resolution
    %Using the FPI timestamps for the beginning and end of the SS, find the equivalent timestamps
    %for srvy B
    %Perform Timing method on
    
    
    
    
    
end