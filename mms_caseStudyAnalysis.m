%Performs basic calculations for Case Study Events

%CS1
event_start= '2019-02-26 00:32:30.000';
event_end = '2019-02-26 00:33:40.000';
leading_leftmost_date = '2019-02-26 00:30:45.000'
leading_rightmost_date = '2019-02-26 00:31:00.000'
trailing_leftmost_date = '2019-02-26 00:33:25.000'
trailing_rightmost_date = '2019-02-26 00:33:33.000'
left_InnerEdge = '2019-02-26 00:33:00.000'
right_InnerEdge= '2019-02-26 00:33:15.000'

%CS2
event_start= '2019-02-22 12:01:00.000';
event_end = '2019-02-22 12:01:40.000';
leading_leftmost_date = '2019-02-22 11:40:00.000'
leading_rightmost_date = '2019-02-22 11:45:00.000'
trailing_leftmost_date = '2019-02-22 12:01:40.000'
trailing_rightmost_date = '2019-02-22 12:02:30.000'
left_InnerEdge = '2019-02-22 12:01:15.000'
right_InnerEdge= '2019-02-22 12:01:25.000'

upstreamMMSprobe = 3;
threshold_std=0;
%% Load Data
%Load MEC data
[mms1_mec_timedata_raw, mms1_mec_rdata_raw] = load_mec(event_start,event_end,1,data_type);
[mms2_mec_timedata_raw, mms2_mec_rdata_raw] = load_mec(event_start,event_end,2,data_type);
[mms3_mec_timedata_raw, mms3_mec_rdata_raw] = load_mec(event_start,event_end,3,data_type);
[mms4_mec_timedata_raw, mms4_mec_rdata_raw] = load_mec(event_start,event_end,4,data_type);

%load FGM data
[mms1_fgm_timedata_raw, mms1_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,1,'brst');
[mms2_fgm_timedata_raw, mms2_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,2,'brst');
[mms3_fgm_timedata_raw, mms3_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,3,'brst');
[mms4_fgm_timedata_raw, mms4_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,4,'brst');

[mmsx_fgm_timedata_srvy, mmsx_fgm_bdata_srvy] = load_fgm(event_start,event_end,upstreamMMSprobe,'srvy'); %For Sliding Window

%Load FPI data
[fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_tparadata,fpi_i_tperpdata,...
    fpi_i_edata,fpi_i_espectdata,fpi_i_pressdata] = load_fpi(event_start,event_end,upstreamMMSprobe,'brst','i');

%% Calculation of Current Sheet Normal
if threshold_std == 0
    %Current Sheet Normal Calculation Plot, manually.
    [n_cs,B_pre_cs,B_post_cs] = manualCurrentSheet(event_start,event_end,...
        leading_leftmost_date,leading_rightmost_date,...
        trailing_leftmost_date,trailing_rightmost_date,mmsx_fgm_timedata_srvy,mmsx_fgm_bdata_srvy);
else
    %Calculate Current Sheet
    [n_cs,B_pre_cs,B_post_cs] = calculateCurrentSheet(event_start,event_end,mmsx_fgm_timedata_srvy,mmsx_fgm_bdata_srvy,threshold_std);
end

%% Calculation of Bow Shock Normal
%Down/Up-stream Averages
[downstreamDensity,upstreamDensity] = calculate_prepostAverages(fpi_i_timedata,fpi_i_ndata,event_start,event_end,'i',5);
[downstreamVelocity,upstreamVelocity] = calculate_prepostAverages(fpi_i_timedata,fpi_i_vdata,event_start,event_end,'i',5);
%Mach Number
mass_ion = 1.6726219e-27; %kg
mu_o = 1.25663706e-6;

VA_pre = (norm(B_pre_cs)*1e-9/(mass_ion*downstreamDensity*1e6*mu_o)^(1/2))/1e3
M_pre = norm(downstreamVelocity)/VA_pre
VA_post = (norm(B_post_cs)*1e-9/(mass_ion*upstreamDensity*1e6*mu_o)^(1/2))/1e3
M_post = norm(upstreamVelocity)/VA_post

[shock_normal,HFAtoBS_Distance,bs_pos,closest_shock_normal,ClosestDistance,closest_point] = calculate_MerkaBS_CS2(...
    mms1_fgm_timedata_raw,...
    mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
    mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
    mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
    mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
    left_InnerEdge,right_InnerEdge,...
    n_cs,downstreamDensity,downstreamVelocity,M_pre);
closest_shock_normal
shock_normal
ClosestDistance
HFAtoBS_Distance
