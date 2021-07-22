%Andrew Vu 9/25/18
%version 2 has fixed the problem of the time interval of the "centered" mva
%and in addition, calculates the total G matrix and averages it over the
%interval -- rather than the average of the normal vectors.
%data=spdfcdfread(filename);
%datainfo=spdfcdfinfo(filename);
%clear
%%
addpath ~/Library/'Mobile Documents'/com~apple~CloudDocs/
addpath ~/data/

tic
close all
clear
figure('Position',[0 0 800 700])
movegui('northwest')
mms_directory = '/Users/andrewvu/data/mms/';
probe_num = '1';
num_plots = 8;
data_type = 'brst';


%%%%%%%%%%%%%%%%%%%%%%%%%%Event Date Input%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Enter 0 for placeholders for dates.
%Overall Event Start and End Times
event_start = '2018-03-12 07:35:50.000';
event_end = '2018-03-12 07:36:30.000';

%Specific Portion of the Event for Analysis
date_start = '2018-03-12 07:36:07.500';
date_end = '2018-03-12 07:36:8.500';


event_start = '2018-01-09 08:34:27.000';
event_end = '2018-01-09 08:35:03.000';

date_start = '2018-01-09 08:34:49.500';
date_end = '2018-01-09 08:34:51.500';


event_start = '2018-04-01 01:10:34.000';
event_end = '2018-04-01 01:11:22.000';

date_start = '2018-04-01 01:10:57.500';
date_end = '2018-04-01 01:10:58.500';


event_start = '2018-04-27 19:47:30.000';
event_end = '2018-04-27 19:48:15.000';

date_start = '2018-04-27 19:47:47.250';
date_end = '2018-04-27 19:47:49.000';


event_start = '2018-03-01 01:03:45.000';
event_end = '2018-03-01 01:04:20.000';

date_start = '2018-03-01 01:03:57.000';
date_end = '2018-03-01 01:03:58.250';



event_start = '2018-01-09 08:34:27.000';
event_end = '2018-01-09 08:35:03.000';

date_start = '2018-01-09 08:34:49.250';
date_end = '2018-01-09 8:34:51.000';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Adjustable Parameters%%%%%%%%%%%%%%%%%%%%%%%%%
avg_i = 0; %for removing fluctuations from calibrationb errors in G
fluctuate_mdd = 0;%1 for force all vectors to stay within their same sign
fluctuate_mva = 0;
MVAB_window_size = 0; %this is the number of data points on either side of t_0; enter 0 for default window size of the time interval

%Plot the minimum, intermediate, or maximum vectors from eigenvalue problem
%1 is max, 2 is mid, 3 is min
mddmm = 1;
mvamm = 3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%Center +/- sides Method%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %center of discontinuity
% center_time = '2018-03-12 07:36:15.000';
% %dt is window size, on each side.
% dt = 1;
% %Create left and Right Initial Bounds
% formatIn='yyyy-MM-dd HH:mm:ss.SSS';
% center_datetime = datetime(center_time,'InputFormat',formatIn);
% date_start = (center_datetime - seconds(dt));
% date_end = (center_datetime + seconds(dt));
%
% formatIn='yyyy-mm-dd HH:MM:SS.FFF';
% date_start = datestr(date_start,formatIn);
% date_end = datestr(date_end,formatIn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Set the Strings for which eigenvector to use
if mddmm==1
    mddl = '{max}';
elseif mddmm==2
    mddl = '{mid}';
elseif mddmm==3
    mddl = '{min}';
end

if mvamm==1
    mval = '{max}';
elseif mvamm ==2
    mval = '{mid}';
elseif mvamm==3
    mval = '{min}';
end

MVA_normal = [0,0,0];

formatIn='yyyy-mm-dd HH:MM:SS.FFF';
tstart = datenum(date_start,formatIn);
tend = datenum(date_end,formatIn);
%%
%%%%%%%%%%%%%%%%%%%%%%%%%Load FGM/Mec Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[mms1_fgm_btimedata_raw, mms1_fgm_bdata_raw, mms1_fgm_rtimedata_raw, mms1_fgm_rdata_raw] = load_fgm(date_start,1,data_type);
[mms2_fgm_btimedata_raw, mms2_fgm_bdata_raw, mms2_fgm_rtimedata_raw, mms2_fgm_rdata_raw] = load_fgm(date_start,2,data_type);
[mms3_fgm_btimedata_raw, mms3_fgm_bdata_raw, mms3_fgm_rtimedata_raw, mms3_fgm_rdata_raw] = load_fgm(date_start,3,data_type);
[mms4_fgm_btimedata_raw, mms4_fgm_bdata_raw, mms4_fgm_rtimedata_raw, mms4_fgm_rdata_raw] = load_fgm(date_start,4,data_type);

%Mec data has fewer data points than FGM Epheremis
[mms1_mec_rtimedata_raw, mms1_mec_rdata_raw] = load_mec(date_start,1,'srvy');
[mms2_mec_rtimedata_raw, mms2_mec_rdata_raw] = load_mec(date_start,2,'srvy');
[mms3_mec_rtimedata_raw, mms3_mec_rdata_raw] = load_mec(date_start,3,'srvy');
[mms4_mec_rtimedata_raw, mms4_mec_rdata_raw] = load_mec(date_start,4,'srvy');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%Interpolate position and bdata with bdata from mms1 %%%%%%%%%%%%%
%Interpolate the Data
% mms1_fgm_rdata_interp(:,1) = interp1(mms1_fgm_rtimedata_raw,mms1_fgm_rdata_raw(:,1),mms1_fgm_btimedata_raw);
% mms1_fgm_rdata_interp(:,2) = interp1(mms1_fgm_rtimedata_raw,mms1_fgm_rdata_raw(:,2),mms1_fgm_btimedata_raw);
% mms1_fgm_rdata_interp(:,3) = interp1(mms1_fgm_rtimedata_raw,mms1_fgm_rdata_raw(:,3),mms1_fgm_btimedata_raw);

[~,mms1_mec_rdata_interp] = interpxyz(mms1_mec_rtimedata_raw,mms1_mec_rdata_raw,mms1_fgm_btimedata_raw);
% mms1_mec_rdata_interp(:,1) = interp1(mms1_mec_rtimedata_raw,mms1_mec_rdata_raw(:,1),mms1_fgm_btimedata_raw);
% mms1_mec_rdata_interp(:,2) = interp1(mms1_mec_rtimedata_raw,mms1_mec_rdata_raw(:,2),mms1_fgm_btimedata_raw);
% mms1_mec_rdata_interp(:,3) = interp1(mms1_mec_rtimedata_raw,mms1_mec_rdata_raw(:,3),mms1_fgm_btimedata_raw);

[~,mms2_fgm_bdata_interp] = interpxyz(mms2_fgm_btimedata_raw,mms2_fgm_bdata_raw,mms1_fgm_btimedata_raw);
% mms2_fgm_bdata_interp(:,1) = interp1(mms2_fgm_btimedata_raw,mms2_fgm_bdata_raw(:,1),mms1_fgm_btimedata_raw);
% mms2_fgm_bdata_interp(:,2) = interp1(mms2_fgm_btimedata_raw,mms2_fgm_bdata_raw(:,2),mms1_fgm_btimedata_raw);
% mms2_fgm_bdata_interp(:,3) = interp1(mms2_fgm_btimedata_raw,mms2_fgm_bdata_raw(:,3),mms1_fgm_btimedata_raw);

% mms2_fgm_rdata_interp(:,1) = interp1(mms2_fgm_rtimedata_raw,mms2_fgm_rdata_raw(:,1),mms1_fgm_btimedata_raw);
% mms2_fgm_rdata_interp(:,2) = interp1(mms2_fgm_rtimedata_raw,mms2_fgm_rdata_raw(:,2),mms1_fgm_btimedata_raw);
% mms2_fgm_rdata_interp(:,3) = interp1(mms2_fgm_rtimedata_raw,mms2_fgm_rdata_raw(:,3),mms1_fgm_btimedata_raw);

[~,mms2_mec_rdata_interp] = interpxyz(mms2_mec_rtimedata_raw,mms2_mec_rdata_raw,mms1_fgm_btimedata_raw);
% mms2_mec_rdata_interp(:,1) = interp1(mms2_mec_rtimedata_raw,mms2_mec_rdata_raw(:,1),mms1_fgm_btimedata_raw);
% mms2_mec_rdata_interp(:,2) = interp1(mms2_mec_rtimedata_raw,mms2_mec_rdata_raw(:,2),mms1_fgm_btimedata_raw);
% mms2_mec_rdata_interp(:,3) = interp1(mms2_mec_rtimedata_raw,mms2_mec_rdata_raw(:,3),mms1_fgm_btimedata_raw);

[~,mms3_fgm_bdata_interp] = interpxyz(mms3_fgm_btimedata_raw,mms3_fgm_bdata_raw,mms1_fgm_btimedata_raw);
% mms3_fgm_bdata_interp(:,1) = interp1(mms3_fgm_btimedata_raw,mms3_fgm_bdata_raw(:,1),mms1_fgm_btimedata_raw);
% mms3_fgm_bdata_interp(:,2) = interp1(mms3_fgm_btimedata_raw,mms3_fgm_bdata_raw(:,2),mms1_fgm_btimedata_raw);
% mms3_fgm_bdata_interp(:,3) = interp1(mms3_fgm_btimedata_raw,mms3_fgm_bdata_raw(:,3),mms1_fgm_btimedata_raw);

% mms3_fgm_rdata_interp(:,1) = interp1(mms3_fgm_rtimedata_raw,mms3_fgm_rdata_raw(:,1),mms1_fgm_btimedata_raw);
% mms3_fgm_rdata_interp(:,2) = interp1(mms3_fgm_rtimedata_raw,mms3_fgm_rdata_raw(:,2),mms1_fgm_btimedata_raw);
% mms3_fgm_rdata_interp(:,3) = interp1(mms3_fgm_rtimedata_raw,mms3_fgm_rdata_raw(:,3),mms1_fgm_btimedata_raw);


[~,mms3_mec_rdata_interp] = interpxyz(mms3_mec_rtimedata_raw,mms3_mec_rdata_raw,mms1_fgm_btimedata_raw);
% mms3_mec_rdata_interp(:,1) = interp1(mms3_mec_rtimedata_raw,mms3_mec_rdata_raw(:,1),mms1_fgm_btimedata_raw);
% mms3_mec_rdata_interp(:,2) = interp1(mms3_mec_rtimedata_raw,mms3_mec_rdata_raw(:,2),mms1_fgm_btimedata_raw);
% mms3_mec_rdata_interp(:,3) = interp1(mms3_mec_rtimedata_raw,mms3_mec_rdata_raw(:,3),mms1_fgm_btimedata_raw);


[~,mms4_fgm_bdata_interp] = interpxyz(mms4_fgm_btimedata_raw,mms4_fgm_bdata_raw,mms1_fgm_btimedata_raw);
% mms4_fgm_bdata_interp(:,1) = interp1(mms4_fgm_btimedata_raw,mms4_fgm_bdata_raw(:,1),mms1_fgm_btimedata_raw);
% mms4_fgm_bdata_interp(:,2) = interp1(mms4_fgm_btimedata_raw,mms4_fgm_bdata_raw(:,2),mms1_fgm_btimedata_raw);
% mms4_fgm_bdata_interp(:,3) = interp1(mms4_fgm_btimedata_raw,mms4_fgm_bdata_raw(:,3),mms1_fgm_btimedata_raw);

% mms4_fgm_rdata_interp(:,1) = interp1(mms4_fgm_rtimedata_raw,mms4_fgm_rdata_raw(:,1),mms1_fgm_btimedata_raw);
% mms4_fgm_rdata_interp(:,2) = interp1(mms4_fgm_rtimedata_raw,mms4_fgm_rdata_raw(:,2),mms1_fgm_btimedata_raw);
% mms4_fgm_rdata_interp(:,3) = interp1(mms4_fgm_rtimedata_raw,mms4_fgm_rdata_raw(:,3),mms1_fgm_btimedata_raw);

[~,mms4_mec_rdata_interp] = interpxyz(mms4_mec_rtimedata_raw,mms4_mec_rdata_raw,mms1_fgm_btimedata_raw);
% mms4_mec_rdata_interp(:,1) = interp1(mms4_mec_rtimedata_raw,mms4_mec_rdata_raw(:,1),mms1_fgm_btimedata_raw);
% mms4_mec_rdata_interp(:,2) = interp1(mms4_mec_rtimedata_raw,mms4_mec_rdata_raw(:,2),mms1_fgm_btimedata_raw);
% mms4_mec_rdata_interp(:,3) = interp1(mms4_mec_rtimedata_raw,mms4_mec_rdata_raw(:,3),mms1_fgm_btimedata_raw);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%Data Cropping and Matrix Formation%%%%%%%%%%%%%%%%%%%%%%%%
%Find the start and end limits of the event in the data, index of, in
%datenum format

%These indices are for plotting the entire event
event_start_index = find(mms1_fgm_btimedata_raw > datenum(event_start,formatIn), 1)-1;
event_end_index = find(mms1_fgm_btimedata_raw > datenum(event_end,formatIn), 1);

%These indices are for the analysis in the specified region
start_index = find(mms1_fgm_btimedata_raw > tstart, 1)-1;
end_index = find(mms1_fgm_btimedata_raw > tend, 1);


%Crop data to specified time range
mms1_fgm_btimedata = mms1_fgm_btimedata_raw(start_index:end_index);
mms2_fgm_btimedata = mms2_fgm_btimedata_raw(start_index:end_index);
mms3_fgm_btimedata = mms3_fgm_btimedata_raw(start_index:end_index);
mms4_fgm_btimedata = mms4_fgm_btimedata_raw(start_index:end_index);
mms_fgm_btimedata = [mms1_fgm_btimedata,mms2_fgm_btimedata,mms3_fgm_btimedata,mms4_fgm_btimedata];

mms1_fgm_bdata_interp = mms1_fgm_bdata_raw(start_index:end_index,1:3);
mms2_fgm_bdata_interp = mms2_fgm_bdata_interp(start_index:end_index,:);
mms3_fgm_bdata_interp = mms3_fgm_bdata_interp(start_index:end_index,:);
mms4_fgm_bdata_interp = mms4_fgm_bdata_interp(start_index:end_index,:);

% mms1_fgm_rdata_interp = mms1_fgm_rdata_interp(start_index:end_index,1:3);
% mms2_fgm_rdata_interp = mms2_fgm_rdata_interp(start_index:end_index,:);
% mms3_fgm_rdata_interp = mms3_fgm_rdata_interp(start_index:end_index,:);
% mms4_fgm_rdata_interp = mms4_fgm_rdata_interp(start_index:end_index,:);

mms1_mec_rdata_interp = mms1_mec_rdata_interp(start_index:end_index,1:3);
mms2_mec_rdata_interp = mms2_mec_rdata_interp(start_index:end_index,:);
mms3_mec_rdata_interp = mms3_mec_rdata_interp(start_index:end_index,:);
mms4_mec_rdata_interp = mms4_mec_rdata_interp(start_index:end_index,:);


%group components of B and r together for faster loops.
bx=[mms1_fgm_bdata_interp(:,1) mms2_fgm_bdata_interp(:,1) mms3_fgm_bdata_interp(:,1) mms4_fgm_bdata_interp(:,1)];
by=[mms1_fgm_bdata_interp(:,2) mms2_fgm_bdata_interp(:,2) mms3_fgm_bdata_interp(:,2) mms4_fgm_bdata_interp(:,2)];
bz=[mms1_fgm_bdata_interp(:,3) mms2_fgm_bdata_interp(:,3) mms3_fgm_bdata_interp(:,3) mms4_fgm_bdata_interp(:,3)];

% rx=[mms1_fgm_rdata_interp(:,1) mms2_fgm_rdata_interp(:,1) mms3_fgm_rdata_interp(:,1) mms4_fgm_rdata_interp(:,1)];
% ry=[mms1_fgm_rdata_interp(:,2) mms2_fgm_rdata_interp(:,2) mms3_fgm_rdata_interp(:,2) mms4_fgm_rdata_interp(:,2)];
% rz=[mms1_fgm_rdata_interp(:,3) mms2_fgm_rdata_interp(:,3) mms3_fgm_rdata_interp(:,3) mms4_fgm_rdata_interp(:,3)];

rx=[mms1_mec_rdata_interp(:,1) mms2_mec_rdata_interp(:,1) mms3_mec_rdata_interp(:,1) mms4_mec_rdata_interp(:,1)];
ry=[mms1_mec_rdata_interp(:,2) mms2_mec_rdata_interp(:,2) mms3_mec_rdata_interp(:,2) mms4_mec_rdata_interp(:,2)];
rz=[mms1_mec_rdata_interp(:,3) mms2_mec_rdata_interp(:,3) mms3_mec_rdata_interp(:,3) mms4_mec_rdata_interp(:,3)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%MDD Calculations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate the Center of the Tetrahedron
rx_0 = mean(rx,2);
ry_0 = mean(ry,2);
rz_0 = mean(rz,2);

%Calculate the relative distance of each probe from the center
rx = rx-rx_0;
ry = ry-ry_0;
rz = rz-rz_0;

%Calculate components Volumetric Tensor
Rxx = (1/4)*sum(rx.^2,2);
Rxy = (1/4)*sum(rx.*ry,2);
Rxz = (1/4)*sum(rx.*rz,2);
Ryy = (1/4)*sum(ry.^2,2);
Ryz = (1/4)*sum(ry.*rz,2);
Rzz = (1/4)*sum(rz.^2,2);

%initialize volumetric tensor matrix
R=zeros(3,3,length(rx));
R_inv = zeros(3,3,length(rx));


%Loop for creating Volumetric Tensor
for i=1:length(rx)
    R(1,1,i) = Rxx(i);
    R(1,2,i) = Rxy(i);
    R(1,3,i) = Rxz(i);
    R(2,1,i) = Rxy(i);
    R(2,2,i) = Ryy(i);
    R(2,3,i) = Ryz(i);
    R(3,1,i) = Rxz(i);
    R(3,2,i) = Ryz(i);
    R(3,3,i) = Rzz(i);
    
    %Inverse Matrix
    R_inv(:,:,i) = inv(R(:,:,i));
end

%For combination of probes
probes = [1,2,3,4];
Summ_probes = nchoosek(probes,2);


% function [grad_B] = gradB(time,Summ_probes,bx,by,bz,rx,ry,rz,R_inv)

%Initialize Matrices and Vectors for G and its eigenvalues/eigenvectors
grad_B = zeros(3,3,length(rx));
eig_lambda = zeros(3,3,length(rx));
eig_vectors = zeros(3,3,length(rx));


%Find the grad_B matrix first, and also its eigenvectors and eigenvalues
for i=1:length(rx)
    %(Grad B)*(Grad_B)T is symmetric
    grad_B(:,:,i) = gradB(i,Summ_probes,probes,bx,by,bz,rx,ry,rz,R_inv);
    L = grad_B(:,:,i)*grad_B(:,:,i)';
    [eig_v, eig_l] = eig(L);
    
    %Order the vectors and eigenvalues from largest to smallest, for all
    %times
    [eig_l,order] = sort(diag(eig_l),'descend');
    eig_v = eig_v(:,order);
    
    %Store each set of eigenvectors and eigenvalues in the 3D array, for
    %each time.
    eig_lambda(:,:,i) = diag(eig_l);
    eig_vectors(:,:,i) = eig_v;
end
%For Calibration Error Elimination, Denton et al. 2010
%Find the largest value for the largest eigenvalue, corresponds to largest G.
Gmax_index = find(squeeze(eig_lambda(1,1,:)) == max(squeeze(eig_lambda(1,1,:))));

%Calculate the components of the average G in the time interval enclosing
%the largest G matrix, or largest eigenvalue
Gxx = mean(squeeze(grad_B(1,1,Gmax_index-avg_i:Gmax_index+avg_i)));
Gyx = mean(squeeze(grad_B(1,2,Gmax_index-avg_i:Gmax_index+avg_i)));
Gzx = mean(squeeze(grad_B(1,3,Gmax_index-avg_i:Gmax_index+avg_i)));

Gxy = mean(squeeze(grad_B(2,1,Gmax_index-avg_i:Gmax_index+avg_i)));
Gyy = mean(squeeze(grad_B(2,2,Gmax_index-avg_i:Gmax_index+avg_i)));
Gzy = mean(squeeze(grad_B(2,3,Gmax_index-avg_i:Gmax_index+avg_i)));

Gxz = mean(squeeze(grad_B(3,1,Gmax_index-avg_i:Gmax_index+avg_i)));
Gyz = mean(squeeze(grad_B(3,2,Gmax_index-avg_i:Gmax_index+avg_i)));
Gzz = mean(squeeze(grad_B(3,3,Gmax_index-avg_i:Gmax_index+avg_i)));

%Average grad_B matrix in the time interval of largest G, or largest eigenvalue
G_0 = abs([Gxx, Gyx, Gzx;...
    Gxy, Gyy, Gzy;...
    Gxz, Gyz, Gzz]);

%Initialize matrices for MDD error calculation
div_B = zeros(length(rx),1);
curl_B = zeros(length(rx),3);


%Calculate the Eigenvalues and Eigenvectors of L = GG', where we have
%subtracted the average value of G in the largest G interval
for i=1:length(rx)
    %Form G*GT and calculate its eigenvalues and eigenvectors
    grad_B(:,:,i) = grad_B(:,:,i) - G_0;
    L = (grad_B(:,:,i))*(grad_B(:,:,i))';
    [eig_v, eig_l] = eig(L);
    
    %Order the eigenvectors and eigenvalues from largest to smallest
    [eig_l,order] = sort(diag(eig_l),'descend');
    eig_v = eig_v(:,order);
    
    %Store the eigenvectors and eigenvalues for each time
    eig_lambda(:,:,i) = diag(eig_l);
    eig_vectors(:,:,i) = eig_v;
    
    %Calculate the MDD error for each time, divB and curlB
    div_B(i,1) = trace(grad_B(:,:,i));
    curl_B(i,1) = grad_B(2,3,i) - grad_B(3,2,i);
    curl_B(i,2) = grad_B(3,1,i) - grad_B(1,3,i);
    curl_B(i,3) = grad_B(1,2,i) - grad_B(2,1,i);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%Sliding Window for MVAB normals%%%%%%%%%%%%%%%%%%%%%
t_0_range = end_index - start_index+1;

%If there is no inputted window size, the sliding window MVAB analysis will default to the same window size
%specified for the MDD analysis
if MVAB_window_size == 0
    MVAB_window_size = end_index - start_index+1;
end

%initialize loop matrices
MVAB_lambda_all=zeros(3,t_0_range);
MVAB_normal_all=zeros(3,t_0_range);
MVAB_time=zeros(1,t_0_range);

%Sliding Window MVAB for every point same as MDD analysis
for i=1:t_0_range
    t_0 = mms1_fgm_btimedata_raw(i+start_index-1); %center point, loop through all of these between the start/end indices
    left_index = i+start_index-MVAB_window_size/2-1; %determine the left boundary from the center point, half of window size.
    right_index = i+start_index+MVAB_window_size/2-1; %determine the right boundary from the center point, half of window size.
    
    bdata_scope = mms1_fgm_bdata_raw(left_index:right_index,:); %Grab the Bdata and splice it to just within the left and right boundaries.
    
    [output, l, v] = mvab(bdata_scope(:,1:3)); %Calculate the MVAB from the spliced time interval
    
    MVAB_lambda_all(:,i) = l; %Store the eigenvalues
    MVAB_normal_all(:,i) = v(:,mvamm);%Store only eigenvector chosen at the beginning of the script
    
    MVAB_time(1,i)=t_0;%Store the center point(time) used for MVAB analysis
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Timing Method%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_diff = zeros(6,3);
t_diff = zeros(6,1);
midpoint_index = floor(t_0_range/2);
%creating the r_diff matrix
for i=1:6
r_diff(i,1)=rx(midpoint_index,Summ_probes(i,1))-rx(midpoint_index,Summ_probes(i,2));
r_diff(i,2)=ry(midpoint_index,Summ_probes(i,1))-ry(midpoint_index,Summ_probes(i,2));
r_diff(i,3)=rz(midpoint_index,Summ_probes(i,1))-rz(midpoint_index,Summ_probes(i,2));
t_diff(i) = mms_fgm_btimedata(midpoint_index,Summ_probes(i,1))-mms_fgm_btimedata(midpoint_index,Summ_probes(i,2));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Unfluctuate%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%If needed and specified in the beginning, flip the equal and opposite
%signs of the eigenvectors. Calculate the mean, and determine the overal
%sign of the eigenvectors

if fluctuate_mdd == 1
    mdd_nx = Unfluctuate(squeeze(eig_vectors(1,mddmm,:)));
    mdd_ny = Unfluctuate(squeeze(eig_vectors(2,mddmm,:)));
    mdd_nz = Unfluctuate(squeeze(eig_vectors(3,mddmm,:)));
elseif fluctuate_mdd == 0
    mdd_nx = (squeeze(eig_vectors(1,mddmm,:)));
    mdd_ny = (squeeze(eig_vectors(2,mddmm,:)));
    mdd_nz = (squeeze(eig_vectors(3,mddmm,:)));
    
end

if fluctuate_mva == 1
    MVAB_normal_all(1,:) = Unfluctuate(MVAB_normal_all(1,:));
    MVAB_normal_all(2,:) = Unfluctuate(MVAB_normal_all(2,:));
    MVAB_normal_all(3,:) = Unfluctuate(MVAB_normal_all(3,:));
elseif fluctuate_mva == 0
    MVAB_normal_all(1,:) = (MVAB_normal_all(1,:));
    MVAB_normal_all(2,:) = (MVAB_normal_all(2,:));
    MVAB_normal_all(3,:) = (MVAB_normal_all(3,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%Average Angle Calculations%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%Mean MDD and Center MVA Angle Calculation %%%%%%%%%%%%%%%%


%Calulation for the mean G_matrix of all time.
AA = mean(squeeze(grad_B(1,1,:)));
BB = mean(squeeze(grad_B(1,2,:)));
CC = mean(squeeze(grad_B(1,3,:)));
DD = mean(squeeze(grad_B(2,1,:)));
EE = mean(squeeze(grad_B(2,2,:)));
FF = mean(squeeze(grad_B(2,3,:)));
GG = mean(squeeze(grad_B(3,1,:)));
HH = mean(squeeze(grad_B(3,2,:)));
II = mean(squeeze(grad_B(3,3,:)));

%Form the mean G_matrix of all time
KK = [AA, BB, CC; DD, EE, FF; GG, HH, II];
%Calculate the eigenvalues and eigenvectors for the mean GGT
[VV,DD] = eig(KK*KK');

%Order the eigenvectors and eigenvalues from largest to smallest
[DD,order] = sort(diag(DD),'descend');
VV = VV(:,order);

%Store the specified eigenvector
Mean_grad_B_MDD_normal = VV(:,mddmm);


%if initial conditions did not have provide for a MVAB normal 
%over the specified time interval, find it.
if MVA_normal == [0, 0, 0] %#ok<BDSCA,BDSCI>
    %Calculate the MVAB for the specified time interval
    [~,mvab_l,mvab_v] = mvab(mms1_fgm_bdata_interp);
    
    %Store the eigenvalues and eigenvectors
    MVA_normal = ([mvab_v(1,mvamm), mvab_v(2,mvamm), mvab_v(3,mvamm)]);
    MVA_values = (mvab_l) %#ok<NOPTS>
    
    MVA_min_normal = ([mvab_v(1,3), mvab_v(2,3), mvab_v(3,3)]) %Store the Minimum Eigenvector direction, specifically.
    
end

%Because MVA is ambiguous for one direction, determine the best normal
%direction to compare with MDD. The best normal direction is the one with
%the smallest angle with MDD.

%Calculate all 4 possibilities
MVA_ambig(1)=angle(MVA_normal.*[-1,1,1],Mean_grad_B_MDD_normal);
MVA_ambig(2)=angle(MVA_normal.*[1,-1,1],Mean_grad_B_MDD_normal);
MVA_ambig(3)=angle(MVA_normal.*[1,1,-1],Mean_grad_B_MDD_normal);
MVA_ambig(4)=angle(MVA_normal.*[1,1,1],Mean_grad_B_MDD_normal);

%Find the normal with the smallest angle with MDD, and replace the edited normal with the old
%normal.
switch find(MVA_ambig==min(MVA_ambig))
    case 1
        MVA_normal = MVA_normal.*[-1,1,1];
    case 2
        MVA_normal = MVA_normal.*[1,-1,1];
    case 3
        MVA_normal = MVA_normal.*[1,1,-1];
    case 4
        MVA_normal = MVA_normal; %#ok<ASGSL>
end

%The smallest angle between MDD and MVA, possible
mean_MDD_dot_MVA_center = min(MVA_ambig);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%Calculation of Angle at every point%%%%%%%%%%%%%%%%%%%%
%Calculations for Point by Point Analysis of MDD and MVA
%Calculate average of each component for the MDDB Normal
mean_x = mean(mdd_nx,'omitnan');
mean_y = mean(mdd_ny,'omitnan');
mean_z = mean(mdd_nz,'omitnan');
% mean_x = mean(mdd_nx,'omitnan');
% mean_y = mean(mdd_ny,'omitnan');
% mean_z = mean(mdd_nz,'omitnan');
MDD_mean_normal = ([mean_x, mean_y, mean_z])./norm([mean_x, mean_y, mean_z]);

%Set up the point by points of MDD and MVA
MDDB_normal_all = [mdd_nx, mdd_ny, mdd_nz]; %MDDB normal at every point
MVAB_normal_all = MVAB_normal_all'; %MVAB normals at every point

%Initialize the matrices
MVA_dot_MDD_angle_all = zeros(length(rx),1);
MVA_dot_MDD_center_all = zeros(length(rx),1);
MVA_ambig_dot_center= zeros(length(rx),4);

for i=1:length(rx)
    %fix the 180 degrees ambiguity by choosing the smallest angle between
    %MDD
    MVA_ambig_all(1)=angle(MVAB_normal_all(i,:).*[-1,1,1],MDDB_normal_all(i,:));
    MVA_ambig_all(2)=angle(MVAB_normal_all(i,:).*[1,-1,1],MDDB_normal_all(i,:));
    MVA_ambig_all(3)=angle(MVAB_normal_all(i,:).*[1,1,-1],MDDB_normal_all(i,:));
    MVA_ambig_all(4)=angle(MVAB_normal_all(i,:),MDDB_normal_all(i,:));
    
    MVA_dot_MDD_angle_all(i,1) = min(MVA_ambig_all); %Store the minimum angle
    
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%Custom Plotting
% plots = [B_all';'B_portion';'MDD_Normal';'MDD_Lambda';'MVA_Normal';'MVA_Lambda';'Angle'];
% 
% for i=1:length(plots)
%     switch plot(i)
%         case 'B_all'
%         case 'B_portion'
%         case 'MDD_Normal'
%         case 'MDD_Lambda'
%         case 'MVA_Normal'
%         case 'MVA_Lambda'
%         case 'Angle'
%     end
%     
%     
% end


%%
%%%%%%%%%%%%%%%%%%%%%%%Plot Magnetic Field Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1 = subplot(num_plots,1,1);
set(gcf,'color','w');
co = [0 0 1;
    0 1 0;
    1 0 0;
    0 0 0];
set(gcf,'defaultAxesColorOrder',co)

plot(mms1_fgm_btimedata_raw(event_start_index:event_end_index,1),mms1_fgm_bdata_raw(event_start_index:event_end_index,1:4))

line([(mms1_fgm_btimedata_raw(start_index)) (mms1_fgm_btimedata_raw(start_index))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
line([(mms1_fgm_btimedata_raw(end_index)) (mms1_fgm_btimedata_raw(end_index))],...
    get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')

box on
datetick
set(gca,'layer','top')
legend({'B_x', 'B_y', 'B_z','B_{tot}'},'FontSize',10)
legend('boxoff');
legend('Location','eastoutside');
ylabel({'B';'[nT]'},'FontSize', 14)
xlim([mms1_fgm_btimedata_raw(event_start_index), mms1_fgm_btimedata_raw(event_end_index)])
set(f1,'XMinorTick','on','YMinorTick','on')

title_name = 'MMS MDDB&MVAB Normals';
title(title_name, 'FontSize', 18, 'FontWeight', 'normal')


plot_pos = get(f1,'Position');
set(f1,'Position',[plot_pos(1), plot_pos(2), plot_pos(3), plot_pos(4)]);

l1= findobj(gcf, 'Type', 'Legend');
l1pos=get(l1, 'Position');
set(l1,'Position',[plot_pos(1)+plot_pos(3)-2, plot_pos(2)+plot_pos(4)/2-l1pos(4)/2, l1pos(3), l1pos(4)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%Plot Portion of B%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f2 = subplot(num_plots,1,2);

plot(mms1_fgm_btimedata_raw(start_index:end_index),mms1_fgm_bdata_raw(start_index:end_index,1:3))
box on

set(gca,'layer','top')
legend({'B_x', 'B_y', 'B_z'},'FontSize',10)
legend('boxoff');
legend('Location','eastoutside');
datetick('x','SS.FFF','keepticks')

ylabel({'B';'[nT]'},'FontSize', 14)
datetick
xlim([mms1_fgm_btimedata(1) mms1_fgm_btimedata(end)])

set(f2,'XMinorTick','on','YMinorTick','on','XTickLabel',[])

set(f2,'Position',[plot_pos(1), plot_pos(2)-1*plot_pos(4)-0.05, plot_pos(3), plot_pos(4)]);

l1= findobj(gcf, 'Type', 'Legend');
l1pos=get(l1(2), 'Position');
set(l1(2),'Position',[plot_pos(1)+plot_pos(3)+0.010, plot_pos(2)+1*plot_pos(4)-0.1, l1pos(3), l1pos(4)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot MDD Normals%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f3 = subplot(num_plots,1,3);

%Minimum
plot(mms1_fgm_btimedata,mdd_nx)
hold on
plot(mms1_fgm_btimedata,mdd_ny)
plot(mms1_fgm_btimedata,mdd_nz)
hold off

box on
datetick
set(gca,'layer','top')
legend({'N_x', 'N_y', 'N_z'},'FontSize',10)
legend('boxoff');
legend('Location','eastoutside');

ylabel(strcat('MDD_',mddl),'FontSize', 14)


xlim([mms1_fgm_btimedata(1) mms1_fgm_btimedata(end)])

set(f3,'XMinorTick','on','YMinorTick','on','XTickLabel',[])

set(f3,'Position',[plot_pos(1), plot_pos(2)-2*plot_pos(4)-0.075, plot_pos(3), plot_pos(4)*1]);

l1= findobj(gcf, 'Type', 'Legend');
l1pos=get(l1(3), 'Position');
%set(l1(3),'Position',[plot_pos(1)+plot_pos(3)+0.075, plot_pos(2)+2*plot_pos(4)+0.15, l1pos(3), l1pos(4)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot MDD Eigenvalues%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f4 = subplot(num_plots,1,4);

plot(mms1_fgm_btimedata,sqrt(squeeze(eig_lambda(1,1,:))./squeeze(eig_lambda(2,2,:))))
hold on
plot(mms1_fgm_btimedata,sqrt(squeeze(eig_lambda(2,2,:))/squeeze(eig_lambda(3,3,:))),'g')
%plot(mms1_fgm_btimedata,sqrt(squeeze(eig_lambda(3,3,:))/squeeze(eig_lambda(1,1,:))),'r')
hold off

datetick
set(gca,'layer','top')
legend({'\lambda_{max}/\lambda_{mid}', '\lambda_{mid}/\lambda_{min}', '\lambda_{min}/\lambda_{max}'},'FontSize',10)
legend('boxoff');
legend('Location','eastoutside');

ylabel(strcat('\lambda_{MDD}'),'FontSize', 14)

xlim([mms1_fgm_btimedata(1) mms1_fgm_btimedata(end)])

set(f4,'XMinorTick','on','YMinorTick','on','XTickLabel',[])

set(f4,'Position',[plot_pos(1), plot_pos(2)-3*plot_pos(4)-0.1, plot_pos(3), plot_pos(4)*1]);

% l1= findobj(gcf, 'Type', 'Legend');
% l1pos=get(l1(4), 'Position');
% set(l1(4),'Position',[plot_pos(1)+plot_pos(3)+0.1, plot_pos(2)+3*plot_pos(4), l1pos(3), l1pos(4)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot MDD Error%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f5 = subplot(num_plots,1,5);

plot(mms1_fgm_btimedata,abs(div_B./(sqrt(curl_B(:,1).^2+curl_B(:,2).^2+curl_B(:,3).^2))))
hold on
plot(mms1_fgm_btimedata,abs(div_B./(norm(G_0))))
line(xlim(), [0.4,0.4], 'LineWidth', 0.25, 'Color', 'k','LineStyle','--');
hold off

datetick
set(gca,'layer','top')
legend({'\nabla\cdotB/\nablaxB', '\nabla\cdotB/max(\nablaB)'},'FontSize',10)
legend('boxoff');
legend('Location','eastoutside');
%xlim([datetime(tstart,'ConvertFrom','datenum') datetime(tend,'ConvertFrom','datenum')]);

ylabel(strcat('MDD_{Error}'),'FontSize', 14)

xlim([mms1_fgm_btimedata(1) mms1_fgm_btimedata(end)])

set(f5,'XMinorTick','on','YMinorTick','on','XTickLabel',[])

set(f5,'Position',[plot_pos(1), plot_pos(2)-4*plot_pos(4)-0.125, plot_pos(3), plot_pos(4)*1]);

% l1= findobj(gcf, 'Type', 'Legend');
% l1pos=get(l1(4), 'Position');
% set(l1(4),'Position',[plot_pos(1)+plot_pos(3)+0.1, plot_pos(2)+3*plot_pos(4), l1pos(3), l1pos(4)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot MVAB Normals%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f6 = subplot(num_plots,1,6);

plot(mms1_fgm_btimedata,MVAB_normal_all')

datetick
set(gca,'layer','top')
legend({'N_x', 'N_y', 'N_z'},'FontSize',10)
legend('boxoff');
legend('Location','eastoutside');

ylabel(strcat('MVAB_',mval),'FontSize', 14)

xlim([mms1_fgm_btimedata(1) mms1_fgm_btimedata(end)])

set(f6,'XMinorTick','on','YMinorTick','on','XTickLabel',[])

set(f6,'Position',[plot_pos(1), plot_pos(2)-5*plot_pos(4)-0.15, plot_pos(3), plot_pos(4)*1]);

% l1= findobj(gcf, 'Type', 'Legend');
% l1pos=get(l1(4), 'Position');
% set(l1(4),'Position',[plot_pos(1)+plot_pos(3)+0.1, plot_pos(2)+3*plot_pos(4), l1pos(3), l1pos(4)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot MVAB Lambdas%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f7 = subplot(num_plots,1,7);
MVAB_lambda_12 = MVAB_lambda_all(1,:)./MVAB_lambda_all(2,:);
% MVAB_lambda_12(MVAB_lambda_12 > 25) = NaN;

MVAB_lambda_23 = MVAB_lambda_all(2,:)./MVAB_lambda_all(3,:);
% MVAB_lambda_23(MVAB_lambda_23 > 25) = NaN;

plot(mms1_fgm_btimedata,MVAB_lambda_12)
hold on
plot(mms1_fgm_btimedata,MVAB_lambda_23)
hold off

datetick
set(gca,'layer','top')
legend({'\lambda_{max}/\lambda_{min}','\lambda_{mid}/\lambda_{min}'},'FontSize',10)
legend('boxoff');
legend('Location','eastoutside');
%xlim([datetime(tstart,'ConvertFrom','datenum') datetime(tend,'ConvertFrom','datenum')]);

ylabel(strcat('\lambda_{MVAB}'),'FontSize', 14)
% ylim([0 25])
xlim([mms1_fgm_btimedata(1) mms1_fgm_btimedata(end)])

set(f7,'XMinorTick','on','YMinorTick','on','XTickLabel',[])

set(f7,'Position',[plot_pos(1), plot_pos(2)-6*plot_pos(4)-0.175, plot_pos(3), plot_pos(4)*1]);

% l1= findobj(gcf, 'Type', 'Legend');
% l1pos=get(l1(4), 'Position');
% set(l1(4),'Position',[plot_pos(1)+plot_pos(3)+0.1, plot_pos(2)+3*plot_pos(4), l1pos(3), l1pos(4)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot Angles%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f8 = subplot(num_plots,1,8);

plot(mms1_fgm_btimedata,MVA_dot_MDD_angle_all)
box on
datetick('x','SS.FFF','keeplimits')
set(gca,'layer','top')
legend({'Sliding MVA'},'FontSize',10)
legend('boxoff');
legend('Location','eastoutside');

ylabel({'Angle'; '[Deg]'},'FontSize', 14)

xlim([mms1_fgm_btimedata(1) mms1_fgm_btimedata(end)])
set(f8,'XMinorTick','on','YMinorTick','on')

set(f8,'Position',[plot_pos(1), plot_pos(2)-7*plot_pos(4)-0.2, plot_pos(3), plot_pos(4)*1]);


% l1= findobj(gcf, 'Type', 'Legend');
% l1pos=get(l1(5), 'Position');
% set(l1(8),'Position',[plot_pos(1)+plot_pos(3)+0.1, plot_pos(2)+4*plot_pos(4), l1pos(3), l1pos(4)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%Text on Plot%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% annotation('textbox',[plot_pos(1), plot_pos(2)-8*plot_pos(4)-.18, plot_pos(3), plot_pos(4)*1],...
%     'String',{date_start(1:10)},...
%     'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);


annotation('textbox',[plot_pos(1), plot_pos(2)-8*plot_pos(4)-.18, plot_pos(3), plot_pos(4)*1],...
    'String',strcat({date_start},' - ',{date_end(11:end)}),...
    'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);

% annotation('textbox',[plot_pos(1), plot_pos(2)-8*plot_pos(4)-.20, plot_pos(3), plot_pos(4)*1],...
%     'String',{date_end},...
%     'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);

% annotation('textbox',[plot_pos(1), plot_pos(2)-8*plot_pos(4)-0.180, plot_pos(3), plot_pos(4)*1],...
%     'String',strcat('(Median MDD)\cdot(Centered MVA) Angle= ', ' ', num2str(median(MVA_dot_MDD_angle_all),'%2.2f'),char(176)),...
%     'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);

annotation('textbox',[plot_pos(1), plot_pos(2)-8*plot_pos(4)-0.20, plot_pos(3), plot_pos(4)*1],...
    'String',strcat('(Average MDD)\cdot(Centered MVA) Angle= ', ' ', num2str(mean_MDD_dot_MVA_center,'%2.2f'),char(176)),...
    'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);


% annotation('textbox',[plot_pos(1), plot_pos(2)-8*plot_pos(4)-0.20, plot_pos(3), plot_pos(4)*1],...
%     'String',strcat('(Mean MDD)\cdot(Centered MVA) Angle= ', ' ', num2str(MDD_dot_MVA_center,'%2.2f'),char(176)),...
%     'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);

annotation('textbox',[plot_pos(1)+0.55*plot_pos(3), plot_pos(2)-8*plot_pos(4)-0.180, plot_pos(3), plot_pos(4)*1],...
    'String',strcat('Average MDD Normal = ', ' ', num2str(Mean_grad_B_MDD_normal',' [%0.4f, %0.4f, %0.4f]')),...
    'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);

% annotation('textbox',[plot_pos(1)+0.55*plot_pos(3), plot_pos(2)-8*plot_pos(4)-0.20, plot_pos(3), plot_pos(4)*1],...
%     'String',strcat('Mean MDD Normal = ', ' ', num2str(MDD_mean_normal,' [%0.4f, %0.4f, %0.4f]')),...
%     'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);

annotation('textbox',[plot_pos(1)+0.55*plot_pos(3), plot_pos(2)-8*plot_pos(4)-0.20, plot_pos(3), plot_pos(4)*1],...
    'String',strcat('Center MVA Normal = ', ' ', num2str(MVA_min_normal,' [%0.4f, %0.4f, %0.4f]')),...
    'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


orient(gcf,'landscape')

plot_name =  strcat('mms','_mdd_',mddl,'_mva_',mval,'_angles2',...
    date_start(1:23),'_avg_i',num2str(avg_i),'_fluctmdd',num2str(fluctuate_mdd),'_fluctmva',num2str(fluctuate_mva),'_MVAB_window',num2str(MVAB_window_size),'.pdf');

print(gcf, '-dpdf', '-opengl',plot_name,'-fillpage');

movefile(plot_name, '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis')



toc


