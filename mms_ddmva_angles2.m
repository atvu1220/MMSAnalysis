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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Data Retrieval%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

event_start = '2018-03-12 07:35:50.000';
event_end = '2018-03-12 07:36:30.000';

%Input Time Range directly
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


event_start = '2018-03-01 01:05:40.000';
event_end = '2018-03-01 01:06:22.000';

date_start = '2018-03-01 01:06:10.500';
date_end = '2018-03-01 01:06:11.250';



% event_start = '2018-01-09 08:34:27.000';
% event_end = '2018-01-09 08:35:03.000';
% 
% date_start = '2018-01-09 08:34:49.000';
% date_end = '2018-01-09 8:34:50.000';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Adjustable Parameters%%%%%%%%%%%%%%%%%%%%%%%%%
avg_i = 6;
fluctuate_mdd = 0;%1 for force all vectors to stay within their same sign
fluctuate_mva = 0;
MVAB_window_size = 1024; %this is the number of data points on either side of t_0; enter 0 for default window size of the time interval

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
% [mms1_mec_rtimedata_raw, mms1_mec_rdata_raw] = load_mec(date_start,1,'srvy');
% [mms2_mec_rtimedata_raw, mms2_mec_rdata_raw] = load_mec(date_start,2,'srvy');
% [mms3_mec_rtimedata_raw, mms3_mec_rdata_raw] = load_mec(date_start,3,'srvy');
% [mms4_mec_rtimedata_raw, mms4_mec_rdata_raw] = load_mec(date_start,4,'srvy');



%%%%%%%%%%%%%%%%%%%%Interpolate position and bdata with bdata from mms1 %%%%%%%%%%%%%%%%%%%

%Interpolate the Data
mms1_fgm_rdata_interp(:,1) = interp1(mms1_fgm_rtimedata_raw,mms1_fgm_rdata_raw(:,1),mms1_fgm_btimedata_raw);
mms1_fgm_rdata_interp(:,2) = interp1(mms1_fgm_rtimedata_raw,mms1_fgm_rdata_raw(:,2),mms1_fgm_btimedata_raw);
mms1_fgm_rdata_interp(:,3) = interp1(mms1_fgm_rtimedata_raw,mms1_fgm_rdata_raw(:,3),mms1_fgm_btimedata_raw);

% mms1_mec_rdata_interp(:,1) = interp1(mms1_mec_rtimedata_raw,mms1_mec_rdata_raw(:,1),mms1_fgm_btimedata_raw);
% mms1_mec_rdata_interp(:,2) = interp1(mms1_mec_rtimedata_raw,mms1_mec_rdata_raw(:,2),mms1_fgm_btimedata_raw);
% mms1_mec_rdata_interp(:,3) = interp1(mms1_mec_rtimedata_raw,mms1_mec_rdata_raw(:,3),mms1_fgm_btimedata_raw);

mms2_fgm_bdata_interp(:,1) = interp1(mms2_fgm_btimedata_raw,mms2_fgm_bdata_raw(:,1),mms1_fgm_btimedata_raw);
mms2_fgm_bdata_interp(:,2) = interp1(mms2_fgm_btimedata_raw,mms2_fgm_bdata_raw(:,2),mms1_fgm_btimedata_raw);
mms2_fgm_bdata_interp(:,3) = interp1(mms2_fgm_btimedata_raw,mms2_fgm_bdata_raw(:,3),mms1_fgm_btimedata_raw);

mms2_fgm_rdata_interp(:,1) = interp1(mms2_fgm_rtimedata_raw,mms2_fgm_rdata_raw(:,1),mms1_fgm_btimedata_raw);
mms2_fgm_rdata_interp(:,2) = interp1(mms2_fgm_rtimedata_raw,mms2_fgm_rdata_raw(:,2),mms1_fgm_btimedata_raw);
mms2_fgm_rdata_interp(:,3) = interp1(mms2_fgm_rtimedata_raw,mms2_fgm_rdata_raw(:,3),mms1_fgm_btimedata_raw);

% mms2_mec_rdata_interp(:,1) = interp1(mms2_mec_rtimedata_raw,mms2_mec_rdata_raw(:,1),mms1_fgm_btimedata_raw);
% mms2_mec_rdata_interp(:,2) = interp1(mms2_mec_rtimedata_raw,mms2_mec_rdata_raw(:,2),mms1_fgm_btimedata_raw);
% mms2_mec_rdata_interp(:,3) = interp1(mms2_mec_rtimedata_raw,mms2_mec_rdata_raw(:,3),mms1_fgm_btimedata_raw);

mms3_fgm_bdata_interp(:,1) = interp1(mms3_fgm_btimedata_raw,mms3_fgm_bdata_raw(:,1),mms1_fgm_btimedata_raw);
mms3_fgm_bdata_interp(:,2) = interp1(mms3_fgm_btimedata_raw,mms3_fgm_bdata_raw(:,2),mms1_fgm_btimedata_raw);
mms3_fgm_bdata_interp(:,3) = interp1(mms3_fgm_btimedata_raw,mms3_fgm_bdata_raw(:,3),mms1_fgm_btimedata_raw);

mms3_fgm_rdata_interp(:,1) = interp1(mms3_fgm_rtimedata_raw,mms3_fgm_rdata_raw(:,1),mms1_fgm_btimedata_raw);
mms3_fgm_rdata_interp(:,2) = interp1(mms3_fgm_rtimedata_raw,mms3_fgm_rdata_raw(:,2),mms1_fgm_btimedata_raw);
mms3_fgm_rdata_interp(:,3) = interp1(mms3_fgm_rtimedata_raw,mms3_fgm_rdata_raw(:,3),mms1_fgm_btimedata_raw);

% mms3_mec_rdata_interp(:,1) = interp1(mms3_mec_rtimedata_raw,mms3_mec_rdata_raw(:,1),mms1_fgm_btimedata_raw);
% mms3_mec_rdata_interp(:,2) = interp1(mms3_mec_rtimedata_raw,mms3_mec_rdata_raw(:,2),mms1_fgm_btimedata_raw);
% mms3_mec_rdata_interp(:,3) = interp1(mms3_mec_rtimedata_raw,mms3_mec_rdata_raw(:,3),mms1_fgm_btimedata_raw);

mms4_fgm_bdata_interp(:,1) = interp1(mms4_fgm_btimedata_raw,mms4_fgm_bdata_raw(:,1),mms1_fgm_btimedata_raw);
mms4_fgm_bdata_interp(:,2) = interp1(mms4_fgm_btimedata_raw,mms4_fgm_bdata_raw(:,2),mms1_fgm_btimedata_raw);
mms4_fgm_bdata_interp(:,3) = interp1(mms4_fgm_btimedata_raw,mms4_fgm_bdata_raw(:,3),mms1_fgm_btimedata_raw);

mms4_fgm_rdata_interp(:,1) = interp1(mms4_fgm_rtimedata_raw,mms4_fgm_rdata_raw(:,1),mms1_fgm_btimedata_raw);
mms4_fgm_rdata_interp(:,2) = interp1(mms4_fgm_rtimedata_raw,mms4_fgm_rdata_raw(:,2),mms1_fgm_btimedata_raw);
mms4_fgm_rdata_interp(:,3) = interp1(mms4_fgm_rtimedata_raw,mms4_fgm_rdata_raw(:,3),mms1_fgm_btimedata_raw);

% mms4_mec_rdata_interp(:,1) = interp1(mms4_mec_rtimedata_raw,mms4_mec_rdata_raw(:,1),mms1_fgm_btimedata_raw);
% mms4_mec_rdata_interp(:,2) = interp1(mms4_mec_rtimedata_raw,mms4_mec_rdata_raw(:,2),mms1_fgm_btimedata_raw);
% mms4_mec_rdata_interp(:,3) = interp1(mms4_mec_rtimedata_raw,mms4_mec_rdata_raw(:,3),mms1_fgm_btimedata_raw);


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

mms1_fgm_bdata_interp = mms1_fgm_bdata_raw(start_index:end_index,1:3);
mms2_fgm_bdata_interp = mms2_fgm_bdata_interp(start_index:end_index,:);
mms3_fgm_bdata_interp = mms3_fgm_bdata_interp(start_index:end_index,:);
mms4_fgm_bdata_interp = mms4_fgm_bdata_interp(start_index:end_index,:);

mms1_fgm_rdata_interp = mms1_fgm_rdata_interp(start_index:end_index,1:3);
mms2_fgm_rdata_interp = mms2_fgm_rdata_interp(start_index:end_index,:);
mms3_fgm_rdata_interp = mms3_fgm_rdata_interp(start_index:end_index,:);
mms4_fgm_rdata_interp = mms4_fgm_rdata_interp(start_index:end_index,:);

% mms1_mec_rdata_interp = mms1_mec_rdata_interp(start_index:end_index,1:3);
% mms2_mec_rdata_interp = mms2_mec_rdata_interp(start_index:end_index,:);
% mms3_mec_rdata_interp = mms3_mec_rdata_interp(start_index:end_index,:);
% mms4_mec_rdata_interp = mms4_mec_rdata_interp(start_index:end_index,:);




%group components of B and r together for faster loops.
bx=[mms1_fgm_bdata_interp(:,1) mms2_fgm_bdata_interp(:,1) mms3_fgm_bdata_interp(:,1) mms4_fgm_bdata_interp(:,1)];
by=[mms1_fgm_bdata_interp(:,2) mms2_fgm_bdata_interp(:,2) mms3_fgm_bdata_interp(:,2) mms4_fgm_bdata_interp(:,2)];
bz=[mms1_fgm_bdata_interp(:,3) mms2_fgm_bdata_interp(:,3) mms3_fgm_bdata_interp(:,3) mms4_fgm_bdata_interp(:,3)];

rx=[mms1_fgm_rdata_interp(:,1) mms2_fgm_rdata_interp(:,1) mms3_fgm_rdata_interp(:,1) mms4_fgm_rdata_interp(:,1)];
ry=[mms1_fgm_rdata_interp(:,2) mms2_fgm_rdata_interp(:,2) mms3_fgm_rdata_interp(:,2) mms4_fgm_rdata_interp(:,2)];
rz=[mms1_fgm_rdata_interp(:,3) mms2_fgm_rdata_interp(:,3) mms3_fgm_rdata_interp(:,3) mms4_fgm_rdata_interp(:,3)];

% rx=[mms1_mec_rdata_interp(:,1) mms2_mec_rdata_interp(:,1) mms3_mec_rdata_interp(:,1) mms4_mec_rdata_interp(:,1)];
% ry=[mms1_mec_rdata_interp(:,2) mms2_mec_rdata_interp(:,2) mms3_mec_rdata_interp(:,2) mms4_mec_rdata_interp(:,2)];
% rz=[mms1_mec_rdata_interp(:,3) mms2_mec_rdata_interp(:,3) mms3_mec_rdata_interp(:,3) mms4_mec_rdata_interp(:,3)];

%Calculate the Center of the Tetrahedron
rx_0 = mean(rx,2);
ry_0 = mean(ry,2);
rz_0 = mean(rz,2);

%Calculate the relative distance of each probe from the center
rx = rx-rx_0;
ry = ry-ry_0;
rz = rz-rz_0;

%Form Volumetric Tensor
Rxx = (1/4)*sum(rx.^2,2);
Rxy = (1/4)*sum(rx.*ry,2);
Rxz = (1/4)*sum(rx.*rz,2);
Ryy = (1/4)*sum(ry.^2,2);
Ryz = (1/4)*sum(ry.*rz,2);
Rzz = (1/4)*sum(rz.^2,2);

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


probes = [1,2,3,4];
Summ_probes = nchoosek(probes,2);


% function [grad_B] = gradB(time,Summ_probes,bx,by,bz,rx,ry,rz,R_inv)
grad_B = zeros(3,3,length(rx));
eig_lambda = zeros(3,3,length(rx));
eig_vectors = zeros(3,3,length(rx));


%Find the grad_B matrix first, and also its eigenvectors and eigenvalues
for i=1:length(rx)
    grad_B(:,:,i) = gradB(i,Summ_probes,probes,bx,by,bz,rx,ry,rz,R_inv);
    L = grad_B(:,:,i)*grad_B(:,:,i)';
    [eig_v, eig_l] = eig(L);
    
    
    [eig_l,order] = sort(diag(eig_l),'descend');
    eig_v = eig_v(:,order);
    
    eig_lambda(:,:,i) = diag(eig_l);
    eig_vectors(:,:,i) = eig_v;
end

%constract multidimensional array for the Gradient of B at every time
% for i=1:length(rx)
%     grad_B(:,:,i) = gradB(i,Summ_probes,bx,by,bz,rx,ry,rz,R_inv);
% end


%Find the largest value for the largest eigenvalue, corresponds to largest G.
Gmax_index = find(squeeze(eig_lambda(1,1,:)) == max(squeeze(eig_lambda(1,1,:))));


%Calculate the average of G in the time interval, Denton et al. 2010
%average_interval

Gxx = mean(squeeze(grad_B(1,1,Gmax_index-avg_i:Gmax_index+avg_i)));
Gyx = mean(squeeze(grad_B(1,2,Gmax_index-avg_i:Gmax_index+avg_i)));
Gzx = mean(squeeze(grad_B(1,3,Gmax_index-avg_i:Gmax_index+avg_i)));

Gxy = mean(squeeze(grad_B(2,1,Gmax_index-avg_i:Gmax_index+avg_i)));
Gyy = mean(squeeze(grad_B(2,2,Gmax_index-avg_i:Gmax_index+avg_i)));
Gzy = mean(squeeze(grad_B(2,3,Gmax_index-avg_i:Gmax_index+avg_i)));

Gxz = mean(squeeze(grad_B(3,1,Gmax_index-avg_i:Gmax_index+avg_i)));
Gyz = mean(squeeze(grad_B(3,2,Gmax_index-avg_i:Gmax_index+avg_i)));
Gzz = mean(squeeze(grad_B(3,3,Gmax_index-avg_i:Gmax_index+avg_i)));

%average grad_B in the time interval of largest G, or largest eigenvalue
G_0 = abs([Gxx, Gyx, Gzx;...
    Gxy, Gyy, Gzy;...
    Gxz, Gyz, Gzz]);
div_B = zeros(length(rx),1);
curl_B = zeros(length(rx),3);


%Calculate the Eigenvalues and Eigenvectors of L = GG', wherewe have
%subtracted the average value of G in the largest G interval
for i=1:length(rx)
    grad_B(:,:,i) = grad_B(:,:,i) - G_0;
    L = (grad_B(:,:,i))*(grad_B(:,:,i))';
    [eig_v, eig_l] = eig(L);
    
    
    [eig_l,order] = sort(diag(eig_l),'descend');
    eig_v = eig_v(:,order);
    
    eig_lambda(:,:,i) = diag(eig_l);
    eig_vectors(:,:,i) = eig_v;
    
    div_B(i,1) = trace(grad_B(:,:,i));
    curl_B(i,1) = grad_B(2,3,i) - grad_B(3,2,i);
    curl_B(i,2) = grad_B(3,1,i) - grad_B(1,3,i);
    curl_B(i,3) = grad_B(1,2,i) - grad_B(2,1,i);
end


% %Filtering data
% %%no longer needed as we have sorted the eigenvalues/eigenvectors
% acceptable_lambda_ratios = squeeze(eig_lambda(1,1,:)) < squeeze(eig_lambda(3,3,:));
% for i=1:length(acceptable_lambda_ratios)
%     if acceptable_lambda_ratios(i) == 0
%         lambda_small_switch = eig_lambda(3,3,i);
%         lambda_large_switch = eig_lambda(1,1,i);
%
%         eig_lambda(1,1,i) = lambda_small_switch;
%         eig_lambda(3,3,i) = lambda_large_switch;
%
%         vectors_small_switch = eig_vectors(:,3,i);
%         vectors_large_switch = eig_vectors(:,1,i);
%
%         eig_vectors(:,1,i) = vectors_small_switch;
%         eig_vectors(:,3,i) = -[vectors_large_switch(1), vectors_large_switch(2), vectors_large_switch(3)];
%
%                 % eig_vectors(:,:,i) = NaN;
%     else
%         %eig_vectors(:,:,i) = eig_vectors(:,:,i);
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%Sliding Window for MVAB normals%%%%%%%%%%%%%%%%%%%%%
%Same window range as MDD
t_0_range = end_index - start_index+1;
%srvy has 16 points per second, so 64 points is 4 seconds total
if MVAB_window_size == 0
    MVAB_window_size = t_0_range;
end

%initialize loop matrices
MVAB_lambda_all=zeros(3,t_0_range);
MVAB_normal_all=zeros(3,t_0_range);
MVAB_time=zeros(1,t_0_range);

for i=1:t_0_range
    
    t_0 = mms1_fgm_btimedata_raw(i+start_index-1);
    left_index = i+start_index-MVAB_window_size/2-1;
    right_index = i+start_index+MVAB_window_size/2-1;
    
    bdata_scope = mms1_fgm_bdata_raw(left_index:right_index,:);
    
    [output, l, v] = mvab(bdata_scope(:,1:3));
    
    MVAB_lambda_all(:,i) = l;
    MVAB_normal_all(:,i) = v(:,mvamm);%Minimum
    
    MVAB_time(1,i)=t_0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

title_name = 'MMS MDDB&MVAB Normals';
title(title_name, 'FontSize', 18, 'FontWeight', 'normal')

set(f1,'XMinorTick','on','YMinorTick','on')
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
if fluctuate_mdd == 1
    mdd_nx = Unfluctuate(squeeze(eig_vectors(1,mddmm,:)));
    mdd_ny = Unfluctuate(squeeze(eig_vectors(2,mddmm,:)));
    mdd_nz = Unfluctuate(squeeze(eig_vectors(3,mddmm,:)));
elseif fluctuate_mdd == 0
    mdd_nx = (squeeze(eig_vectors(1,mddmm,:)));
    mdd_ny = (squeeze(eig_vectors(2,mddmm,:)));
    mdd_nz = (squeeze(eig_vectors(3,mddmm,:)));
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot MDD Normals%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f3 = subplot(num_plots,1,3);

%Maxmimum
% plot(mms1_fgm_btimedata,squeeze(eig_vectors(1,mddmm,:)))
% hold on
% plot(mms1_fgm_btimedata,squeeze(eig_vectors(2,mddmm,:)))
% plot(mms1_fgm_btimedata,squeeze(eig_vectors(3,mddmm,:)))
% hold off

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
%xlim([datetime(tstart,'ConvertFrom','datenum') datetime(tend,'ConvertFrom','datenum')]);

ylabel(strcat('MDD_',mddl),'FontSize', 14)


xlim([mms1_fgm_btimedata(1) mms1_fgm_btimedata(end)])

set(f3,'XMinorTick','on','YMinorTick','on','XTickLabel',[])

set(f3,'Position',[plot_pos(1), plot_pos(2)-2*plot_pos(4)-0.075, plot_pos(3), plot_pos(4)*1]);

l1= findobj(gcf, 'Type', 'Legend');
l1pos=get(l1(3), 'Position');
%set(l1(3),'Position',[plot_pos(1)+plot_pos(3)+0.075, plot_pos(2)+2*plot_pos(4)+0.15, l1pos(3), l1pos(4)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot MDD Dimensionality%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f4 = subplot(num_plots,1,4);

plot(mms1_fgm_btimedata,sqrt(squeeze(eig_lambda(1,1,:))./squeeze(eig_lambda(2,2,:))))
hold on
plot(mms1_fgm_btimedata,sqrt(squeeze(eig_lambda(2,2,:))/squeeze(eig_lambda(3,3,:))),'g')

%plot(mms1_fgm_btimedata,sqrt(squeeze(eig_lambda(3,3,:))/squeeze(eig_lambda(1,1,:))),'r')
%plot(mms1_fgm_btimedata,squeeze(eig_lambda(3,3,:)))
hold off

datetick
set(gca,'layer','top')
legend({'\lambda_{max}/\lambda_{mid}', '\lambda_{mid}/\lambda_{min}', '\lambda_{min}/\lambda_{max}'},'FontSize',10)
legend('boxoff');
legend('Location','eastoutside');
%xlim([datetime(tstart,'ConvertFrom','datenum') datetime(tend,'ConvertFrom','datenum')]);

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

if fluctuate_mva == 1
    
    MVAB_normal_all(1,:) = Unfluctuate(MVAB_normal_all(1,:));
    MVAB_normal_all(2,:) = Unfluctuate(MVAB_normal_all(2,:));
    MVAB_normal_all(3,:) = Unfluctuate(MVAB_normal_all(3,:));
elseif fluctuate_mva == 0
    MVAB_normal_all(1,:) = (MVAB_normal_all(1,:));
    MVAB_normal_all(2,:) = (MVAB_normal_all(2,:));
    MVAB_normal_all(3,:) = (MVAB_normal_all(3,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot MVAB Normals%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f6 = subplot(num_plots,1,6);

plot(mms1_fgm_btimedata,MVAB_normal_all')

datetick
set(gca,'layer','top')
legend({'N_x', 'N_y', 'N_z'},'FontSize',10)
legend('boxoff');
legend('Location','eastoutside');
%xlim([datetime(tstart,'ConvertFrom','datenum') datetime(tend,'ConvertFrom','datenum')]);

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
%plot(mms1_fgm_btimedata,MVAB_lambda_all(3,:))
hold off

datetick
set(gca,'layer','top')
%legend({'\lambda_{max}', '\lambda_{mid}', '\lambda_{min}'},'FontSize',12)
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Average Angle%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate average of each component for the MDDB Normal
mean_x = mean(mdd_nx,'omitnan');
mean_y = mean(mdd_ny,'omitnan');
mean_z = mean(mdd_nz,'omitnan');
% mean_x = mean(mdd_nx,'omitnan');
% mean_y = mean(mdd_ny,'omitnan');
% mean_z = mean(mdd_nz,'omitnan');
MDD_mean_normal = ([mean_x, mean_y, mean_z])./norm([mean_x, mean_y, mean_z])


%mean grad_b of all terms
AA = mean(squeeze(grad_B(1,1,:)));
BB = mean(squeeze(grad_B(1,2,:)));
CC = mean(squeeze(grad_B(1,3,:)));
DD = mean(squeeze(grad_B(2,1,:)));
EE = mean(squeeze(grad_B(2,2,:)));
FF = mean(squeeze(grad_B(2,3,:)));
GG = mean(squeeze(grad_B(3,1,:)));
HH = mean(squeeze(grad_B(3,2,:)));
II = mean(squeeze(grad_B(3,3,:)));

KK = [AA, BB, CC; DD, EE, FF; GG, HH, II];
[VV,DD] = eig(KK*KK');

[DD,order] = sort(diag(DD),'descend');
VV = VV(:,order);

Mean_grad_B_MDD_normal = VV(:,1);


%if initial conditions did not have provide for a MVAB normal, find it.
if MVA_normal == [0, 0, 0]
    %Normal from MVA
    [~,mvab_l,mvab_v] = mvab(mms1_fgm_bdata_interp);
    
    MVA_normal = ([mvab_v(1,mvamm), mvab_v(2,mvamm), mvab_v(3,mvamm)]);%For the center t_0
    MVA_min_normal = ([mvab_v(1,3), mvab_v(2,3), mvab_v(3,3)])
    MVA_values = (mvab_l)
end


%fix 180 degrees ambiguity for one direction
% MVA_ambig(1)=rad2deg(acos(dot(MVA_normal.*[-1,1,1],Mean_grad_B_MDD_normal)/(norm(MVA_normal)*norm(Mean_grad_B_MDD_normal))));
% MVA_ambig(2)=rad2deg(acos(dot(MVA_normal.*[1,-1,1],Mean_grad_B_MDD_normal)/(norm(MVA_normal)*norm(Mean_grad_B_MDD_normal))));
% MVA_ambig(3)=rad2deg(acos(dot(MVA_normal.*[1,1,-1],Mean_grad_B_MDD_normal)/(norm(MVA_normal)*norm(Mean_grad_B_MDD_normal))));
% MVA_ambig(4)=rad2deg(acos(dot(MVA_normal,Mean_grad_B_MDD_normal)/(norm(MVA_normal)*norm(Mean_grad_B_MDD_normal))));

MVA_ambig(1)=angle(MVA_normal.*[-1,1,1],Mean_grad_B_MDD_normal);
MVA_ambig(2)=angle(MVA_normal.*[1,-1,1],Mean_grad_B_MDD_normal);
MVA_ambig(3)=angle(MVA_normal.*[1,1,-1],Mean_grad_B_MDD_normal);
MVA_ambig(4)=angle(MVA_normal.*[1,1,1],Mean_grad_B_MDD_normal);



switch find(MVA_ambig==min(MVA_ambig))
    case 1
        MVA_normal = MVA_normal.*[-1,1,1]
    case 2
        MVA_normal = MVA_normal.*[1,-1,1]
    case 3
        MVA_normal = MVA_normal.*[1,1,-1]
    case 4
        MVA_normal = MVA_normal
end



% %fix 180 degrees ambiguity for one direction
% MVA_ambig(1)=rad2deg(acos(dot(MVA_normal.*[-1,1,1],MDD_mean_normal)/(norm(MVA_normal)*norm(MDD_mean_normal))));
% MVA_ambig(2)=rad2deg(acos(dot(MVA_normal.*[1,-1,1],MDD_mean_normal)/(norm(MVA_normal)*norm(MDD_mean_normal))));
% MVA_ambig(3)=rad2deg(acos(dot(MVA_normal.*[1,1,-1],MDD_mean_normal)/(norm(MVA_normal)*norm(MDD_mean_normal))));
% MVA_ambig(4)=rad2deg(acos(dot(MVA_normal,MDD_mean_normal)/(norm(MVA_normal)*norm(MDD_mean_normal))));
% 
% switch find(MVA_ambig==min(MVA_ambig))
%     case 1
%         MVA_normal = MVA_normal.*[-1,1,1]
%     case 2
%         MVA_normal = MVA_normal.*[1,-1,1]
%     case 3
%         MVA_normal = MVA_normal.*[1,1,-1]
%     case 4
%         MVA_normal = MVA_normal
% end


%This is the mean MDD dotted with the center MVA
MDD_dot_MVA_center = min(MVA_ambig)
% MVA_cross_MDD_center = cross(MVA_normal,MDD_mean_normal)
% MVA_min_normal;
% MVA_cross_MDD_dot_MVA_min = rad2deg(acos(dot(MVA_cross_MDD_center,MVA_min_normal)))

%%%Calculation of Angle at every point
MDDB_normal = [mdd_nx, mdd_ny, mdd_nz]; %MDDB normal at every point
MVAB_normal_all = MVAB_normal_all'; %MVAB normals at every point

MVA_dot_MDD_angle_all = zeros(length(rx),1);
MVA_dot_MDD_center_all = zeros(length(rx),1);
MVA_ambig_dot_center= zeros(length(rx),4);

for i=1:length(rx)
    %fix 180 degrees ambiguity for one direction
    MVA_ambig_dot(1)=angle(MVAB_normal_all(i,:).*[-1,1,1],MDD_mean_normal);
    MVA_ambig_dot(2)=angle(MVAB_normal_all(i,:).*[1,-1,1],MDD_mean_normal);
    MVA_ambig_dot(3)=angle(MVAB_normal_all(i,:).*[1,1,-1],MDD_mean_normal);
    MVA_ambig_dot(4)=angle(MVAB_normal_all(i,:),MDD_mean_normal);
    %     MVA_dot_MDD_angle_all(i,1) = rad2deg(acos(dot(MVAB_normal_all(i,:),MDDB_normal(i,:)))/(norm(MVAB_normal_all(i,:))*norm(MDDB_normal(i,:))));
    MVA_dot_MDD_angle_all(i,1) = min(MVA_ambig_dot); %find the minimum dot product value.
    
    
%     MVA_ambig_dot_center(i,1) = rad2deg(acos(dot(MVA_min_normal.*[-1,1,1],MDDB_normal(i,:)))/(norm(MVA_min_normal)*norm(MDDB_normal(i,:))));
%     MVA_ambig_dot_center(i,2) = rad2deg(acos(dot(MVA_min_normal.*[1,-1,1],MDDB_normal(i,:)))/(norm(MVA_min_normal)*norm(MDDB_normal(i,:))));
%     MVA_ambig_dot_center(i,3) = rad2deg(acos(dot(MVA_min_normal.*[1,1,-1],MDDB_normal(i,:)))/(norm(MVA_min_normal)*norm(MDDB_normal(i,:))));
%     MVA_ambig_dot_center(i,4) = rad2deg(acos(dot(MVA_min_normal,MDDB_normal(i,:)))/(norm(MVA_min_normal)*norm(MDDB_normal(i,:))));
    
    %     MVA_dot_MDD_center_all(i,1) = rad2deg(acos(dot(MVA_min_normal,MDDB_normal(i,:)))/(norm(MVA_min_normal)*norm(MDDB_normal(i,:))));
end
% ambig_column = find(mean(MVA_ambig_dot_center)==min(mean(MVA_ambig_dot_center)));
% MVA_dot_MDD_center_all(:,1) = MVA_ambig_dot_center(:,ambig_column);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot Angles%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f8 = subplot(num_plots,1,8);

plot(mms1_fgm_btimedata,MVA_dot_MDD_angle_all)
% hold on
% plot(mms1_fgm_btimedata,MVA_dot_MDD_center_all)
% hold off
box on
datetick('x','SS.FFF','keeplimits')
set(gca,'layer','top')
legend({'Sliding MVA'},'FontSize',10)
legend('boxoff');
legend('Location','eastoutside');
%xlim([datetime(tstart,'ConvertFrom','datenum') datetime(tend,'ConvertFrom','datenum')]);
ylabel({'Angle'; '[Deg]'},'FontSize', 14)

xlim([mms1_fgm_btimedata(1) mms1_fgm_btimedata(end)])
% ylim([-45 45])
%line(xlim(), [0,0], 'LineWidth', 0.5, 'Color', 'k','LineStyle','--');
%grid on
% grid on;
set(f8,'XMinorTick','on','YMinorTick','on')

set(f8,'Position',[plot_pos(1), plot_pos(2)-7*plot_pos(4)-0.2, plot_pos(3), plot_pos(4)*1]);

%
% l1= findobj(gcf, 'Type', 'Legend');
% l1pos=get(l1(5), 'Position');
% set(l1(8),'Position',[plot_pos(1)+plot_pos(3)+0.1, plot_pos(2)+4*plot_pos(4), l1pos(3), l1pos(4)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%Text on Plotot%%%%%%%%%%%%%%%%%%%%%%%
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
    'String',strcat('(Average MDD)\cdot(Centered MVA) Angle= ', ' ', num2str(MDD_dot_MVA_center,'%2.2f'),char(176)),...
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


