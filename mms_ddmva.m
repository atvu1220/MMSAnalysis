%Andrew Vu 9/25/18
%data=spdfcdfread(filename);
%datainfo=spdfcdfinfo(filename);
%clear
addpath ~/Library/'Mobile Documents'/com~apple~CloudDocs/
addpath ~/data/

tic
close all
clear
figure('Position',[1 1 1600 600])
mms_directory = '/Users/andrewvu/data/mms/';
probe_num = '1';
num_plots = 5;
data_type = 'brst';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Data Retrieval%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
event_start = '2018-03-12 07:35:50.000';
event_end = '2018-03-12 07:36:30.000';
date_start = '2018-03-12 07:36:04.000';
date_end = '2018-03-12 07:36:06.000';
MVA_normal = [0,0,0];

formatIn='yyyy-mm-dd HH:MM:SS.FFF';
tstart = datenum(date_start,formatIn);
tend = datenum(date_end,formatIn);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[mms1_fgm_btimedata_raw, mms1_fgm_bdata_raw, mms1_fgm_rtimedata_raw, mms1_fgm_rdata_raw] = load_fgm(date_start,1,data_type);
[mms2_fgm_btimedata_raw, mms2_fgm_bdata_raw, mms2_fgm_rtimedata_raw, mms2_fgm_rdata_raw] = load_fgm(date_start,2,data_type);
[mms3_fgm_btimedata_raw, mms3_fgm_bdata_raw, mms3_fgm_rtimedata_raw, mms3_fgm_rdata_raw] = load_fgm(date_start,3,data_type);
[mms4_fgm_btimedata_raw, mms4_fgm_bdata_raw, mms4_fgm_rtimedata_raw, mms4_fgm_rdata_raw] = load_fgm(date_start,4,data_type);



%%%%%%%%%%%%%%%%%%%%Interpolate position and bdata with bdata from mms1 %%%%%%%%%%%%%%%%%%%

%Interpolate Sata
mms1_fgm_rdata_interp(:,1) = interp1(mms1_fgm_rtimedata_raw,mms1_fgm_rdata_raw(:,1),mms1_fgm_btimedata_raw);
mms1_fgm_rdata_interp(:,2) = interp1(mms1_fgm_rtimedata_raw,mms1_fgm_rdata_raw(:,2),mms1_fgm_btimedata_raw);
mms1_fgm_rdata_interp(:,3) = interp1(mms1_fgm_rtimedata_raw,mms1_fgm_rdata_raw(:,3),mms1_fgm_btimedata_raw);

mms2_fgm_bdata_interp(:,1) = interp1(mms2_fgm_btimedata_raw,mms2_fgm_bdata_raw(:,1),mms1_fgm_btimedata_raw);
mms2_fgm_bdata_interp(:,2) = interp1(mms2_fgm_btimedata_raw,mms2_fgm_bdata_raw(:,2),mms1_fgm_btimedata_raw);
mms2_fgm_bdata_interp(:,3) = interp1(mms2_fgm_btimedata_raw,mms2_fgm_bdata_raw(:,3),mms1_fgm_btimedata_raw);

mms2_fgm_rdata_interp(:,1) = interp1(mms2_fgm_rtimedata_raw,mms2_fgm_rdata_raw(:,1),mms1_fgm_btimedata_raw);
mms2_fgm_rdata_interp(:,2) = interp1(mms2_fgm_rtimedata_raw,mms2_fgm_rdata_raw(:,2),mms1_fgm_btimedata_raw);
mms2_fgm_rdata_interp(:,3) = interp1(mms2_fgm_rtimedata_raw,mms2_fgm_rdata_raw(:,3),mms1_fgm_btimedata_raw);

mms3_fgm_bdata_interp(:,1) = interp1(mms3_fgm_btimedata_raw,mms3_fgm_bdata_raw(:,1),mms1_fgm_btimedata_raw);
mms3_fgm_bdata_interp(:,2) = interp1(mms3_fgm_btimedata_raw,mms3_fgm_bdata_raw(:,2),mms1_fgm_btimedata_raw);
mms3_fgm_bdata_interp(:,3) = interp1(mms3_fgm_btimedata_raw,mms3_fgm_bdata_raw(:,3),mms1_fgm_btimedata_raw);

mms3_fgm_rdata_interp(:,1) = interp1(mms3_fgm_rtimedata_raw,mms3_fgm_rdata_raw(:,1),mms1_fgm_btimedata_raw);
mms3_fgm_rdata_interp(:,2) = interp1(mms3_fgm_rtimedata_raw,mms3_fgm_rdata_raw(:,2),mms1_fgm_btimedata_raw);
mms3_fgm_rdata_interp(:,3) = interp1(mms3_fgm_rtimedata_raw,mms3_fgm_rdata_raw(:,3),mms1_fgm_btimedata_raw);

mms4_fgm_bdata_interp(:,1) = interp1(mms4_fgm_btimedata_raw,mms4_fgm_bdata_raw(:,1),mms1_fgm_btimedata_raw);
mms4_fgm_bdata_interp(:,2) = interp1(mms4_fgm_btimedata_raw,mms4_fgm_bdata_raw(:,2),mms1_fgm_btimedata_raw);
mms4_fgm_bdata_interp(:,3) = interp1(mms4_fgm_btimedata_raw,mms4_fgm_bdata_raw(:,3),mms1_fgm_btimedata_raw);

mms4_fgm_rdata_interp(:,1) = interp1(mms4_fgm_rtimedata_raw,mms4_fgm_rdata_raw(:,1),mms1_fgm_btimedata_raw);
mms4_fgm_rdata_interp(:,2) = interp1(mms4_fgm_rtimedata_raw,mms4_fgm_rdata_raw(:,2),mms1_fgm_btimedata_raw);
mms4_fgm_rdata_interp(:,3) = interp1(mms4_fgm_rtimedata_raw,mms4_fgm_rdata_raw(:,3),mms1_fgm_btimedata_raw);


%Find the start and end limits of the event in the data, index of, in
%datenum format
event_start_index = find(mms1_fgm_btimedata_raw > datenum(event_start,formatIn), 1)-1;
event_end_index = find(mms1_fgm_btimedata_raw > datenum(event_end,formatIn), 1);

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


%group compoinents of B and r together for faster loops.
bx=[mms1_fgm_bdata_interp(:,1) mms2_fgm_bdata_interp(:,1) mms3_fgm_bdata_interp(:,1) mms4_fgm_bdata_interp(:,1)];
by=[mms1_fgm_bdata_interp(:,2) mms2_fgm_bdata_interp(:,2) mms3_fgm_bdata_interp(:,2) mms4_fgm_bdata_interp(:,2)];
bz=[mms1_fgm_bdata_interp(:,3) mms2_fgm_bdata_interp(:,3) mms3_fgm_bdata_interp(:,3) mms4_fgm_bdata_interp(:,3)];

rx=[mms1_fgm_rdata_interp(:,1) mms2_fgm_rdata_interp(:,1) mms3_fgm_rdata_interp(:,1) mms4_fgm_rdata_interp(:,1)];
ry=[mms1_fgm_rdata_interp(:,2) mms2_fgm_rdata_interp(:,2) mms3_fgm_rdata_interp(:,2) mms4_fgm_rdata_interp(:,2)];
rz=[mms1_fgm_rdata_interp(:,3) mms2_fgm_rdata_interp(:,3) mms3_fgm_rdata_interp(:,3) mms4_fgm_rdata_interp(:,3)];

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
    
    
    [eig_l,order] = sort(diag(eig_l),'ascend');
    eig_v = eig_v(:,order);
    
    eig_lambda(:,:,i) = diag(eig_l);
    eig_vectors(:,:,i) = eig_v;
end

%constract multidimensional array for the Gradient of B at every time
% for i=1:length(rx)
%     grad_B(:,:,i) = gradB(i,Summ_probes,bx,by,bz,rx,ry,rz,R_inv);
% end


%Find the largest eigenvalue, corresponds to largest G.
Gmax_index = find(squeeze(eig_lambda(3,3,:)) == max(squeeze(eig_lambda(3,3,:))));


%Calculate the average of G in the time interval, Denton et al. 2010
Gxx = mean(squeeze(grad_B(1,1,Gmax_index-0:Gmax_index+0)));
Gyx = mean(squeeze(grad_B(1,2,Gmax_index-0:Gmax_index+0)));
Gzx = mean(squeeze(grad_B(1,3,Gmax_index-0:Gmax_index+0)));

Gxy = mean(squeeze(grad_B(2,1,Gmax_index-0:Gmax_index+0)));
Gyy = mean(squeeze(grad_B(2,2,Gmax_index-0:Gmax_index+0)));
Gzy = mean(squeeze(grad_B(2,3,Gmax_index-0:Gmax_index+0)));

Gxz = mean(squeeze(grad_B(3,1,Gmax_index-0:Gmax_index+0)));
Gyz = mean(squeeze(grad_B(3,2,Gmax_index-0:Gmax_index+0)));
Gzz = mean(squeeze(grad_B(3,3,Gmax_index-0:Gmax_index+0)));

%average grad_B in the time interval of largest G, or largest eigenvalue
G_0 = abs([Gxx, Gyx, Gzx;...
    Gxy, Gyy, Gzy;...
    Gxz, Gyz, Gzz]);

%Calculate the Eigenvalues and Eigenvectors of L = GG', wherewe have
%subtracted the average value of G in the largest G interval
for i=1:length(rx)
    L = (grad_B(:,:,i)-G_0)*(grad_B(:,:,i)-G_0)';
        [eig_v, eig_l] = eig(L);
    
    
    [eig_l,order] = sort(diag(eig_l),'ascend');
    eig_v = eig_v(:,order);
    
    eig_lambda(:,:,i) = diag(eig_l);
    eig_vectors(:,:,i) = eig_v;
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
%arbitrary scope set dt_0
window_range = end_index - start_index+1;
%srvy has 16 points per second, so 64 points is 4 seconds total
window_size = length(rx);


%initialize loop matrices
MVAB_lambda=zeros(3,window_range);
MVAB_normal=zeros(3,window_range);
MVAB_time=zeros(1,window_range);

for i=1:window_range

    t_0 = mms1_fgm_btimedata_raw(i+start_index-1);
    left_index = i+start_index-window_size/2-1;
    right_index = i+start_index+window_size/2-1;
    
    bdata_scope = mms1_fgm_bdata_raw(left_index:right_index,:);
    
    [output, l, v] = mvab(bdata_scope(:,1:3));
    
    MVAB_lambda(:,i) = l;
    MVAB_normal(:,i) = v(:,3);
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
legend({'B_x', 'B_y', 'B_z','B_{tot}'},'FontSize',14)
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
set(l1,'Position',[plot_pos(1)+plot_pos(3)+0.025, plot_pos(2)+plot_pos(4)/2-l1pos(4)/2, l1pos(3), l1pos(4)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot Lambdas%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f2 = subplot(num_plots,1,2);

plot(mms1_fgm_btimedata,mms1_fgm_bdata_raw(start_index:end_index,1:4))
box on
datetick
set(gca,'layer','top')
legend({'B_x', 'B_y', 'B_z','B_{tot}'},'FontSize',14)
legend('boxoff');
legend('Location','eastoutside');
datetick('x','HH:MM:ss')
xlim([mms1_fgm_btimedata(1) mms1_fgm_btimedata(end)])
ylabel({'B';'[nT]'},'FontSize', 14)

xlim([mms1_fgm_btimedata(1) mms1_fgm_btimedata(end)])

set(f2,'XMinorTick','on','YMinorTick','on','XTickLabel',[])

set(f2,'Position',[plot_pos(1), plot_pos(2)-1*plot_pos(4)-0.05, plot_pos(3), plot_pos(4)]);

l1= findobj(gcf, 'Type', 'Legend');
l1pos=get(l1(2), 'Position');
set(l1(2),'Position',[plot_pos(1)+plot_pos(3)+0.025, plot_pos(2)+1*plot_pos(4)-0.15, l1pos(3), l1pos(4)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot Normals%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f3 = subplot(num_plots,1,3);

plot(mms1_fgm_btimedata,squeeze(eig_vectors(1,3,:)))
hold on
plot(mms1_fgm_btimedata,squeeze(eig_vectors(2,3,:)))
plot(mms1_fgm_btimedata,squeeze(eig_vectors(3,3,:)))
hold off

box on
datetick
set(gca,'layer','top')
legend({'N_x', 'N_y', 'N_z'},'FontSize',14)
legend('boxoff');
legend('Location','eastoutside');
%xlim([datetime(tstart,'ConvertFrom','datenum') datetime(tend,'ConvertFrom','datenum')]);
ylabel('N_{MDD}','FontSize', 14)

xlim([mms1_fgm_btimedata(1) mms1_fgm_btimedata(end)])

set(f3,'XMinorTick','on','YMinorTick','on','XTickLabel',[])

set(f3,'Position',[plot_pos(1), plot_pos(2)-2*plot_pos(4)-0.075, plot_pos(3), plot_pos(4)*1]);

l1= findobj(gcf, 'Type', 'Legend');
l1pos=get(l1(3), 'Position');
%set(l1(3),'Position',[plot_pos(1)+plot_pos(3)+0.075, plot_pos(2)+2*plot_pos(4)+0.15, l1pos(3), l1pos(4)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot Normals%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f4 = subplot(num_plots,1,4);

plot(mms1_fgm_btimedata,MVAB_normal')

datetick
set(gca,'layer','top')
legend({'N_x', 'N_y', 'N_z'},'FontSize',14)
legend('boxoff');
legend('Location','eastoutside');
%xlim([datetime(tstart,'ConvertFrom','datenum') datetime(tend,'ConvertFrom','datenum')]);
ylabel('N_{MVAB}','FontSize', 14)

xlim([mms1_fgm_btimedata(1) mms1_fgm_btimedata(end)])

set(f4,'XMinorTick','on','YMinorTick','on','XTickLabel',[])

set(f4,'Position',[plot_pos(1), plot_pos(2)-3*plot_pos(4)-0.10, plot_pos(3), plot_pos(4)*1]);

% l1= findobj(gcf, 'Type', 'Legend');
% l1pos=get(l1(4), 'Position');
% set(l1(4),'Position',[plot_pos(1)+plot_pos(3)+0.1, plot_pos(2)+3*plot_pos(4), l1pos(3), l1pos(4)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% %Calculate average of each component for the Normal
MDDB_Nx = (squeeze(eig_vectors(1,3,:)));
MDDB_Ny = (squeeze(eig_vectors(2,3,:)));
MDDB_Nz = (squeeze(eig_vectors(3,3,:)));

MDDB_normal = [MDDB_Nx, MDDB_Ny, MDDB_Nz];
%%%Calculation of Angle at every step
MVAB_normal = MVAB_normal';
MVA_dot_MDD_angle = zeros(length(rx),1);
for i=1:length(rx)
    disp(i)
MVA_dot_MDD_angle(i,1) = rad2deg(acos(dot(MVAB_normal(i,:),MDDB_normal(i,:)))/(norm(MVAB_normal(i,:))*norm(MDDB_normal(i,:))));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot Normals%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f5 = subplot(num_plots,1,5);








plot(mms1_fgm_btimedata,MVA_dot_MDD_angle)

box on
datetick
set(gca,'layer','top')
% legend({'N_x', 'N_y', 'N_z'},'FontSize',14)
% legend('boxoff');
% legend('Location','eastoutside');
%xlim([datetime(tstart,'ConvertFrom','datenum') datetime(tend,'ConvertFrom','datenum')]);
ylabel({'Angle'; '[Deg]'},'FontSize', 14)

xlim([mms1_fgm_btimedata(1) mms1_fgm_btimedata(end)])

set(f5,'XMinorTick','on','YMinorTick','on')

set(f5,'Position',[plot_pos(1), plot_pos(2)-4*plot_pos(4)-0.125, plot_pos(3), plot_pos(4)*1]);
% 
% l1= findobj(gcf, 'Type', 'Legend');
% l1pos=get(l1(5), 'Position');
% set(l1(5),'Position',[plot_pos(1)+plot_pos(3)+0.1, plot_pos(2)+4*plot_pos(4), l1pos(3), l1pos(4)]);



% mean_x = mean(Ax,'omitnan');
% mean_y = mean(Ay,'omitnan');
% mean_z = mean(Az,'omitnan');
% 
% MMD_normal = [mean_x, mean_y, mean_z]
% 
% if initial conditions did not have provide for a normal, find it.
% if MVA_normal == [0, 0, 0]
%     %Normal from MVA
%     [~,mvab_l,mvab_v] = mvab(mms1_fgm_bdata_interp);
%     MVA_normal = [mvab_v(1,3), mvab_v(2,3), mvab_v(3,3)]
%     MVA_values = (mvab_l)
% end
% 
% MVA_dot_MMD_angle = rad2deg(acos(dot(MVA_normal,MMD_normal)/(norm(MVA_normal)*norm(MMD_normal))))


annotation('textbox',[plot_pos(1), plot_pos(2)-5*plot_pos(4)-0.06, plot_pos(3), plot_pos(4)*1],...
    'String',{date_start(1:10)},...
    'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

orient(gcf,'landscape')

plot_name =  strcat('mms','_mdd_mva_angles',...
    date_start(1:19),'.pdf');

print(gcf, '-dpdf', '-opengl',plot_name,'-fillpage');

movefile(plot_name, '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis')



toc

function [grad_B] = gradB(time,Summ_probes,probes,bx,by,bz,rx,ry,rz,R_inv)
    
    %dbx/dx
    brx = sum((bx(time,Summ_probes(:,1))-bx(time,Summ_probes(:,2)))...
        .*(rx(time,Summ_probes(:,1))-rx(time,Summ_probes(:,2))));
    bry = sum((bx(time,Summ_probes(:,1))-bx(time,Summ_probes(:,2)))...
        .*(ry(time,Summ_probes(:,1))-ry(time,Summ_probes(:,2))));
    brz = sum((bx(time,Summ_probes(:,1))-bx(time,Summ_probes(:,2)))...
        .*(rz(time,Summ_probes(:,1))-rz(time,Summ_probes(:,2))));
    
%         brx=sum(bx(time,probes).*rx(time,probes));
%         bry=sum(bx(time,probes).*ry(time,probes));
%         brz=sum(bx(time,probes).*rz(time,probes));
    
    br = [brx bry brz];
    dBxx = br*R_inv(:,1,time)/4^2;
    
    %dby/dx
    brx = sum((by(time,Summ_probes(:,1))-by(time,Summ_probes(:,2)))...
        .*(rx(time,Summ_probes(:,1))-rx(time,Summ_probes(:,2))));
    bry = sum((by(time,Summ_probes(:,1))-by(time,Summ_probes(:,2)))...
        .*(ry(time,Summ_probes(:,1))-ry(time,Summ_probes(:,2))));
    brz = sum((by(time,Summ_probes(:,1))-by(time,Summ_probes(:,2)))...
        .*(rz(time,Summ_probes(:,1))-rz(time,Summ_probes(:,2))));
    
%         brx=sum(by(time,probes).*rx(time,probes));
%         bry=sum(by(time,probes).*ry(time,probes));
%         brz=sum(by(time,probes).*rz(time,probes));
    
    br = [brx bry brz];
    dByx = br*R_inv(:,1,time)/4^2;
    
    %dbz/dx
    brx = sum((bz(time,Summ_probes(:,1))-bz(time,Summ_probes(:,2)))...
        .*(rx(time,Summ_probes(:,1))-rx(time,Summ_probes(:,2))));
    bry = sum((bz(time,Summ_probes(:,1))-bz(time,Summ_probes(:,2)))...
        .*(ry(time,Summ_probes(:,1))-ry(time,Summ_probes(:,2))));
    brz = sum((bz(time,Summ_probes(:,1))-bz(time,Summ_probes(:,2)))...
        .*(rz(time,Summ_probes(:,1))-rz(time,Summ_probes(:,2))));
    
%         brx=sum(bz(time,probes).*rx(time,probes));
%         bry=sum(bz(time,probes).*ry(time,probes));
%         brz=sum(bz(time,probes).*rz(time,probes));
    
    br = [brx bry brz];
    dBzx = br*R_inv(:,1,time)/4^2;
    
    
    %dbx/dy
    brx = sum((bx(time,Summ_probes(:,1))-bx(time,Summ_probes(:,2)))...
        .*(rx(time,Summ_probes(:,1))-rx(time,Summ_probes(:,2))));
    bry = sum((bx(time,Summ_probes(:,1))-bx(time,Summ_probes(:,2)))...
        .*(ry(time,Summ_probes(:,1))-ry(time,Summ_probes(:,2))));
    brz = sum((bx(time,Summ_probes(:,1))-bx(time,Summ_probes(:,2)))...
        .*(rz(time,Summ_probes(:,1))-rz(time,Summ_probes(:,2))));
    
        brx=sum(bx(time,probes).*rx(time,probes));
        bry=sum(bx(time,probes).*ry(time,probes));
        brz=sum(bx(time,probes).*rz(time,probes));
    
    br = [brx bry brz];
    dBxy = br*R_inv(:,2,time)/4^2;
    
    %dby/dy
    brx = sum((by(time,Summ_probes(:,1))-by(time,Summ_probes(:,2)))...
        .*(rx(time,Summ_probes(:,1))-rx(time,Summ_probes(:,2))));
    bry = sum((by(time,Summ_probes(:,1))-by(time,Summ_probes(:,2)))...
        .*(ry(time,Summ_probes(:,1))-ry(time,Summ_probes(:,2))));
    brz = sum((by(time,Summ_probes(:,1))-by(time,Summ_probes(:,2)))...
        .*(rz(time,Summ_probes(:,1))-rz(time,Summ_probes(:,2))));
    
    
%         brx=sum(by(time,probes).*rx(time,probes));
%         bry=sum(by(time,probes).*ry(time,probes));
%         brz=sum(by(time,probes).*rz(time,probes));
    
    br = [brx bry brz];
    dByy = br*R_inv(:,2,time)/4^2;
    
    %dbz/dy
    brx = sum((bz(time,Summ_probes(:,1))-bz(time,Summ_probes(:,2)))...
        .*(rx(time,Summ_probes(:,1))-rx(time,Summ_probes(:,2))));
    bry = sum((bz(time,Summ_probes(:,1))-bz(time,Summ_probes(:,2)))...
        .*(ry(time,Summ_probes(:,1))-ry(time,Summ_probes(:,2))));
    brz = sum((bz(time,Summ_probes(:,1))-bz(time,Summ_probes(:,2)))...
        .*(rz(time,Summ_probes(:,1))-rz(time,Summ_probes(:,2))));
    
%         brx=sum(bz(time,probes).*rx(time,probes));
%         bry=sum(bz(time,probes).*ry(time,probes));
%         brz=sum(bz(time,probes).*rz(time,probes));
   
    br = [brx bry brz];
    dBzy = br*R_inv(:,2,time)/4^2;
    
    
    %dbx/dz
    brx = sum((bx(time,Summ_probes(:,1))-bx(time,Summ_probes(:,2)))...
        .*(rx(time,Summ_probes(:,1))-rx(time,Summ_probes(:,2))));
    bry = sum((bx(time,Summ_probes(:,1))-bx(time,Summ_probes(:,2)))...
        .*(ry(time,Summ_probes(:,1))-ry(time,Summ_probes(:,2))));
    brz = sum((bx(time,Summ_probes(:,1))-bx(time,Summ_probes(:,2)))...
        .*(rz(time,Summ_probes(:,1))-rz(time,Summ_probes(:,2))));
    
%         brx=sum(bx(time,probes).*rx(time,probes));
%         bry=sum(bx(time,probes).*ry(time,probes));
%         brz=sum(bx(time,probes).*rz(time,probes));
    
    br = [brx bry brz];
    dBxz = br*R_inv(:,3,time)/4^2;
    
    %dby/dz
    brx = sum((by(time,Summ_probes(:,1))-by(time,Summ_probes(:,2)))...
        .*(rx(time,Summ_probes(:,1))-rx(time,Summ_probes(:,2))));
    bry = sum((by(time,Summ_probes(:,1))-by(time,Summ_probes(:,2)))...
        .*(ry(time,Summ_probes(:,1))-ry(time,Summ_probes(:,2))));
    brz = sum((by(time,Summ_probes(:,1))-by(time,Summ_probes(:,2)))...
        .*(rz(time,Summ_probes(:,1))-rz(time,Summ_probes(:,2))));
    
%         brx=sum(by(time,probes).*rx(time,probes));
%         bry=sum(by(time,probes).*ry(time,probes));
%         brz=sum(by(time,probes).*rz(time,probes));
    
    br = [brx bry brz];
    dByz = br*R_inv(:,3,time)/4^2;
    
    %dbz/dz
    brx = sum((bz(time,Summ_probes(:,1))-bz(time,Summ_probes(:,2)))...
        .*(rx(time,Summ_probes(:,1))-rx(time,Summ_probes(:,2))));
    bry = sum((bz(time,Summ_probes(:,1))-bz(time,Summ_probes(:,2)))...
        .*(ry(time,Summ_probes(:,1))-ry(time,Summ_probes(:,2))));
    brz = sum((bz(time,Summ_probes(:,1))-bz(time,Summ_probes(:,2)))...
        .*(rz(time,Summ_probes(:,1))-rz(time,Summ_probes(:,2))));
    
%         brx=sum(bz(time,probes).*rx(time,probes));
%         bry=sum(bz(time,probes).*ry(time,probes));
%         brz=sum(bz(time,probes).*rz(time,probes));
    
    br = [brx bry brz];
    dBzz = br*R_inv(:,3,time)/4^2;
    
    %form the gradient matrix
    grad_B = [dBxx, dByx, dBzx;...
        dBxy, dByy, dBzy;...
        dBxz, dByz, dBzz];
    
    %Constraint for delB = 0 ;
    lambda = trace(grad_B)/(trace(R_inv(:,:,time)));
    
    grad_B = grad_B - lambda*R_inv(:,:,time);
end

