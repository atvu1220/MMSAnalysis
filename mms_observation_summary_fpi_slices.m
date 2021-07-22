%MATLAB MMS Observation Summary with FPI Slices from IDL
%Andrew Vu 4/04/2019
%data=spdfcdfread(filename);
%datainfo=spdfcdfinfo(filename);
clear
close all
figure('Position',[1 1 1200 600])
co = [0 0 1;
    0 1 0;
    1 0 0;
    0 0 0];
set(gcf,'defaultAxesColorOrder',co)
set(gcf,'color','w');
plot_gap = 1.00;
cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis'
mms_directory = '/Users/andrewvu/data/mms/';
probe_num = '1';
num_plots = 5;
data_type = 'brst';
Rotation = 'BE';
probe = 1;

EventNumber = 176;
if EventNumber == 0 
    trange = {'2017-11-22 04:56:20.000','2017-11-22 04:57:50.000'};
    event_start = cell2mat(trange(1));
    event_end = cell2mat(trange(2));
else
[~, ~,~,~,...
    event_start,event_end,...
    ~, ~,...
    ~, ~,...
    ~, ~,...
    ~, ~,...
    ~, ~,...
    ~, ~] = get_eventTimes(EventNumber);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%Load Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%load fgm data
[mms1_fgm_timedata_raw, mms1_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,probe,data_type);
% [mms2_fgm_timedata_raw, mms2_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,2,data_type);
% [mms3_fgm_timedata_raw, mms3_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,3,data_type);
% [mms4_fgm_timedata_raw, mms4_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,4,data_type);

%Load Mec Data from Epheremis
[mms1_mec_timedata_raw, mms1_mec_rdata_raw] = load_mec(event_start,probe,'srvy');
% [mms2_mec_timedata_raw, mms2_mec_rdata_raw] = load_mec(event_start,2,'srvy');
% [mms3_mec_timedata_raw, mms3_mec_rdata_raw] = load_mec(event_start,3,'srvy');
% [mms4_mec_timedata_raw, mms4_mec_rdata_raw] = load_mec(event_start,4,'srvy');

%Load FPI_e
[fpi_e_timedata,fpi_e_ndata,fpi_e_vdata,fpi_e_tparadata,fpi_e_tperpdata,...
    fpi_e_edata,fpi_e_espectdata] = load_fpi(event_start,event_end,probe_num,data_type,'e');
%Load FPI_i
[fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_tparadata,fpi_i_tperpdata,...
    fpi_i_edata,fpi_i_espectdata] = load_fpi(event_start,event_end,probe_num,data_type,'i');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%Plot Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Time Indices
%Find the start and end limits of the event in the data
formatIn='yyyy-mm-dd HH:MM:SS.FFF';
tstart = datenum(event_start,formatIn);
tend = datenum(event_end,formatIn);

B_start_index = find(mms1_fgm_timedata_raw >= tstart, 1);
B_end_index = find(mms1_fgm_timedata_raw >= tend, 1);

i_start_index = find(fpi_i_timedata > tstart, 1);
i_end_index = find(fpi_i_timedata > tend, 1);

e_start_index = find(fpi_e_timedata > tstart, 1);
e_end_index = find(fpi_e_timedata > tend, 1);

%crop data
fgm_timedata = mms1_fgm_timedata_raw(B_start_index:B_end_index,1);
fpi_i_timedata = fpi_i_timedata(i_start_index:i_end_index,1);
fpi_e_timedata = fpi_e_timedata(e_start_index:e_end_index,1);


%% Magnetic Field
plot_order = 1;
fgm_bdata = mms1_fgm_bdata_raw(B_start_index:B_end_index,:);

Bplot = subplot(num_plots,2,plot_order);
plot(fgm_timedata,fgm_bdata(:,:),'LineWidth',1)
legend({'B_x', 'B_y', 'B_z','B_t'},'FontSize',8)
legend('boxoff')
legend('Location','eastoutside','AutoUpdate','off')
colormap('winter');
datetick
xlim([fgm_timedata(1) fgm_timedata(end)])
ylabel({'B';'[nT]'},'FontSize', 14)
set(gca, 'XTickLabel', [],'XMinorTick','on','YMinorTick','on','linewidth',1.25)

title('MMS1 Observatory FPI Summary', 'FontSize', 18, 'FontWeight', 'normal')
plot_pos = get(gca,'Position');
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order/2-1+0.5), plot_pos(3), plot_pos(4)]);
plot_order = plot_order+2;
%% Density
fpi_i_ndata = fpi_i_ndata(i_start_index:i_end_index,:);
%fpi_e_ndata = fpi_e_ndata(e_start_index:e_end_index,:);

nplot = subplot(num_plots,2,plot_order);
plot(fpi_i_timedata,fpi_i_ndata,'LineWidth',1)
% hold on
% plot(fpi_e_timedata,fpi_e_ndata,'LineWidth',1)
%legend({'N^i','N^e'},'FontSize',10)
legend({'N^i'},'FontSize',10)
legend('boxoff')
legend('Location','eastoutside','AutoUpdate','off')
colormap('winter');
datetick
xlim([fpi_i_timedata(1) fpi_i_timedata(end)])
ylabel({'n';'[cm^{-3}]'},'FontSize', 14)
set(gca, 'XTickLabel', [],'XMinorTick','on','YMinorTick','on','linewidth',1.25)
hold off

set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order/2-1+0.5), plot_pos(3), plot_pos(4)]);
plot_order = plot_order+2;
%% Velocity
fpi_vxdata=(fpi_i_vdata(i_start_index:i_end_index,1));
fpi_vydata=(fpi_i_vdata(i_start_index:i_end_index,2));
fpi_vzdata=(fpi_i_vdata(i_start_index:i_end_index,3));

vplot = subplot(num_plots,2,plot_order);
plot(fpi_i_timedata,[fpi_vxdata,fpi_vydata,fpi_vzdata],'LineWidth',1)
legend({'V_x','V_y','V_z'},'FontSize',10)
legend('boxoff')
legend('Location','eastoutside','AutoUpdate','off')
colormap('winter');
datetick
xlim([fpi_i_timedata(1) fpi_i_timedata(end)])
ylabel({'v_{ion}';'[km/s]'},'FontSize', 14)
set(gca, 'XTickLabel', [],'XMinorTick','on','YMinorTick','on','linewidth',1.25)


set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order/2-1+0.5), plot_pos(3), plot_pos(4)]);
plot_order = plot_order+2;
%% Temp
fpi_tparadata=(fpi_i_tparadata(i_start_index:i_end_index));
fpi_tperpdata=(fpi_i_tperpdata(i_start_index:i_end_index));

Tplot = subplot(num_plots,2,plot_order);
plot(fpi_i_timedata,[fpi_tparadata,fpi_tperpdata])
legend({'T^{i}_{\mid\mid}','T^{i}_{\perp}'},'FontSize',8)
legend('boxoff')
legend('Location','eastoutside','AutoUpdate','off')
colormap('winter');
datetick
xlim([fpi_i_timedata(1) fpi_i_timedata(end)])
ylabel({'Temp';'[eV]'},'FontSize', 14)
set(gca, 'YScale','log','XTickLabel', [],'XMinorTick','on','YMinorTick','on','Ytick',[0 10 100 1000],'linewidth',1.25)

set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order/2-1+0.5), plot_pos(3), plot_pos(4)]);
plot_order = plot_order+2;
%% Ion Energy
fpi_edata=fpi_i_edata(i_start_index:i_end_index,:);
fpi_espectdata=log10(fpi_i_espectdata(i_start_index:i_end_index,:));

Eplot = subplot(num_plots,2,plot_order);
imagesc(fpi_i_timedata,log10(fpi_edata(1,:)'),fpi_espectdata')
color_bar = colorbar('Ticks', [3, 4, 5, 6, 7,8],...
    'TickLabels', {'10^3', '10^4', '10^5', '10^6', '10^7','10^8'},'FontSize', 10);
ylabel(color_bar,{'keV/cm^2 s sr keV'},'FontSize', 8)
shading interp
whitejet = [1 1 1; jet];
colormap(whitejet);
ylabel({'Energy';'[eV]'},'FontSize', 14)
set(gca,'Ydir','normal', 'XTickLabel', [],'YMinorTick','on','XMinorTick','on','Yticklabels',[10 100 1000 10000],'layer','top','linewidth',1.25)
datetick
xlim([fpi_i_timedata(1) fpi_i_timedata(end)])

set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order/2-1+0.5), plot_pos(3), plot_pos(4)]);
plot_order = plot_order+2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

annotation('textbox',[plot_pos(1)-plot_pos(4)/2-0.005, plot_pos(2)-plot_pos(4)*plot_gap*(plot_order/2-1+0.5), plot_pos(3), plot_pos(4)],...
    'String',{datestr(event_start,'YYYY mmm dd'),'X-GSE (Re):', 'Y-GSE (Re):', 'Z-GSE (Re):'},...
    'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);

r_cols = length(mms1_mec_rdata_raw);
n_cs = tdnormal(event_start,event_end,mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,'event');
annotation('textbox',[plot_pos(1)+plot_pos(3)/r_cols+0.05, plot_pos(2)-plot_pos(4)*plot_gap*(plot_order/2-1+0.5), plot_pos(3), plot_pos(4)],...
    'String',{'',num2str(round(mms1_mec_rdata_raw(1,:)/6371.2*10)/10,'%2g\n')},...
    'VerticalAlignment','Top','Edgecolor','none','FontSize', 14);


slices_directory = strcat('~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis/Velocity Distributions/' ,num2str(EventNumber));
slices_directory = strcat('~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis/Velocity Distributions/EE1-MMS' ,num2str(probe));

cd(slices_directory)
directoryProperties = dir(pwd);
all_filenames = {directoryProperties.name};
all_filenames = all_filenames(contains(all_filenames,Rotation));
numberOfSlices = length(all_filenames);
sliceplot = subplot(num_plots,2,2:2:num_plots+2);
% set(gca,'Position',[plot_pos(1)-0.23, plot_pos(2)-plot_pos(4)*plot_gap*((plot_order-2)/2-1+0.75), 5.5*plot_pos(3), 5.5*plot_pos(4)]);
formatIn = 'yyyymmdd_HHMMSS.FFF';

vidfile = VideoWriter(strcat(num2str(EventNumber),'_VelocityDistribution','.mp4'),'MPEG-4');
vidfile.FrameRate = 4;
open(vidfile);

for i=1:numberOfSlices
    filename = char(all_filenames(i));
    timeMark = datenum(filename(1:19),formatIn);
    axes(Bplot)
    Bline = line([timeMark,timeMark],get(gca,'YLim'),'Color','r','LineWidth',2,'LineStyle','--');
    axes(nplot)
    nline = line([timeMark,timeMark],get(gca,'YLim'),'Color','r','LineWidth',2,'LineStyle','--');
    axes(vplot)
    vline = line([timeMark,timeMark],get(gca,'YLim'),'Color','r','LineWidth',2,'LineStyle','--');
    axes(Tplot)
    Tline = line([timeMark,timeMark],get(gca,'YLim'),'Color','r','LineWidth',2,'LineStyle','--');
    axes(Eplot)
    Eline = line([timeMark,timeMark],get(gca,'YLim'),'Color','r','LineWidth',2,'LineStyle','--');
    axes(sliceplot)
    imshow(filename)
    drawnow
    
    F(i) = getframe(gcf); 
    writeVideo(vidfile,F(i));

    delete(Bline);delete(nline);delete(vline);delete(Tline);delete(Eline);
    
end
close(vidfile)