%Timing Method for 4-spacecraft data
%date_start = '2018-03-01 01:05:40.000';
%date_end = '2018-03-01 01:06:22.000';
clear
close all
figure('Position',[0 0 800 600])
co = [0 0 1;
    0 1 0;
    1 0 0;
    0 0 0];
set(gcf,'defaultAxesColorOrder',co)
set(gcf,'color','w');
plot_order = 1;
plot_gap = 1.4;
num_plots = 3;
data_type='brst';

%Dates for Event Boundaries
date_start = '2018-01-09 08:34:28.000';
date_end = '2018-01-09 08:34:59.000';


date_start = '2018-01-09 08:34:49.250';
date_end = '2018-01-09 08:34:51.000';

event_start = '2018-03-01 01:03:45.000';
event_end = '2018-03-01 01:04:20.000';
date_start = '2018-03-01 01:04:11.000';
date_end = '2018-03-01 01:04:11.500';



%%
%%%%%%%%%%%%%%%%%%%%%%%%%Load % Interpolate % Crop Data%%%%%%%%%%%%%%%%%%%%
%load mec data
[mms1_mec_timedata_raw,mms1_mec_r_gsedata_raw] = load_mec(event_start,'1',data_type);
[mms2_mec_timedata_raw,mms2_mec_r_gsedata_raw] = load_mec(event_start,'2',data_type);
[mms3_mec_timedata_raw,mms3_mec_r_gsedata_raw] = load_mec(event_start,'3',data_type);
[mms4_mec_timedata_raw,mms4_mec_r_gsedata_raw] = load_mec(event_start,'4',data_type);

% %load fgm data
% [mms1_fpi_timedata_raw, mms1_fpi_ndata_raw,~,~]= load_fgm(event_start,'1',data_type);
% [mms2_fpi_timedata_raw, mms2_fpi_ndata_raw,~,~]= load_fgm(event_start,'2',data_type);
% [mms3_fpi_timedata_raw, mms3_fpi_ndata_raw,~,~]= load_fgm(event_start,'3',data_type);
% [mms4_fpi_timedata_raw, mms4_fpi_ndata_raw,~,~]= load_fgm(event_start,'4',data_type);


%load fpi data
[mms1_fpi_timedata_raw, mms1_fpi_ndata_raw,~,~]= load_fpi(event_start,'1',data_type,'e');
[mms2_fpi_timedata_raw, mms2_fpi_ndata_raw,~,~]= load_fpi(event_start,'2',data_type,'e');
[mms3_fpi_timedata_raw, mms3_fpi_ndata_raw,~,~]= load_fpi(event_start,'3',data_type,'e');
[mms4_fpi_timedata_raw, mms4_fpi_ndata_raw,~,~]= load_fpi(event_start,'4',data_type,'e');

%interpolate with mms1_fpi_time
[~,mms1_mec_r_gsedata_interp] = interpxyz(mms1_mec_timedata_raw,mms1_mec_r_gsedata_raw,mms1_fpi_timedata_raw);
[~,mms2_fpi_ndata_interp] = interpxyz(mms2_fpi_timedata_raw,mms2_fpi_ndata_raw,mms1_fpi_timedata_raw);
[~,mms2_mec_r_gsedata_interp] = interpxyz(mms2_mec_timedata_raw,mms2_mec_r_gsedata_raw,mms1_fpi_timedata_raw);
[~,mms3_fpi_ndata_interp] = interpxyz(mms3_fpi_timedata_raw,mms3_fpi_ndata_raw,mms1_fpi_timedata_raw);
[~,mms3_mec_r_gsedata_interp] = interpxyz(mms3_mec_timedata_raw,mms3_mec_r_gsedata_raw,mms1_fpi_timedata_raw);
[~,mms4_fpi_ndata_interp] = interpxyz(mms4_fpi_timedata_raw,mms4_fpi_ndata_raw,mms1_fpi_timedata_raw);
[~,mms4_mec_r_gsedata_interp] = interpxyz(mms4_mec_timedata_raw,mms4_mec_r_gsedata_raw,mms1_fpi_timedata_raw);


%Crop mec data to specific time period
[~,mms1_mec_r_gsedata,~,~] = crop(mms1_fpi_timedata_raw,mms1_mec_r_gsedata_interp,date_start,date_end);
[~,mms2_mec_r_gsedata,~,~] = crop(mms1_fpi_timedata_raw,mms2_mec_r_gsedata_interp,date_start,date_end);
[~,mms3_mec_r_gsedata,~,~] = crop(mms1_fpi_timedata_raw,mms3_mec_r_gsedata_interp,date_start,date_end);
[~,mms4_mec_r_gsedata,~,~] = crop(mms1_fpi_timedata_raw,mms4_mec_r_gsedata_interp,date_start,date_end);

%Crop fpi data to specific time period
[mms1_fpi_timedata,mms1_fpi_ndata,~,~] = crop(mms1_fpi_timedata_raw,mms1_fpi_ndata_raw,date_start,date_end);
[~,mms2_fpi_ndata,~,~] = crop(mms1_fpi_timedata_raw,mms2_fpi_ndata_interp,date_start,date_end);
[~,mms3_fpi_ndata,~,~] = crop(mms1_fpi_timedata_raw,mms3_fpi_ndata_interp,date_start,date_end);
[~,mms4_fpi_ndata,~,~] = crop(mms1_fpi_timedata_raw,mms4_fpi_ndata_interp,date_start,date_end);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%Correlation of Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Separate the data into quadrants,
%we really only test the correlation in the 2 & 3rd quadriles of mms1 bdata
all_data = length(mms1_fpi_timedata);
half_data = floor(all_data/2); %index for half data
quarter_data = floor(half_data/2); %index for quarter data

%Initialize correlation vectors, xyz, for each spacecraft pairings
%We only correlate half the data's worth of points
cc12 = zeros(half_data,1);
cc13 = zeros(half_data,1);
cc14 = zeros(half_data,1);


%loop over i=0 to data_points/2 for bxyz for mms234
for i=1:half_data %loop for each subsequent data point.
    %timing window is 1/4 - 3/4 of time range of mms1_btime
    
    %Correlation for this point,i, for this component,j, for MMS1 and MMS2
    temp_matrix = corrcoef(mms1_fpi_ndata(quarter_data:half_data+quarter_data),...
        mms2_fpi_ndata(i:half_data+i));
    cc12(i) =  temp_matrix(1,2);
    
    %Correlation for this point,i, for this component,j, for MMS1 and MMS3
    temp_matrix = corrcoef(mms1_fpi_ndata(quarter_data:half_data+quarter_data),...
        mms3_fpi_ndata(i:half_data+i));
    cc13(i) = temp_matrix(1,2);
    
    %Correlation for this point,i, for this component,j, for MMS1 and MMS4
    temp_matrix = corrcoef(mms1_fpi_ndata(quarter_data:half_data+quarter_data),...
        mms4_fpi_ndata(i:half_data+i));
    cc14(i) = temp_matrix(1,2);
    
    
end

%find the best index for the starting point of the time range for the most
%correlation, highest correlation coefficient "most same" after sliding.
%for each component and each spacecraft, find the largest correlation
%coefficient, closer to 1 is better. we can only choose the normal
%components from one direction, x y or z, so should choose the direction that has the
%mean highest correlation coefficient for all spacecraft pairs.


%find the max index and value for each spacecraft pairing correlation
%coefficients
[cc12max,cc12max_index] = max(cc12);
[cc13max,cc13max_index] = max(cc13);
[cc14max,cc14max_index] = max(cc14);

%the index is the starting point and the ending point is i+half_data, this
%max index has the largest correlation coefficient, thus we should use this
%for the timing method


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%Calculation through Timing Method%%%%%%%%%%%%%%%%%%%%%%%
%Initialize matrices for all three components' parameters
%time_boundary = zeros(4,3);
n_boundary = zeros(3,3);
v_boundary = zeros(1,3);



    %subscript is the starting point for btime and analysis, chosen because
    %of its highest correlation coefficient.
    index1=quarter_data;
    index2=cc12max_index;
    index3=cc13max_index;
    index4=cc14max_index;
    
    %time of crossing of the discontinity by each spacecraft
    time1=mms1_fpi_timedata(index1);
    time2=mms1_fpi_timedata(index2);
    time3=mms1_fpi_timedata(index3);
    time4=mms1_fpi_timedata(index4);
    
    
    %calculation of normal and speed
    %Calculate the time difference between spacecrafts
    time_vector = 86400*[time1-time2;
        time1-time3;
        time1-time4;
        time2-time3;
        time2-time4;
        time3-time4]'; %86400 is from datenum to seconds, datenum is in days.
    
    %Distance vector of each spacecraft
    r1=mms1_mec_r_gsedata(index1,1:3);
    r2=mms2_mec_r_gsedata(index2,1:3);
    r3=mms3_mec_r_gsedata(index3,1:3);
    r4=mms4_mec_r_gsedata(index4,1:3);
    
    %Calculate the distance difference between spacecrafts
    r_matrix = [r1-r2;
        r1-r3;
        r1-r4;
        r2-r3;
        r2-r4;
        r3-r4];
    
    %Calculate volumetric matrix
    %rearrange positions for easier calculation of R_alphabeta
    rx = [r1(1),r2(1),r3(1),r4(1)];
    ry = [r1(2),r2(2),r3(2),r4(2)];
    rz = [r1(3),r2(3),r3(3),r4(3)];
    
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
    R=zeros(3,3);
    %Set values of Volumetric tensor
    R(1,1) = Rxx;
    R(1,2) = Rxy;
    R(1,3) = Rxz;
    R(2,1) = Rxy;
    R(2,2) = Ryy;
    R(2,3) = Ryz;
    R(3,1) = Rxz;
    R(3,2) = Ryz;
    R(3,3) = Rzz;
    
    
    m_l= (1/4^2)*(time_vector*r_matrix)/R; %Section 12.1.2 Analysis Methods of Multispacecrafts
    mag_m_l = norm(m_l);
    
    v = 1/mag_m_l; %speed
    n = m_l*v; %normal
    
    %Save variables, columns are for each component
    time_boundary = [datetime(time1,'ConvertFrom','datenum','format','yyyy-MM-dd/HH:mm:ss.SSSS')
        datetime(time2,'ConvertFrom','datenum','format','yyyy-MM-dd/HH:mm:ss.SSSS')
        datetime(time3,'ConvertFrom','datenum','format','yyyy-MM-dd/HH:mm:ss.SSSS')
        datetime(time4,'ConvertFrom','datenum','format','yyyy-MM-dd/HH:mm:ss.SSSS')];
    
    n_boundary = n;
    v_boundary = v;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%Print Normal And Speed%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%From one CC plot, find the three highest cc for each spacecraft pairing,
%the closer to 1 these numbers are, the better.
%the component that has all three numbers closest to 1 is the best for the
%timing method and whose normal and speeds we should use.

%find the best component to use based on mean correlation coefficients for
%each spacecraft pairing
% [cc_max_component_value,cc_max_component_index] = ...
%     max([mean([cc12max(1),cc13max(1),cc14max(1)]);...
%     mean([cc12max(2),cc13max(2),cc14max(2)]);...
%     mean([cc12max(3),cc13max(3),cc14max(3)])]);

%display normal and speeds
n = sprintf(strcat(num2str(n_boundary','[%.4f %.4f %.4f]')))

v = sprintf(strcat(num2str(v_boundary','%3.1f'),'km/s'))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plotting%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Plot the total magnetic field for each spacecraft
plot_fgm_magnetic(date_start,date_end,mms1_fpi_timedata,mms1_fpi_ndata,num_plots,plot_order)
hold on
plot_fgm_magnetic(date_start,date_end,mms1_fpi_timedata,mms2_fpi_ndata,num_plots,plot_order)
plot_fgm_magnetic(date_start,date_end,mms1_fpi_timedata,mms3_fpi_ndata,num_plots,plot_order)
plot_fgm_magnetic(date_start,date_end,mms1_fpi_timedata,mms4_fpi_ndata,num_plots,plot_order)
hold off
ylabel({'N_e';'[cm^-3]'},'FontSize', 14)
legend('off')

title('MMS1 Timing Method: Correlation Coefficients', 'FontSize', 18, 'FontWeight', 'normal')
legend({'1', '2', '3','4'},'FontSize',14)
legend('boxoff')
legend('Location','eastoutside')
datetick
xlim([mms1_fpi_timedata(1) mms1_fpi_timedata(end)])
set(gca,'XMinorTick','on','YMinorTick','on','linewidth',1.25)
plot_pos = get(gca,'Position');
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);

plot_order = plot_order+1;

time_delay = 0.0078*linspace(-quarter_data,quarter_data,half_data);%7.8ms for brst data


%Plot the time delay
subplot(num_plots,1,plot_order)
plot(time_delay,[cc12,cc13,cc14])
legend({'1-2', '1-3', '1-4'},'FontSize',14)
legend('boxoff')
legend('Location','eastoutside')
set(gca,'XMinorTick','on','YMinorTick','on','linewidth',1.25)
ylabel({'cc'},'FontSize', 14)
xlabel({'Time Delay [s]'},'FontSize', 14)
set(gca,'Position',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1), plot_pos(3), plot_pos(4)]);

plot_order = plot_order+1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%Annotations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Date
date_range = {date_start;date_end};
annotation('textbox',[plot_pos(1), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1.2), plot_pos(3), plot_pos(4)],...
    'String',date_range,...
    'VerticalAlignment','Top','Edgecolor','none','FontSize', 12);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Anotations%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Mean cc
annotation('textbox',[plot_pos(1)+plot_pos(4), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1.1), plot_pos(3), plot_pos(4)],...
    'String',strcat('cc_{mean}: ', num2str(mean([cc12max,cc13max,cc14max]))),...
    'VerticalAlignment','Top','Edgecolor','none','FontSize', 12);
%Timings for each spacecraft
annotation('textbox',[plot_pos(1)+plot_pos(4), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1.1)-0.04, plot_pos(3), plot_pos(4)],...
    'String',strcat('Timings:', datestr(time_boundary(:,1),'HHMMSS.FFF')),...
    'VerticalAlignment','Top','Edgecolor','none','FontSize', 12);
%Normal
annotation('textbox',[plot_pos(1)+2*plot_pos(4), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1.1)-0.04, plot_pos(3), plot_pos(4)],...
    'String',strcat('Normal: ',' ', num2str(n_boundary,'[%0.4f, %0.4f, %0.4f]')),...
    'VerticalAlignment','Top','Edgecolor','none','FontSize', 12);
%speed
annotation('textbox',[plot_pos(1)+2*plot_pos(4), plot_pos(2)-plot_pos(4)*plot_gap*(plot_order-1.1)-0.10, plot_pos(3), plot_pos(4)],...
    'String',strcat('Speed: ', num2str(v_boundary,'%3.2f'),'km/s'),...
    'VerticalAlignment','Top','Edgecolor','none','FontSize', 12);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%Save Plots
orient(gcf,'landscape')
plot_name =  strcat('mms','_Timing_Ne_',date_start(1:19),'.pdf');
print(gcf, '-dpdf', '-opengl',plot_name,'-fillpage');
movefile(plot_name, '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis')
