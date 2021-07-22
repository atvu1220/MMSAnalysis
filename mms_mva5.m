%Andrew Vu 9/25/18
%data=spdfcdfread(filename);
%datainfo=spdfcdfinfo(filename);
%clear
close all
clear
fig = figure('Position',[1 1 750 1200])
%figure('PaperPositionMode', 'auto');
mms_directory = '/Users/andrewvu/data/mms/';
probe_num = '1';
num_plots = 4;
data_type = 'srvy';
window_lower_limit = 16*10;
data_points = 16*30; %total on both sides.
minRatio = 10; %mid/low eigenvalue ratio
% date = '2018/03/18';
% date = '2018/03/12';
% %date=['2018-03-12/07:35:50']
%
% date_start = '2018-03-18 00:10:15.000';
% date_end = '2018-03-18 00:11:15.000';
%
% date_start = '2018-03-18 00:10:36.000';
% date_end = '2018-03-18 00:10:40.000';
%
% date_start = '2018-03-18 00:10:43.000';
% date_end = '2018-03-18 00:10:44.500';
%
% date_start = '2018-03-18 00:10:44.500';
% date_end = '2018-03-18 00:10:46.500';
%
% date_start = '2018-03-18 00:10:46.500';
% date_end = '2018-03-18 00:10:49.000';
%
% date_start = '2018-03-18 00:10:36.000';
% date_end = '2018-03-18 00:10:40.000';
%
%
% month = 03;
% day = 12;
% hour = 7;
% min = 36;


% date_start = '2018-03-12 07:35:50.000';
% date_start_vec = datevec(date_start);
% year = date_start_vec(1);
% month = date_start_vec(2);
% day = date_start_vec(3);
% hour = date_start_vec(4);
% min_start = date_start_vec(5);
% sec_start = date_start_vec(6);
%
% date_end = '2018-03-12 07:36:30.000';
% date_end_vec = datevec(date_start);
% min_end = date_end_vec(5);
% sec_end = date_end_vec(6);

%seconds
% center = 17;
% left = center -2;
% right = center+2;
%
% date_start = strcat('2018-03-12 07:36:',num2str(left),'.0');
% date_end = strcat('2018-03-12 07:36:',num2str(right),'.0');

date_start = strcat('2018-03-12 07:35:59.0');
date_end = strcat('2018-03-12 07:36:09.0');

date_start = strcat('2017-11-04 05:15:45.0');
date_end = strcat('2017-11-04 05:17:00.0');

duration=120;

%datetime(2018,03,18,00,[10 11],15,000)
formatIn='yyyy-mm-dd HH:MM:SS.FFF';
tstart = datenum(date_start,formatIn);
tend = datenum(date_end,formatIn);


%%%%%%%%%%%%%%%%%%%%%%%%%%Magnetic Field Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fgm_timedata, fgm_bdata] = load_fgm(date_start,date_end,1,data_type);
fgm_bdata(:,1) = smooth(fgm_bdata(:,1));
fgm_bdata(:,2) = smooth(fgm_bdata(:,2));
fgm_bdata(:,3) = smooth(fgm_bdata(:,3));
fgm_bdata(:,4) = smooth(fgm_bdata(:,4));
%Find the start and end limits of the event in the data
start_index = find(fgm_timedata > tstart, 1)-1;
end_index = find(fgm_timedata > tend, 1);
% fgm_timedata = fgm_timedata(start_index:end_index,1);
% bdata = bdata(start_index:end_index,:);




%%%%%%%%%%%%%%%%%%%%%%%Loop over sliding window%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%arbitrary scope set dt_0
window_range = end_index - start_index+1;
%srvy has 16 points per second, so 64 points is 4 seconds total
%data poitns is total, so each side is data_points/2
window_size = data_points;

nested_lambda=cell(window_range,ceil(window_size/2));
nested_normal=cell(window_range,ceil(window_size/2));
nested_time=cell(window_range,ceil(window_size/2));
nested_lmidmin=zeros(window_range,ceil(window_size/2));
nested_lmaxmid=zeros(window_range,ceil(window_size/2));

nested_n_min_x=zeros(window_range,ceil(window_size/2));
nested_n_min_y=zeros(window_range,ceil(window_size/2));
nested_n_min_z=zeros(window_range,ceil(window_size/2));

nested_n_max_x=zeros(window_range,ceil(window_size/2));
nested_n_max_y=zeros(window_range,ceil(window_size/2));
nested_n_max_z=zeros(window_range,ceil(window_size/2));
nested_matrix=cell(window_range,ceil(window_size/2));

t_0_vector = zeros(window_range,1);
data_point_vector = zeros(1,ceil(window_size/2),1);

for i=1:window_range
    %i is the starting point index
    t_0 = fgm_timedata(i+start_index-1);
    t_0_vector(i) = t_0;
    %window_lower_limit is how many of the first points to skip
    for j=window_lower_limit:ceil(window_size/2)
        data_point_vector(j) = j;
        %j is the number of data points to one side of the center point.
        left_index = (i+start_index)-j;
        right_index = (i+start_index)+j;
        
        bdata_scope = fgm_bdata(left_index:right_index,:);
        
        [output, l, v] = mvab(bdata_scope(:,1:3));
        
        nested_lmidmin(i,j)= l(2)/l(3);
        nested_lmaxmid(i,j) = l(1)/l(3);
        
        nested_n_min_x(i,j)=v(1,3);
        nested_n_min_y(i,j)=v(2,3);
        nested_n_min_z(i,j)=v(3,3);
        
        nested_n_max_x(i,j)=v(1,1);
        nested_n_max_y(i,j)=v(2,1);
        nested_n_max_z(i,j)=v(3,1);
        
        nested_lambda{i,j} = l;
        nested_normal{i,j} = v(:,3);
        nested_time{i,j}=t_0;
        nested_matrix{i,j} = v;
        
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[~,indexOfMaxEigenvalueRatio] = max(nested_lmidmin.*nested_lmaxmid,[],'all','linear');
output = fgm_bdata(:,1:3)*cell2mat(nested_matrix(indexOfMaxEigenvalueRatio));
lambda = cell2mat(nested_lambda(indexOfMaxEigenvalueRatio));
normal = cell2mat(nested_normal(indexOfMaxEigenvalueRatio))
% [output, lambda, v] = mvab(bdata(:,1:3));

v = VideoWriter('wave_analysis_mva','Motion JPEG AVI');
v.FrameRate= 16;
open(v);

%%%%%%%%%%%%%%%%%%%%%%%Plot Original Magnetic Field%%%%%%%%%%%%%%%%%%%%%%%
for i=1:8:end_index-start_index
    
    f1 = subplot(num_plots,2,[1 2]);
    set(gcf,'color','w');
    co = [1 0 0;
        0 1 0;
        0 0 1;
        0 0 0];
    set(gcf,'defaultAxesColorOrder',co)
    plot(fgm_timedata(start_index:start_index+i),fgm_bdata(start_index:start_index+i,1:4))
    set(gca,'XMinorTick','off','TickDir','out','YMinorTick','on','linewidth',2)
    legend({'B_x', 'B_y', 'B_z','B_{tot}'},'FontSize',14)
    legend('boxoff')
    legend('Location','eastoutside')
    colormap('jet');
    % xlim([fgm_timedata(1) fgm_timedata(end)])
    xlim([fgm_timedata(start_index) fgm_timedata(end_index)])
    datetick
    ylabel({'B';'[nT]'},'FontSize', 14)
    set(f1, 'XTickLabel', [],'XMinorTick','on','YMinorTick','on')
    title('MMS1 Minimum Variance Analysis on B', 'FontSize', 18, 'FontWeight', 'normal')
    
    plot_pos = get(f1,'Position');
    set(f1,'Position',[plot_pos(1), plot_pos(2), plot_pos(3), plot_pos(4)]);
    
    
    % l1= findobj(gcf, 'Type', 'Legend');
    % l1pos=get(l1, 'Position');
    % set(l1,'Position',[plot_pos(1)+plot_pos(3)+0.025,...
    %     plot_pos(2)+plot_pos(4)/2-l1pos(4)/2, l1pos(3), l1pos(4)]);
    
    % hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%Transformed Magnetic Field%%%%%%%%%%%%%%%%%%%%%%%%
    %plot transformed magnetic field
    f2 = subplot(num_plots,2,[3 4]);
    
    set(gcf,'color','w');
    plot(fgm_timedata(start_index:start_index+i),output(start_index:start_index+i,:))
    set(gca,'XMinorTick','off','TickDir','out','YMinorTick','on','linewidth',2)
    legend({'B_i', 'B_j', 'B_k'},'FontSize',14)
    legend('boxoff','Location','eastoutside')
    colormap('jet');
    
    xlim([fgm_timedata(start_index) fgm_timedata(end_index)])
    datetick
    ylabel({'B_{Transformed}';'[nT]'},'FontSize', 14)
    set(f2,'XMinorTick','on','YMinorTick','on','Position',[plot_pos(1), plot_pos(2)-plot_pos(4), plot_pos(3), plot_pos(4)]);
    
    l2= findobj(gcf, 'Type', 'Legend');
    l2pos=get(l2(2), 'Position');
    set(l2(1),'Position',[plot_pos(1)+plot_pos(3)+0.025,...
        plot_pos(2)-plot_pos(4)+plot_pos(4)/2-l2pos(4)/2, l2pos(3), l2pos(4)]);
    % hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Hodograms%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot transformed magnetic field
    
    f3 = subplot(num_plots,2,[5 8]);
    axis([-1.5 1.5 -3 3])
    %set(gcf,'color','w','MarkerFaceColor','k');
    % plot(smooth(output(1,2)),smooth(output(1,1)),'color',[ 1 0 0]);
    plot(smooth(output(start_index:(2*16):start_index+i,2)),smooth(output(start_index:(2*16):start_index+i,1)),'color',[ 0 0 0]);
    set(gca,'XMinorTick','off','TickDir','out','YMinorTick','on','linewidth',2)
   % axis([-1.5 1.5 -3 3])
    ylabel({'B_i'},'FontSize', 14)
    xlabel({'B_j'},'FontSize', 14)
        title('Intermediate vs Maxmimum Eigenvector Hodogram', 'FontSize', 18, 'FontWeight', 'normal')

    % hold on
    %
    % f4 = subplot(num_plots,2,[6 8]);
    % %set(gcf,'color','w','MarkerFaceColor','k');
    % plot(smooth(output(start_index:start_index+i,3)),smooth(output(start_index:start_index+i,1)),'color',[ 0 0 0]);
    % xlabel({'B_k [nT]'},'FontSize', 14)
    % set(f4,'YTick',[]);
    %
    % set(f3,'Position',[plot_pos(1), plot_pos(2)-4*plot_pos(4)-0.025, plot_pos(3)/2, 3*plot_pos(4)]);
    % set(f4,'Position',[plot_pos(1)+plot_pos(3)/2, plot_pos(2)-4*plot_pos(4)-0.025, plot_pos(3)/2, 3*plot_pos(4)]);
    % hold on
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    drawnow
    g(i) = getframe(fig);
    writeVideo(v,g(i));
    
end
close(v)

% %%%%%%%%%%%%%%%%%%%%%%%%%%Normals and Eigenvalues%%%%%%%%%%%%%%%%%%%%%%%%%%
% annotation('textbox',[plot_pos(1)+plot_pos(3)+0.0125, plot_pos(2)-2.5* plot_pos(4), plot_pos(3), plot_pos(4)/2],...
%     'String',{'Boundary N:',v(1,3),v(2,3),v(3,3)},...
%     'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 16);
% annotation('textbox',[plot_pos(1)+plot_pos(3)+0.0125, plot_pos(2)-3.5*plot_pos(4), plot_pos(3), plot_pos(4)/2],...
%     'String',{'Eigenvalues:',strcat('tmax = ',num2str(lambda(1))),...
%     strcat('tmid = ',num2str(lambda(2))),strcat('tmin = ',num2str(lambda(3)))},...
%     'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 16);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% annotation('textbox',[plot_pos(1), plot_pos(2)-plot_pos(4), plot_pos(3), plot_pos(4)*1],...
%     'String',{date_start(1:10)},...
%     'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);
%
%
plot_name =  strcat('mms',probe_num,'_mvab_',date_start(1:19),'_',num2str(duration),'.pdf');
% %
% % print(gcf, '-dpdf', '-opengl',plot_name,'-fillpage');

movefile(plot_name, '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis')


