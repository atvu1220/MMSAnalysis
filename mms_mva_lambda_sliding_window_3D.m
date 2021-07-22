%Andrew Vu 9/25/18
%data=spdfcdfread(filename);
%datainfo=spdfcdfinfo(filename);
%clear
tic
close all
clear
figure('Position',[1 1 1100 700])
mms_directory = '/Users/andrewvu/data/mms/';
probe_num = '1';
num_plots = 3;
addpath ~/Library/'Mobile Documents'/com~apple~CloudDocs/
addpath ~/data/

data_type = 'srvy';
window_lower_limit = 16*10;
data_points = 16*30; %total on both sides.
minRatio = 10; %mid/low eigenvalue ratio

date_start = '2018-03-12 07:35:50.000';
date_end = '2018-03-12 07:36:30.000';

% date_start = '2018-01-09 08:34:27.000';
% date_end = '2018-01-09 08:35:03.000';
%
% date_start = '2018-04-01 01:10:47.000';
% date_end = '2018-04-01 01:11:05.000';
%
% date_start = '2018-04-27 19:47:30.000';
% date_end = '2018-04-27 19:48:15.000';
%
% date_start = '2018-03-01 01:05:40.000';
% date_end = '2018-03-01 01:06:22.000';

date_start = '2018-01-09 08:34:27.000';
date_end = '2018-01-09 08:35:03.000';

date_start = '2018-03-01 01:03:44.000';
date_end = '2018-03-01 01:06:22.000';

date_start = '2018-01-29 03:36:10.000';
date_end = '2018-01-29 03:36:20.000';

date_start = strcat('2017-11-04 05:15:45.0');
date_end = strcat('2017-11-04 05:17:45.0');

%0 for 10 plus, 1 for max, 3 is for min
flag = 0;
if flag == 0
    plot_type = '10plus';
    colorbarString = strcat('log_{10}(\lambda>',num2str(minRatio),')');
elseif flag == 1
    plot_type = 'N_{max}';
    colorbarString = plot_type;
elseif flag == 3
    plot_type = 'N_{min}';
    colorbarString = plot_type;
elseif flag == 4
    plot_type = 'angle';
    colorbarString = 'Angle with N_{cs}';
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Data Retrieval%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[fgm_timedata, fgm_bdata] = load_fgm(date_start,date_end,1,data_type);

formatIn='yyyy-mm-dd HH:MM:SS.FFF';
tstart = datenum(date_start,formatIn);
tend = datenum(date_end,formatIn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%Magnetic Field Data%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find the start and end limits of the event in the data
start_index = find(fgm_timedata > tstart, 1)-1;
end_index = find(fgm_timedata > tend, 1);


f1 = subplot(num_plots,1,1);
set(gcf,'color','w');
co = [0 0 1;
    0 1 0;
    1 0 0;
    0 0 0];
set(gcf,'defaultAxesColorOrder',co)
hold on
plot(fgm_timedata(start_index:end_index,:),fgm_bdata(start_index:end_index,1:4),'LineWidth',2)
set(gca,'linewidth',1.25)
%Display borders for largest window size on figure
% mid = start_index+round(length(fgm_timedata(start_index:end_index))/2);
% line([(fgm_timedata(mid-data_points/2)) (fgm_timedata(mid-data_points/2))],...
%     get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
% line([(fgm_timedata(mid+data_points/2)) (fgm_timedata(mid+data_points/2))],...
%     get(gca,'YLim'),'Color','k','LineWidth',0.5,'LineStyle','--')
box on
hold off

datetick
set(gca,'layer','top')
legend({'B_x', 'B_y', 'B_z','B_{tot}'},'FontSize',14)
legend('boxoff');
legend('Location','eastoutside');
%xlim([datetime(tstart,'ConvertFrom','datenum') datetime(tend,'ConvertFrom','datenum')]);
ylabel({'B';'[nT]'},'FontSize', 14)

title_name = strcat('MMS1 MVAB Eigenvalue Ratios: \Deltat=', num2str(data_points/16),'s');
title(title_name, 'FontSize', 18, 'FontWeight', 'normal')

xlim([fgm_timedata(start_index) fgm_timedata(end_index)])
axis tight

set(f1,'XMinorTick','on','YMinorTick','on')
plot_pos = get(f1,'Position');
set(f1,'Position',[plot_pos(1), plot_pos(2), plot_pos(3), plot_pos(4)]);

l1= findobj(gcf, 'Type', 'Legend');
l1pos=get(l1, 'Position');
set(l1,'Position',[plot_pos(1)+plot_pos(3)+0.025, plot_pos(2)+plot_pos(4)/2-l1pos(4)/2, l1pos(3), l1pos(4)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


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

nested_n_min_x=zeros(window_range,ceil(window_size/2));
nested_n_min_y=zeros(window_range,ceil(window_size/2));
nested_n_min_z=zeros(window_range,ceil(window_size/2));

nested_n_max_x=zeros(window_range,ceil(window_size/2));
nested_n_max_y=zeros(window_range,ceil(window_size/2));
nested_n_max_z=zeros(window_range,ceil(window_size/2));

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
        
    end
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%lambda ratios %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f2 = subplot(num_plots,1,2);
%just in case for negatives
% nested_lmidmin(nested_lmidmin < 0) = NaN;


%%%%%%%%%%%for normal case
% flag = '';
% t_0_axis = linspace(1,length(t_0_vector),length(t_0_vector));
% pcolor(t_0_vector,2*data_point_vector,log(nested_lmidmin'))


%%%%%%%%%%for just >10 case
if flag == 0
    nested_lmidmin(nested_lmidmin < minRatio ) = NaN;
    % pcolor(t_0_vector,2*data_point_vector,log10(nested_lmidmin'))
    imagesc(fgm_timedata(start_index:end_index,:),2*data_point_vector/16,log10(nested_lmidmin'))
    ylabel(strcat({'\lambda >'},num2str(minRatio)),'FontSize',14)
    
    %%%%%%%%%for nmin comparsions case
elseif  flag == 3
    nested_lmidmin(nested_lmidmin < minRatio ) = NaN;
    nested_n_min_x_10plus = ~isnan(nested_lmidmin).*nested_n_min_x;
    nested_n_min_x_10plus(nested_n_min_x_10plus == 0) = NaN;
    imagesc(fgm_timedata(start_index:end_index,:),2*data_point_vector/16,abs(nested_n_min_x_10plus'))
    c=colorbar;
    ylabel(c,{strcat('N_i for \lambda >',num2str(minRatio))},'FontSize',14)
    
    
    %%%%%%%%%for nmin comparsions case
elseif  flag == 1
    
    nested_lmaxmid(nested_lmaxmid < minRatio ) = NaN;
    nested_n_max_x_10plus = ~isnan(nested_lmaxmid).*nested_n_max_x;
    nested_n_max_x_10plus(nested_n_max_x_10plus == 0) = NaN;
    imagesc(fgm_timedata(start_index:end_index,:),2*data_point_vector/16,abs(nested_n_max_x_10plus'))
    c=colorbar;
    ylabel(c,{strcat('N_i for \lambda >',num2str(minRatio))},'FontSize',14)
    
elseif  flag == 4
    
    nested_lmidmin(nested_lmidmin < minRatio ) = NaN;
    nested_n_min_x_10plus = ~isnan(nested_lmidmin).*nested_n_min_x;
    nested_n_min_x_10plus(nested_n_min_x_10plus == 0) = NaN;
    nested_n_min_y_10plus = ~isnan(nested_lmidmin).*nested_n_min_x;
    nested_n_min_y_10plus(nested_n_min_y_10plus == 0) = NaN;
    nested_n_min_z_10plus = ~isnan(nested_lmidmin).*nested_n_min_x;
    nested_n_min_z_10plus(nested_n_min_z_10plus == 0) = NaN;
    nested_angles = zeros(window_range,ceil(window_size/2));
    
    n_td = tdnormal(date_start,date_end,fgm_timedata,fgm_bdata,'all');
    for i=1:window_range
        for j=window_lower_limit:ceil(window_size/2)
            point_normal = [nested_n_min_x_10plus(i,j), nested_n_min_y_10plus(i,j), nested_n_min_z_10plus(i,j)];
            %nested_angles(i,j) = rad2deg(acos(dot(point_normal,n_cs)/(norm(point_normal)*norm(n_cs))));
            if isnan(point_normal) == [1,1,1]
                nested_angles(i,j) = NaN;
            else
                [angle,~] = ambig_angle(point_normal,n_td);
                nested_angles(i,j) = angle;
            end
            
        end
    end
    %     for i=1:window_range
    %         for j=window_lower_limit:ceil(window_size/2)
    %             point_normal = [nested_n_min_x_10plus(i,j), nested_n_min_y_10plus(i,j), nested_n_min_z_10plus(i,j)];
    %             nested_angles(i,j) = rad2deg(acos(dot(point_normal,n_td)/(norm(point_normal)*norm(n_td))));
    %         end
    %     end
    
    imagesc(fgm_timedata(start_index:end_index,:),2*data_point_vector/16,abs(nested_angles'))
    c=colorbar;
    ylabel(c,{strcat('Angle with n_{td} for \lambda >',num2str(minRatio))},'FontSize',14)
end




set(gca,'linewidth',1.25);
c=colorbar;
ylabel(c,colorbarString,'FontSize',14)
set(gca,'Ydir','reverse')
set(gca,'layer','top')
set(gcf,'color','w');
box on
% c=colorbar;
cmap = jet;
cmap = [1 1 1; cmap];
colormap(cmap)
%set(gca,'colorscale','log')
xlim([fgm_timedata(start_index) fgm_timedata(end_index)])
datetick('x','keeplimits')
set(gca,'TickDir','out');
xtickangle(-40)
%ylabel(c,{'\lambda_{mid}/\lambda_{min} > 10'},'FontSize', 14)

%set(c,'TickLabels',[5 10 15 20 25])
ylabel({'Window Size';'[Seconds]'},'FontSize', 14)
xlabel(strcat('t_0 position for \lambda > ',num2str(minRatio)))

set(f2,'XTickLabel',[],'XMinorTick','on','YMinorTick','on','Position',...
    [plot_pos(1), plot_pos(2)-1*plot_pos(4)-0.05, plot_pos(3), 1*plot_pos(4)]);

l2= findobj(gcf, 'Type', 'ColorBar');




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% annotation('textbox',[plot_pos(1), plot_pos(2)-2*plot_pos(4)-.1, plot_pos(3), plot_pos(4)*1],...
%     'String',{date_start(1:10)},...
%     'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);
orient(gcf,'landscape')
plot_name =  strcat('mms',probe_num,'_mvab_slidingwindow3D_',flag,...
    date_start(1:19),'_',num2str(data_points),'.pdf');

print(gcf, '-dpdf', '-opengl', plot_name,'-fillpage');

movefile(plot_name, '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis')

toc