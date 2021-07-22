%Andrew Vu 9/25/18
%This program calculates MVAB multiple times from a FIXED center point, t_0,
%with increasing window size
tic
close all
clear
figure('Position',[1 1 1050 800])
mms_directory = '/Users/andrewvu/data/mms/';
probe_num = '1';
num_plots = 4;
data_type = 'srvy';



%center is t_0. left and right are the first (initial) bounds for t_0,
%program will continue outward from 'left' and 'right'.

%entire event times
date_start = '2018-03-12 07:35:50.000';
date_end = '2018-03-12 07:36:30.000';

%center of discontinuity
center_time = '2018-03-12 07:36:04.000';

%dt is window size, on each side.
dt = 1;
duration=2*dt;

%Create left and Right Initial Bounds
formatIn='yyyy-MM-dd HH:mm:ss.SSS';
center_datetime = datetime(center_time,'InputFormat',formatIn);
left = center_datetime - seconds(dt);
right = center_datetime + seconds(dt);

% %Create start date string for load_fgm
% window_start = datestr(left,'yyyy-mm-dd HH:MM:SS.FFF');
% window_end = datestr(right,'yyyy-mm-dd HH:MM:SS.FFF');




%load FGM data
[fgm_timedata, bdata] = load_fgm(date_start,1,data_type);




%%%%%%%%%%%%%%%%%%%%%%%%%%Plot Magnetic Field%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Find the start and end limits of the entire event in the data
event_start = find(fgm_timedata > datenum(date_start), 1)-1;
event_end = find(fgm_timedata > datenum(date_end), 1);

f1 = subplot(num_plots,1,1);
set(gcf,'color','w');
co = [0 0 1;
    0 1 0;
    1 0 0;
    0 0 0];
set(gcf,'defaultAxesColorOrder',co)
plot(fgm_timedata(event_start:event_end),bdata(event_start:event_end,1:4))
legend({'B_x', 'B_y', 'B_z','B_{tot}'},'FontSize',14)
legend('boxoff');
legend('Location','eastoutside');
datetick('x','HH:MM:ss')
xlim([fgm_timedata(event_start) fgm_timedata(event_end)])
ylabel({'B';'[nT]'},'FontSize', 14)
set(f1,'XMinorTick','on','YMinorTick','on')
title('MMS1 MVAB Nested Boundary Normals and Eigenvalues', 'FontSize', 18, 'FontWeight', 'normal')

plot_pos = get(f1,'Position');
set(f1,'Position',[plot_pos(1), plot_pos(2), plot_pos(3), plot_pos(4)]);

l1= findobj(gcf, 'Type', 'Legend');
l1pos=get(l1, 'Position');
set(l1,'Position',[plot_pos(1)+plot_pos(3)+0.025, plot_pos(2)+plot_pos(4)/2-l1pos(4)/2, l1pos(3), l1pos(4)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%Loop over nested scopes%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%For indexing timedata,
tstart = datenum(left);
tend = datenum(right);

%Find the indices of the left and right initial window (between left and
%right bounds)
left_index = find(fgm_timedata > tstart, 1)-1;
right_index = find(fgm_timedata > tend, 1)-1;

%arbitrary scope set dt_0, initial window range in # indices, minimum
%number of indices.
deltat_0 = right_index - left_index;

%center point index
t_0 = floor(left_index+right_index)/2;

%multiple of the minimal-window-sizes to perform loops over.
factor=2;

%Display only the B data where we actually enclose the increasing interval.
leftleft = fgm_timedata(t_0-factor*deltat_0);
rightright = fgm_timedata(t_0+factor*deltat_0);
xlim(([leftleft rightright]))
hold on
line(([(fgm_timedata(t_0)) (fgm_timedata(t_0))]),get(gca,'YLim'),'Color','k','LineWidth',2)
legend({'B_x', 'B_y', 'B_z','B_{tot}','t_0'},'FontSize',14)
hold off

%initialize loop matrices
nested_lambda=zeros(3,factor*deltat_0);
nested_normal=zeros(3,factor*deltat_0);
nested_time=zeros(1,factor*deltat_0);
left_index = t_0-1;
right_index =t_0+1;

for i=1:factor*deltat_0
    %plot transformed magnetic field
    bdata_scope = bdata(left_index:right_index,:);
    
    [output, l, v] = mvab(bdata_scope(:,1:3));
    nested_lambda(:,i) = l;
    nested_normal(:,i) = v(:,3);
    nested_time(1,i)=fgm_timedata(right_index,:)-fgm_timedata(t_0,:);
    
    left_index = t_0-i;
    right_index = t_0+i;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%normals %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f2 = subplot(num_plots,1,2);
set(gcf,'color','w');

plot(datetime(nested_time,'ConvertFrom','datenum','Format','ss'),nested_normal')
legend({'N_x', 'N_y', 'N_z'},'FontSize',14)
legend('boxoff')
legend('Location','eastoutside')
ylabel({'Normal';'Components'},'FontSize', 14)
set(f2,'XTickLabel', [], 'XMinorTick','on','YMinorTick','on','Position',...
    [plot_pos(1), plot_pos(2)-plot_pos(4)-0.025, plot_pos(3), plot_pos(4)]);
l2= findobj(gcf, 'Type', 'Legend');
l2pos=get(l2(2), 'Position');
set(l2(1),'Position',[plot_pos(1)+plot_pos(3)+0.025, ...
    plot_pos(2)-plot_pos(4)+plot_pos(4)/2-l2pos(4)/2-0.025, l2pos(3), l2pos(4)]);




%%%%%%%%%%%%%%%%%%%%%%%%%nested eigenvalues%%%%%%%%%%%%%%%%%%%%%%%%
%plot nested eigenvalues
f3 = subplot(num_plots,1,3);
set(gcf,'color','w');

plot(datetime(nested_time,'ConvertFrom','datenum','Format','ss'),nested_lambda')
legend({'\lambda_{max}', '\lambda_{mid}', '\lambda_{min}'},'FontSize',14)
legend('boxoff')
legend('Location','eastoutside')
ylabel({'Eigenvalues'},'FontSize', 14)
set(f3,'XTickLabel', [], 'XMinorTick','on','YMinorTick','on','Position',...
    [plot_pos(1), plot_pos(2)-2*plot_pos(4)-0.05, plot_pos(3), plot_pos(4)]);
l3= findobj(gcf, 'Type', 'Legend');
l3pos=get(l3(1), 'Position');
set(l3(1),'Position',[plot_pos(1)+plot_pos(3)+0.035, ...
    plot_pos(2)-2*plot_pos(4)+plot_pos(4)/2-l2pos(4)/2-0.05, l2pos(3), l2pos(4)]);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%eigenvalue ratios%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f4 = subplot(num_plots,1,4);
set(gcf,'color','w');
lam12=abs((nested_lambda(1,:)./nested_lambda(2,:)));
lam23=abs((nested_lambda(2,:)./nested_lambda(3,:)));
lam31=abs((nested_lambda(3,:)./nested_lambda(1,:)));

lam12(lam12>25) = NaN;
lam23(lam23>25) = NaN;
lam31(lam31>25) = NaN;

plot(datetime(nested_time,'ConvertFrom','datenum','Format',['ss']),[lam12',lam23',lam31'])
legend({'\lambda_{max}/\lambda_{mid}', '\lambda_{mid}/\lambda_{min}', '\lambda_{min}/\lambda_{max}'},'FontSize',14)
legend('boxoff')
legend('Location','eastoutside')
ylabel({'Eigenvalues Ratios'},'FontSize', 14)
xlabel('Seconds from t_0')
set(f4,'XMinorTick','on','YMinorTick','on','Position',...
    [plot_pos(1), plot_pos(2)-3*plot_pos(4)-0.075, plot_pos(3), plot_pos(4)]);
xtickformat('ss.SSS')
l4= findobj(gcf, 'Type', 'Legend');
l4pos=get(l4(1), 'Position');
set(l4(1),'Position',[plot_pos(1)+plot_pos(3)+0.065, ...
    plot_pos(2)-3*plot_pos(4)+plot_pos(4)/2-l2pos(4)/2-0.075, l2pos(3), l2pos(4)]);
annotation('textbox',[plot_pos(1), plot_pos(2)-3*plot_pos(4)-0.125, plot_pos(3), plot_pos(4)*1],...
    'String',{date_start(1:10)},...
    'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 14);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%Normals and Eigenvalues%%%%%%%%%%%%%%%%%%%%%%%%%%
% annotation('textbox',[plot_pos(1), plot_pos(2)-plot_pos(4)*1.5, plot_pos(3), plot_pos(4)/2],...
%     'String',{'Boundary Normal:',v(1,3),v(2,3),v(3,3)},...
%     'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 18);
% annotation('textbox',[plot_pos(1)+plot_pos(3)/2, plot_pos(2)-plot_pos(4)*1.5, plot_pos(3), plot_pos(4)/2],...
%     'String',{'Eigenvalues:',strcat('tmax = ',num2str(lambda(1))),...
%     strcat('tmid = ',num2str(lambda(2))),strcat('tmin = ',num2str(lambda(3)))},...
%     'VerticalAlignment','Bottom','Edgecolor','none','FontSize', 18);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





plot_name =  strcat('mms',probe_num,'_mvab_nested_',date_start(1:19),'_',num2str(duration),'.pdf');

print(gcf, '-dpdf', '-opengl',plot_name,'-fillpage');

movefile(plot_name, '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis')

toc
