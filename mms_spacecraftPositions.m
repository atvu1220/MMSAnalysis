%Plots MMS Spacecraft Positiosns from Mec Data
clear all
close all
event_start= '2019-02-26 00:32:30.000';
event_end = '2019-02-26 00:33:40.000';

event_start= '2019-02-22 12:01:00.000';
event_end = '2019-02-22 12:02:00.000';

event_start= '2019-02-02 18:09:00.000';
event_end = '2019-02-02 18:09:40.000';

event_start= '2019-01-28 04:45:40.000';
event_end = '2019-01-28 04:46:30.000';

event_start= '2019-02-23 06:04:50.000';
event_end = '2019-02-23 06:05:20.000';

[mms1_mec_timedata_raw, mms1_mec_rdata_raw] = load_mec(event_start,event_end,1,'srvy');
[mms2_mec_timedata_raw, mms2_mec_rdata_raw] = load_mec(event_start,event_end,2,'srvy');
[mms3_mec_timedata_raw, mms3_mec_rdata_raw] = load_mec(event_start,event_end,3,'srvy');
[mms4_mec_timedata_raw, mms4_mec_rdata_raw] = load_mec(event_start,event_end,4,'srvy');

mms1pos = (mms1_mec_rdata_raw);%/6371.2)
mms2pos = (mms2_mec_rdata_raw);%/6371.2)
mms3pos = (mms3_mec_rdata_raw);%/6371.2)
mms4pos = (mms4_mec_rdata_raw);%/6371.2)

% for i=1:length(mms1pos)

scatter(mms1pos(1,1),mms1pos(1,2),200,'ko'); hold on
scatter(mms2pos(1,1),mms2pos(1,2),200,'bp');
scatter(mms3pos(1,1),mms3pos(1,2),200,'gs');
scatter(mms4pos(1,1),mms4pos(1,2),200,'rd');
drawnow

scatter(mms1pos(end,1),mms1pos(end,2),200,'ko','filled'); hold on
scatter(mms2pos(end,1),mms2pos(end,2),200,'bp','filled');
scatter(mms3pos(end,1),mms3pos(end,2),200,'gs','filled');
scatter(mms4pos(end,1),mms4pos(end,2),200,'rd','filled');
drawnow

w = waitforbuttonpress;
title('MMS Position','fontsize',20)
legend({'MMS1';'MMS2';'MMS3';'MMS4'},'fontsize',18,'location','northwest')
box on
set(gca,'XMinorTick','on','TickDir','out','YMinorTick','on','linewidth',2)
set(gcf,'color','w');
xlabel({'X';'[km]'},'fontsize',18)
ylabel({'Y';'[km]'},'fontsize',18)
set(gca, 'xdir','reverse')
% end

