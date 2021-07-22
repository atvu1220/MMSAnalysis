%Calculates the distances between the 4spacecraft that are in a linear formation.
%Enter linear order of spacecraft from most sunward to most earthward
scOrder = [3,4,1,2];
trange = {'2019-02-16 05:40:25.000','2019-02-16 05:42:25.000'}
event_start = cell2mat(trange(1));
event_end = cell2mat(trange(2));

[mms1_mec_timedata_raw, mms1_mec_rdata_raw] = load_mec(event_start,1,'srvy');
[mms2_mec_timedata_raw, mms2_mec_rdata_raw] = load_mec(event_start,2,'srvy');
[mms3_mec_timedata_raw, mms3_mec_rdata_raw] = load_mec(event_start,3,'srvy');
[mms4_mec_timedata_raw, mms4_mec_rdata_raw] = load_mec(event_start,4,'srvy');


[~,rdata1,~,~] = crop(mms1_mec_timedata_raw,mms1_mec_rdata_raw,event_start,event_end);
[~,rdata2,~,~] = crop(mms2_mec_timedata_raw,mms2_mec_rdata_raw,event_start,event_end);
[~,rdata3,~,~] = crop(mms3_mec_timedata_raw,mms3_mec_rdata_raw,event_start,event_end);
[~,rdata4,~,~] = crop(mms4_mec_timedata_raw,mms4_mec_rdata_raw,event_start,event_end);


rdata1 = mean(rdata1);
rdata2 = mean(rdata2);
rdata3 = mean(rdata3);
rdata4 = mean(rdata4);

rdata1_mag = norm(rdata1)
rdata2_mag = norm(rdata2)
rdata3_mag = norm(rdata3)
rdata4_mag = norm(rdata4)

distancefrom21 = norm(rdata2-rdata1)
distancefrom14 = norm(rdata1-rdata4)
distancefrom43 = norm(rdata4-rdata3)

[~,maxIndex] = max([rdata1(1),rdata2(1),rdata3(1),rdata4(1)])
[~,minIndex] = min([rdata1(1),rdata2(1),rdata3(1),rdata4(1)])

