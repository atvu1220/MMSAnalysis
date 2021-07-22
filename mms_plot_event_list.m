%Plot Events in list
cd '/Users/andrewvu/Library/Mobile Documents/com~apple~CloudDocs/Research/MATLAB Analysis'
%MMS Master Plot
%Plots summary and boundary analysis entirely on an event
warning('off','MATLAB:polyfit:RepeatedPointsOrRescale')
clear
tic
Typical = [ 2,3,5,6,8,11,12,14,16,17,18,20,21,24,29,31,32,35,37,38,41,42,43,48,50,52,54,55,58,60,65,67,68,76,79,81,82,83,84,85,93,94,98,100,...
        105,111,116,117,121,122,123,124,127,132,133,135,136,138,142,153,154,157,158,159,160,161,162,164,165,166,167,168,169,170];
Typical2 = [2;3;5;6;7;8;9;11;12;14;15;16;17;18;20;21;22;24;25;28;29;30;31;32;34;35;36;37;38;39;41;42;43;47;48;50;51;52;54;55;56;57;58;60;61;65;66;67;68;76;77;78;79;80;81;82;83;84;85;87;93;94;97;98;99;100;104;105;106;111;116;117;121;122;123;124;127;130;132;133;135;136;138;142;    152;153;154;157;158;159;160;161;162;164;165;166;167;168;169;170]';
Typical7=  [2;3;5;6;7;8;9;11;12;14;15;16;17;18;20;21;22;24;25;28;29;30;31;32;34;35;36;37;38;39;41;42;43;47;48;50;51;52;54;55;56;57;58;60;61;65;66;67;68;76;77;78;79;80;81;82;83;84;85;87;93;94;97;98;99;100;104;105;106;111;116;117;121;122;123;124;127;130;132;133;135;136;138;142;143;152;153;154;157;158;159;160;161;162;164;165;166;167;168;169;170;...
    179;183;184;185;187;188;189;191;199;200;201;202;203;204;205]';    % Loop for n events.
All_Events = [2:174,179,180,181,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,199,200,201,202,203,204,205,207];
All_Events2 = [2:48,50:68,70:72,74:90,92:108,111:112,115:118,121:125,127,130:136,138,142:143,150:154,156:172,174,179:181,183:197,199:205,207];
All_Events3 = [2:26,28:48,50:68,74:89,92:108,111:112,115:118,121:125,127,130:136,138,142:143,150:154,156:172,174,180:181,183:197,199:205,207];
%questionable = [16*,27,43,70,71,72,75,90,179]
for Event_number=All_Events3;%:206%Typical7%'CS1';%48%7%154%Typical7 %Typical7 ;%Typical[1:176,179,180,181,183,184,185,186,187,188,189,190,191,192,193,194,195,196,197,199,200,201,202,203,204,205,207]
   %2:174
   %Not
   %49,69:73,91,109,110,113,114,119,120,126,128,129,137,139,140,141,144,145,146,147,148,149,155,173
Event_number
close all
   [Event_Type, Substructure,threshold_std,data_type,...
          event_start,event_end,...
          left_OuterEdge, left_InnerEdge,...
          right_InnerEdge, right_OuterEdge,...
          leading_leftmost_date, leading_rightmost_date,...
          trailing_leftmost_date, trailing_rightmost_date,...
          leading_start, leading_end,...
          trailing_start, trailing_end,...
          downstream_geometry,upstream_geometry] = get_eventTimes(Event_number);
    

%% Perform Analysis    
    %What kind of analysis?
    %-2 is to get magnetic fields before and after 15 minutes.
    %-1 is for current normal and summary only
    %0 is for sliding window & Timing
    %0.5 is for sliding window and boundary analysis
    %1 is for boundary analysis
    %1.5 is for current normal, summary, MVAB sliding, boundary, and save
    %parameters to excel.
    %2 is for boundary analysis and save parameters to excel
    %3 is just for summary plot with n_cs, shear angle and duration in all events folder
    %4 is just for summary plot, with no timing analysis
    
    event_analysis = 2; 
    timing_window = 0;
    sliding_minratio = 10;
    n=1; %number of data points on each side for additional averaging for OMNI data. 1 data point = 5minutes
    
    plot_master(event_analysis,threshold_std,timing_window,sliding_minratio,data_type,...
        Event_number,Event_Type, Substructure,...
        event_start, event_end,...
        left_OuterEdge, left_InnerEdge,...
        right_InnerEdge, right_OuterEdge,...
        leading_leftmost_date, leading_rightmost_date,...
        trailing_leftmost_date,trailing_rightmost_date,...
        leading_start,leading_end,...
        trailing_start, trailing_end,...
        downstream_geometry,upstream_geometry,...
        n)
    

    
end
toc