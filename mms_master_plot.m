%MMS Master Plot
%Plots summary and boundary analysis entirely on an event
clear
close all


%% Event Number 0 ; blue tags are events already entered in excel.
%add 5 seconds before and after outer boundaries for TDnormal means.
Event_number = 0; %For row entry in Excel Spreadsheet
Substructure = 0; %0 for no substructures, 1 for yes substructures, for differ sheets on excel
%-2 is to get magnetic fields before and after 15 minutes.
%-1 is for current normal and summary only
%0 is for sliding window
%0.5 is for sliding window and boundary analysis
%1 is for boundary analysis
%1.5 is for current normal, summary, MVAB sliding, boundary, and save
%parameters to excel.
%2 is for boundary analysis and save parameters to excel
event_analysis = 1;

%timing window %estimate it based on the time length of the gradient.
timing_window = 256*0.25;

%threshold for magnetic field mean for n_cs
threshold_std = 2.75;
%t'',''['','']trange = ['',''] ;FB
%trange = ['',''] ;
%% % %Event 76 %HFA
Event_number = 76;
Substructure = 1;
threshold_std = 0;

event_start    = '2018-02-17 19:45:35.000';
event_end      = '2018-02-17 19:46:50.000';

leading_leftmost_date   = '2018-02-17 19:45:35.000';
leading_rightmost_date  = '2018-02-17 19:45:43.000';
trailing_leftmost_date  = '2018-02-17 19:46:45.000';
trailing_rightmost_date = '2018-02-17 19:47:30.000';

leading_start  = '2018-02-17 19:46:07.594';
leading_end    = '2018-02-17 19:46:10.594'; 
					
trailing_start = '2018-02-17 19:46:22.852';		
trailing_end   = '2018-02-17 19:46:26.852';
		
%% % %Event 75 %HFA
% Event_number = 75;
% Substructure = 1;
% threshold_std = 2.75;
% 
% event_start    = '2018-03-29 03:25:45.000';
% event_end      = '2018-03-29 03:26:35.000';
% 
% leading_leftmost_date   = '2018-03-29 03:17:40.000';
% leading_rightmost_date  = '2018-03-29 03:24:40.000';
% trailing_leftmost_date  = '2018-03-29 03:26:40.000';
% trailing_rightmost_date = '2018-03-29 03:28:20.000';
% 
% leading_start  = '2018-03-29 03:25:53.481';
% leading_end    = '2018-03-29 03:25:54.481'; 
% 				
% trailing_start = '2018-03-29 03:26:23.013';		
% trailing_end   = '2018-03-29 03:26:25.013';
	
%% % %Event 74 %FB
% Event_number = 74;
% Substructure = 0;
% threshold_std = 0;
% 
% event_start    = '2018-01-29 03:35:45.000';
% event_end      = '2018-01-29 03:36:45.000';
% 
% leading_leftmost_date   = '2018-01-29 03:20:45.000';
% leading_rightmost_date  = '2018-01-29 03:35:50.000';
% trailing_leftmost_date  = '2018-01-29 03:36:25.000';
% trailing_rightmost_date = '2018-01-29 03:36:35.000';
% 
% leading_start  = '2018-01-29 03:35:51.472';
% leading_end    = '2018-01-29 03:35:54.472'; 
% 			
% trailing_start = '2018-01-29 03:36:15.027';		
% trailing_end   = '2018-01-29 03:36:17.027';

%% % %Event 73 %FB
% Event_number = 73;
% Substructure = 0;
% threshold_std = 2.75;
% 
% event_start    = '2018-01-15 16:54:30.000';
% event_end      = '2018-01-15 16:55:10.000';
% 
% leading_leftmost_date   = '2018-01-15 16:54:30.000';
% leading_rightmost_date  = '2018-01-15 16:54:30.000';
% trailing_leftmost_date  = '2018-01-15 16:55:10.000';
% trailing_rightmost_date = '2018-01-15 16:55:10.000';
% 
% leading_start  = '2018-01-15 16:54:35.726';
% leading_end    = '2018-01-15 16:54:39.726'; 
% 		
% trailing_start = '2018-01-15 16:54:54.546';		
% trailing_end   = '2018-01-15 16:54:55.546';

%% % %Event 72 %FB
% Event_number = 72;
% Substructure = 0;
% threshold_std = 0;
% 
% event_start    = '2018-01-12 02:59:15.000';
% event_end      = '2018-01-12 03:00:10.000';
% 
% leading_leftmost_date   = '2018-01-12 02:45:30.000';
% leading_rightmost_date  = '2018-01-12 02:59:20.000';
% trailing_leftmost_date  = '2018-01-12 03:00:10.000';
% trailing_rightmost_date = '2018-01-12 03:07:00.000';
% 
% leading_start  = '2018-01-12 02:59:23.250';
% leading_end    = '2018-01-12 02:59:34.250'; 
% 		
% trailing_start = '2018-01-12 02:59:51.027';		
% trailing_end   = '2018-01-12 02:59:54.027';

%% % %Event 71 %HFA
% Event_number = 71;
% Substructure = 0;
% threshold_std = 0;
% 
% event_start    = '2017-12-18 12:08:17.000';
% event_end      = '2017-12-18 12:10:55.000';
% 
% leading_leftmost_date   = '2017-12-18 12:08:28.000';
% leading_rightmost_date  = '2017-12-18 12:08:32.000';
% trailing_leftmost_date  = '2017-12-18 12:12:30.000';
% trailing_rightmost_date = '2017-12-18 12:21:00.000';
% 
% leading_start  = '2017-12-18 12:08:58.127';
% leading_end    = '2017-12-18 12:09:00.127'; 
% 	
% trailing_start = '2017-12-18 12:09:50.065';		
% trailing_end   = '2017-12-18 12:09:53.065';

%% % %Event 70
% Event_number = 70;
% Substructure = 0;
% threshold_std = 0;
% 
% event_start    = '2018-01-07 05:50:30.000';
% event_end      = '2018-01-07 05:53:15.000';
% 
% leading_leftmost_date   = '2018-01-07 05:35:00.000';
% leading_rightmost_date  = '2018-01-07 05:45:00.000';
% trailing_leftmost_date  = '2018-01-07 05:52:50.000';
% trailing_rightmost_date = '2018-01-07 05:53:14.000';
% 
% leading_start  = '2018-01-07 05:52:08.216';
% leading_end    = '2018-01-07 05:52:11.216'; 
% 	
% trailing_start = '2018-01-07 05:52:33.036';		
% trailing_end   = '2018-01-07 05:52:35.036';	
	
%% % %Event 69 
% Event_number = 69;
% Substructure = 0;
% threshold_std = 0;
% 
% event_start    = '2018-01-07 05:50:30.000';
% event_end      = '2018-01-07 05:53:15.000';
% 
% leading_leftmost_date   = '2018-01-07 05:35:00.000';
% leading_rightmost_date  = '2018-01-07 05:45:00.000';
% trailing_leftmost_date  = '2018-01-07 05:52:50.000';
% trailing_rightmost_date = '2018-01-07 05:53:14.000';
% 
% leading_start  = '2018-01-07 05:51:10.652';
% leading_end    = '2018-01-07 05:51:13.653'; 
% 	
% trailing_start = '2018-01-07 05:51:49.520';		
% trailing_end   = '2018-01-07 05:51:51.520';	
				
%% % %Event 68 %Potential Reconnection
% Event_number = 68;
% Substructure = 1;
% threshold_std = 0;
% 
% event_start    = '2017-12-15 14:24:50.000';
% event_end      = '2017-12-15 14:26:10.000';
% 
% leading_leftmost_date   = '2017-12-15 14:18:20.000';
% leading_rightmost_date  = '2017-12-15 14:20:30.000';
% trailing_leftmost_date  = '2017-12-15 14:26:40.000';
% trailing_rightmost_date = '2017-12-15 14:30:20.000';
% 
% leading_start  = '2017-12-15 14:25:24.286';
% leading_end    = '2017-12-15 14:25:26.286'; 
% 
% trailing_start = '2017-12-15 14:25:47.286';		
% trailing_end   = '2017-12-15 14:25:48.286';	
% 		
%% % %Event 67
% Event_number = 67;
% Substructure = 1;
% threshold_std = 0;
% 
% event_start    = '2018-02-20 13:42:00.000';
% event_end      = '2018-02-20 13:46:30.000';
% 
% leading_leftmost_date   = '2018-02-20 13:30:00.000';
% leading_rightmost_date  = '2018-02-20 13:42:00.000';
% trailing_leftmost_date  = '2018-02-20 13:46:00.000';
% trailing_rightmost_date = '2018-02-20 13:52:55.000';
% 
% leading_start  = '2018-02-20 13:44:44.731';
% leading_end    = '2018-02-20 13:44:45.731'; 
% 
% trailing_start = '2018-02-20 13:45:45.904';		
% trailing_end   = '2018-02-20 13:45:48.904';	

%% % %Event 66
% Event_number = 66;
% Substructure = 0;
% threshold_std = 0;
% 
% event_start    = '2018-02-20 13:42:00.000';
% event_end      = '2018-02-20 13:46:30.000';
% 
% leading_leftmost_date   = '2018-02-20 13:30:00.000';
% leading_rightmost_date  = '2018-02-20 13:42:00.000';
% trailing_leftmost_date  = '2018-02-20 13:46:00.000';
% trailing_rightmost_date = '2018-02-20 13:52:55.000';
% 								
% leading_start  = '2018-02-20 13:43:49.332';
% leading_end    = '2018-02-20 13:43:53.332'; 			
% 
% trailing_start = '2018-02-20 13:44:29.309';		
% trailing_end   = '2018-02-20 13:44:31.309';	

%% % %Event 65
% Event_number = 65;
% Substructure = 0;
% threshold_std = 0;
% 
% event_start    = '2018-02-20 13:42:00.000';
% event_end      = '2018-02-20 13:46:30.000';
% 
% leading_leftmost_date   = '2018-02-20 13:30:00.000';
% leading_rightmost_date  = '2018-02-20 13:42:00.000';
% trailing_leftmost_date  = '2018-02-20 13:46:00.000';
% trailing_rightmost_date = '2018-02-20 13:52:55.000';
% 								
% leading_start  = '2018-02-20 13:42:56.566';
% leading_end    = '2018-02-20 13:42:57.066'; 			
% 	
% trailing_start = '2018-02-20 13:43:37.215';		
% trailing_end   = '2018-02-20 13:43:39.215';	
			
%% % %Event 64 %FB
% Event_number = 64;
% Substructure = 0;
% threshold_std = 0;
% 
% event_start    = '2018-01-07 05:53:00.000';
% event_end      = '2018-01-07 05:54:10.000';
% 
% leading_leftmost_date   = '2018-01-07 05:52:50.000';
% leading_rightmost_date  = '2018-01-07 05:53:14.000';
% trailing_leftmost_date  = '2018-01-07 05:53:50.000';
% trailing_rightmost_date = '2018-01-07 06:09:10.000';
% trailing_rightmost_date = '2018-01-07 05:54:06.000';
% 								
% leading_start  = '2018-01-07 05:53:13.826';
% leading_end    = '2018-01-07 05:53:15.826'; 			
% 
% trailing_start = '2018-01-07 05:53:23.365';		
% trailing_end   = '2018-01-07 05:53:26.365';	

%% % %Event 63
% Event_number = 63;
% Substructure = 0;
% threshold_std = 0;
% 
% event_start    = '2018-01-07 05:50:30.000';
% event_end      = '2018-01-07 05:53:15.000';
% 
% leading_leftmost_date   = '2018-01-07 05:40:00.000';
% leading_rightmost_date  = '2018-01-07 05:45:00.000';
% trailing_leftmost_date  = '2018-01-07 05:52:50.000';
% trailing_rightmost_date = '2018-01-07 05:53:14.000';
% 								
% leading_start  = '2018-01-07 05:52:08.216';		
% leading_end    = '2018-01-07 05:52:11.216'; 			
% 								
% trailing_start = '2018-01-07 05:52:33.036';		
% trailing_end   = '2018-01-07 05:52:35.036';	

%% % %Event 62
% Event_number = 62;
% Substructure = 0;
% threshold_std = 0;
% 
% event_start    = '2018-01-07 05:50:30.000';
% event_end      = '2018-01-07 05:53:15.000';
% 
% leading_leftmost_date   = '2018-01-07 05:40:00.000';
% leading_rightmost_date  = '2018-01-07 05:45:00.000';
% trailing_leftmost_date  = '2018-01-07 05:52:50.000';
% trailing_rightmost_date = '2018-01-07 05:53:14.000';
% 								
% leading_start  = '2018-01-07 05:51:04.082';		
% leading_end    = '2018-01-07 05:51:08.082'; 				
% 								
% trailing_start = '2018-01-07 05:51:33.309';		
% trailing_end   = '2018-01-07 05:51:35.309';

%% % %Event 61
% Event_number = 61;
% Substructure = 1;
% threshold_std = 0;
% 
% event_start    = '2017-12-15 14:24:50.000';
% event_end      = '2017-12-15 14:26:10.000';
% 
% leading_leftmost_date   = '2017-12-15 14:23:30.000';
% leading_rightmost_date  = '2017-12-15 14:24:40.000';
% trailing_leftmost_date  = '2017-12-15 14:26:40.000';
% trailing_rightmost_date = '2017-12-15 14:28:40.000';
% 								
% leading_start  = '2017-12-15 14:25:37.216';		
% leading_end    = '2017-12-15 14:25:38.216'; 		
% 								
% trailing_start = '2017-12-15 14:25:44.223';		
% trailing_end   = '2017-12-15 14:25:47.224';	

%% % %Event 60
% Event_number = 60;
% Substructure = 1;
% threshold_std = 0;
% 
% event_start    = '2017-12-15 14:24:50.000';
% event_end      = '2017-12-15 14:26:10.000';
% 
% leading_leftmost_date   = '2017-12-15 14:23:30.000';
% leading_rightmost_date  = '2017-12-15 14:24:40.000';
% trailing_leftmost_date  = '2017-12-15 14:26:40.000';
% trailing_rightmost_date = '2017-12-15 14:28:40.000';
% 
% leading_start  = '2017-12-15 14:25:31.403';		
% leading_end    = '2017-12-15 14:25:31.903'; 	
% 
% trailing_start = '2017-12-15 14:25:36.809';		
% trailing_end   = '2017-12-15 14:25:37.809';

%% % %Event 59
% Event_number = 59;
% Substructure = 0;
% threshold_std = 0;
% 
% event_start    = '2017-12-15 14:24:50.000';
% event_end      = '2017-12-15 14:26:10.000';
% 
% leading_leftmost_date   = '2017-12-15 14:23:30.000';
% leading_rightmost_date  = '2017-12-15 14:24:40.000';
% trailing_leftmost_date  = '2017-12-15 14:26:40.000';
% trailing_rightmost_date = '2017-12-15 14:28:40.000';
% 								
% leading_start  = '2017-12-15 14:25:21.754';	
% leading_end    = '2017-12-15 14:25:24.754'; 	
% 							
% trailing_start = '2017-12-15 14:25:30.559';		
% trailing_end   = '2017-12-15 14:25:31.559';

%% % %Event 58
% Event_number = 58;
% Substructure = 1;
% threshold_std = 0;
% 
% event_start    = '2018-02-21 11:34:30.000';
% event_end      = '2018-02-21 11:37:20.000';
% 
% leading_leftmost_date   = '2018-02-21 11:25:00.000';
% leading_rightmost_date  = '2018-02-21 11:30:00.000';
% trailing_leftmost_date  = '2018-02-21 11:37:15.000';
% trailing_rightmost_date = '2018-02-21 11:40:00.000';
% 								
% leading_start  = '2018-02-21 11:35:03.498';
% leading_end    = '2018-02-21 11:35:04.498'; 	
% 							
% trailing_start = '2018-02-21 11:36:21.054';	
% trailing_end   = '2018-02-21 11:36:25.054';

%% % %Event 57
% Event_number = 57;
% Substructure = 0;
% threshold_std = 1.75;
% 
% event_start    = '2018-03-01 00:29:50.000';
% event_end      = '2018-03-01 00:31:10.000';
% 								
% leading_start  = '2018-03-01 00:30:23.776';
% leading_end    = '2018-03-01 00:30:27.776';
% 							
% trailing_start = '2018-03-01 00:30:52.550';	
% trailing_end   = '2018-03-01 00:30:55.550';	

%% % %Event 56 %SHFA
% Event_number = 56;
% Substructure = 1;
% threshold_std = 2.75;
% 
% event_start    = '2018-03-01 00:31:40.000';
% event_end      = '2018-03-01 00:32:25.000';
% 								
% leading_start  = '2018-03-01 00:31:51.590';
% leading_end    = '2018-03-01 00:31:54.590';
% 								
% trailing_start = '2018-03-01 00:32:09.020';		
% trailing_end   = '2018-03-01 00:32:10.020';
	
%% % %Event 55
% Event_number = 55;
% Substructure = 1;
% threshold_std = 4.75;
% 
% 
% event_start    = '2018-03-18 00:10:10.000';
% event_end      = '2018-03-18 00:11:05.000';
% 								
% leading_start  = '2018-03-18 00:10:38.597';
% leading_end    = '2018-03-18 00:10:39.597';
% 							
% trailing_start = '2018-03-18 00:10:47.699';		
% trailing_end   = '2018-03-18 00:10:49.699';

%% % %Event 54 %foreshock bubble, maybe RC
% Event_number = 54;
% Substructure = 1;
% threshold_std = 3.75;
% 
% event_start    = '2018-03-30 08:40:45.000';
% event_end      = '2018-03-30 08:42:45.000';
% 									
% leading_start  = '2018-03-30 08:41:00.583';
% leading_end    = '2018-03-30 08:41:04.583';
% 							
% trailing_start = '2018-03-30 08:42:18.459';		
% trailing_end   = '2018-03-30 08:42:21.460';
% 
% trailing_start = '2018-03-30 08:42:25.569';		
% trailing_end   = '2018-03-30 08:42:27.569';

%% % %Event 53
% Event_number = 53;
% Substructure = 1;
% threshold_std = 2.75;
% 
% event_start    = '2017-12-18 12:04:50.000';
% event_end      = '2017-12-18 12:06:35.000';
% 								
% leading_start  = '2017-12-18 12:05:24.420';
% leading_end    = '2017-12-18 12:05:25.420';
% 							
% trailing_start = '2017-12-18 12:06:01.827';		
% trailing_end   = '2017-12-18 12:06:03.827';

%% % %Event 53 %Foreshock Bubble
% Event_number = 53;
% Substructure = 0;	
% %threshold_std = 5;
% threshold_std = -0.5;
% 
% event_start    = '2018-04-03 13:57:10.000';
% event_end      = '2018-04-03 13:58:40.000';
% 							
% leading_start  = '2018-04-03 13:57:32.259';
% leading_end    = '2018-04-03 13:57:34.259';
% 									
% trailing_start = '2018-04-03 13:58:01.447';		
% trailing_end   = '2018-04-03 13:58:04.447';
	
%% % %Event 52
% Event_number = 52;
% Substructure = 0;	
% %threshold_std = 5;
% threshold_std = 2.00;
% 
% event_start    = '2018-03-12 07:35:50.000';
% event_end      = '2018-03-12 07:36:30.000';
% 						
% leading_start  = '2018-03-12 07:36:05.102';
% leading_end    = '2018-03-12 07:36:08.102';
% 									
% trailing_start = '2018-03-12 07:36:18.727';		
% trailing_end   = '2018-03-12 07:36:19.727';

%% % %Event 51, connected with Event 50,49,48,2
% Event_number = 51;
% Substructure = 1;	
% %threshold_std = 5;
% threshold_std = 2.75;
% 
% event_start    = '2018-03-01 01:06:15.000';
% event_end      = '2018-03-01 01:07:36.000';
% 					
% leading_start  = '2018-03-01 01:06:44.360';
% leading_end    = '2018-03-01 01:06:48.360';
% 
% trailing_start = '2018-03-01 01:07:01.517';		
% trailing_end   = '2018-03-01 01:07:03.517';

%% % %Event 50, connected with Event 49
% Event_number = 50;
% Substructure = 1;	
% %threshold_std = 5;
% threshold_std = 3;
% 
% %Overall Events
% event_start    = '2018-03-01 01:02:30.000';
% event_end      = '2018-03-01 01:08:40.000';
% 
% event_start    = '2018-03-01 01:05:40.000';
% event_end      = '2018-03-01 01:06:20.000';
% threshold_std = -0.25;
% 				
% leading_start  = '2018-03-01 01:05:53.860';
% leading_end    = '2018-03-01 01:05:54.860';
% 							
% trailing_start = '2018-03-01 01:06:10.352';	
% trailing_end   = '2018-03-01 01:06:11.352';
	
%% % %Event 49, connected with Event 50	
% Event_number = 49;
% Substructure = 0;	
% %threshold_std = 5;
% threshold_std = 3;
% 
% %Overall Events
% event_start    = '2018-03-01 01:02:30.000';
% event_end      = '2018-03-01 01:08:40.000';
% 
% event_start    = '2018-03-01 01:05:40.000';
% event_end      = '2018-03-01 01:06:20.000';
% threshold_std = -0.25;
% 
% leading_start  = '2018-03-01 01:05:50.469';
% leading_end    = '2018-03-01 01:05:50.969';
% 							
% trailing_start = '2018-03-01 01:05:53.860';	
% trailing_end   = '2018-03-01 01:05:54.860';
	
%% % %Event 48, connected with Event 2 	
% Event_number = 48;
% Substructure = 0;	
% %threshold_std = 5;
% threshold_std = 3;
% 
% %Overall Events
% event_start    = '2018-03-01 01:02:30.000';
% event_end      = '2018-03-01 01:08:40.000';
% 
% event_start    = '2018-03-01 01:04:30.000';
% event_end      = '2018-03-01 01:05:05.000';
% threshold_std = 3.00;
% 
% leading_start  = '2018-03-01 01:04:42.624';
% leading_end    = '2018-03-01 01:04:44.624';	
% 							
% trailing_start = '2018-03-01 01:04:51.359';	
% trailing_end   = '2018-03-01 01:04:53.359';

%% % %Event 47
% Event_number = 47;
% Substructure = 0;	
% %threshold_std = 5;
% threshold_std = -2.65;
% 
% event_start    = '2018-01-23 09:11:50.000';
% event_end      = '2018-01-23 09:12:50.000';
% 			
% leading_start  = '2018-01-23 09:12:05.643';
% leading_end    = '2018-01-23 09:12:06.643';
% 						
% trailing_start = '2018-01-23 09:12:22.557';	
% trailing_end   = '2018-01-23 09:12:26.557';
%% % %Event 46
% Event_number = 46;
% Substructure = 1;	
% threshold_std = 4;
% % threshold_std = -2;
% 
% event_start    = '2018-01-09 09:20:51.000';
% event_end      = '2018-01-09 09:21:45.000';
% 		
% leading_start  = '2018-01-09 09:21:15.868';
% leading_end    = '2018-01-09 09:21:16.868';
% 					
% trailing_start = '2018-01-09 09:21:26.384';	
% trailing_end   = '2018-01-09 09:21:28.384';
%% % %Event 45, connected 44-45
% Event_number = 45;
% Substructure = 1;
% threshold_std = 4;
% % threshold_std = -2;
% 
% event_start    = '2018-01-09 09:18:45.000';
% event_end      = '2018-01-09 09:20:00.000';
% 		
% leading_start  = '2018-01-09 09:19:18.953';
% leading_end    = '2018-01-09 09:19:22.953';
% 				
% trailing_start = '2018-01-09 09:19:43.469';	
% trailing_end   = '2018-01-09 09:19:45.469';
	
%% % %Event 44, connected 44-45
% Event_number = 44;
% Substructure = 0;
% threshold_std = 4;
% %threshold_std = -2;
% 
% event_start    = '2018-01-09 09:18:45.000';
% event_end      = '2018-01-09 09:20:00.000';
% 		
% leading_start  = '2018-01-09 09:18:52.976';
% leading_end    = '2018-01-09 09:18:55.976';
% 				
% trailing_start = '2018-01-09 09:19:12.593';	
% trailing_end   = '2018-01-09 09:19:16.593';

%% % %Event 43 %SHFA
% Event_number = 43;
% Substructure = 0;
% threshold_std = 2.25;
% % threshold_std = -2;
% 
% event_start    = '2017-12-29 19:22:40.000';
% event_end      = '2017-12-29 19:23:15.000';
% 		
% leading_start  = '2017-12-29 19:22:51.828';
% leading_end    = '2017-12-29 19:22:55.828';
% 			
% trailing_start = '2017-12-29 19:23:07.625';	
% trailing_end   = '2017-12-29 19:23:08.625';
	
%% % %Event 42
% Event_number = 42;
% Substructure = 0;
% threshold_std = 2.75;
% % threshold_std = -2;
% 
% event_start    = '2017-12-29 17:29:45.000';
% event_end      = '2017-12-29 17:30:35.000';
% 		
% leading_start  = '2017-12-29 17:30:03.874';
% leading_end    = '2017-12-29 17:30:05.874';
% 			
% trailing_start = '2017-12-29 17:30:15.570';	
% trailing_end   = '2017-12-29 17:30:17.570';
			
%% % %Event 41
% Event_number = 41;
% Substructure = 0;
% threshold_std = 2.85;
% threshold_std = -2;
% 
% event_start    = '2017-12-18 14:28:40.000';
% event_end      = '2017-12-18 14:30:35.000';
% 		
% leading_start  = '2017-12-18 14:30:02.318';
% leading_end    = '2017-12-18 14:30:05.318';
% 			
% trailing_start = '2017-12-18 14:30:16.655';	
% trailing_end   = '2017-12-18 14:30:19.655';
	
%% % %Event 40
% Event_number = 40;
% Substructure = 0;
% threshold_std = 2.85;
% threshold_std = -2;
% 
% event_start    = '2017-12-18 14:28:40.000';
% event_end      = '2017-12-18 14:30:35.000';
% 	
% leading_start  = '2017-12-18 14:29:24.583';
% leading_end    = '2017-12-18 14:29:28.584';
% 			
% trailing_start = '2017-12-18 14:29:51.435';	
% trailing_end   = '2017-12-18 14:29:54.435';

%% % %Event 39
% Event_number = 39;
% Substructure = 1;
% threshold_std = 2.85;
% threshold_std = -2;
% 
% event_start    = '2017-12-18 14:28:40.000';
% event_end      = '2017-12-18 14:30:35.000';
% 	
% leading_start  = '2017-12-18 14:29:10.779';
% leading_end    = '2017-12-18 14:29:12.779';
% 			
% trailing_start = '2017-12-18 14:29:25.334';
% trailing_end   = '2017-12-18 14:29:28.334';
	
%% % %Event 38
% Event_number = 38;
% Substructure = 1;
% threshold_std = 2.85;
% threshold_std = -2;
% 
% event_start    = '2017-12-18 14:28:40.000';
% event_end      = '2017-12-18 14:30:35.000';
% 
% leading_start  = '2017-12-18 14:28:57.458';
% leading_end    = '2017-12-18 14:28:58.458';
% 
% trailing_start = '2017-12-18 14:29:06.193';
% trailing_end   = '2017-12-18 14:29:10.193';

%% % %Event 37, connected 37-41
% Event_number = 37;
% Substructure = 0;
% threshold_std = 2.85;
% threshold_std = -2;
% 
% event_start    = '2017-12-18 14:28:40.000';
% event_end      = '2017-12-18 14:30:35.000';
% 
% leading_start  = '2017-12-18 14:28:46.880';
% leading_end    = '2017-12-18 14:28:50.880';
% 		
% trailing_start = '2017-12-18 14:28:57.458';
% trailing_end   = '2017-12-18 14:28:58.458';

%% % %Event 36, connected with Event 35
% Event_number = 36;
% Substructure = 0;
% threshold_std = 1.5;
% 
% event_start    = '2017-12-18 14:10:20.000';
% event_end      = '2017-12-18 14:11:30.000';
% 
% leading_start  = '2017-12-18 14:10:50.692';
% leading_end    = '2017-12-18 14:10:54.692';
% 		
% trailing_start = '2017-12-18 14:10:56.466';
% trailing_end   = '2017-12-18 14:10:58.466';
		
%% % %Event 35, connected with Event 36
% Event_number = 35;
% Substructure = 0;
% threshold_std = 1.5;
% 
% event_start    = '2017-12-18 14:10:20.000';
% event_end      = '2017-12-18 14:11:30.000';
% 
% leading_start  = '2017-12-18 14:10:57.318';
% leading_end    = '2017-12-18 14:11:00.318';
% 			
% trailing_start = '2017-12-18 14:11:07.443';
% trailing_end   = '2017-12-18 14:11:09.443';

%% % %Event 34 trange = ['','']
% Event_number = 34;
% Substructure = 0;
% threshold_std = 2.75;
% 
% event_start    = '2017-12-18 12:35:30.000';
% event_end      = '2017-12-18 12:36:30.000';
% 
% leading_start  = '2017-12-18 12:35:54.962';
% leading_end    = '2017-12-18 12:35:55.962';
% 	
% trailing_start = '2017-12-18 12:36:12.478';
% trailing_end   = '2017-12-18 12:36:14.478';
	
%% % %Event 33
% Event_number = 33;
% Substructure = 1;
% threshold_std = 2.00;
% 
% event_start    = '2017-12-15 14:22:05.000';
% event_end      = '2017-12-15 14:24:00.000';
% 
% leading_start  = '2017-12-15 14:22:38.947';
% leading_end    = '2017-12-15 14:22:41.947';
% 
% trailing_start = '2017-12-15 14:23:07.338';
% trailing_end   = '2017-12-15 14:23:09.338';

%% % %Event 32
% Event_number = 32;
% Substructure = 0;
% threshold_std = 3.75;
% 
% event_start    = '2017-12-10 00:42:10.000';
% event_end      = '2017-12-10 00:44:50.000';
% 
% leading_start  = '2017-12-10 00:42:24.729';
% leading_end    = '2017-12-10 00:42:28.729';
% 		
% trailing_start = '2017-12-10 00:44:11.441';
% trailing_end   = '2017-12-10 00:44:13.441';
	
%% % %Event 31
% Event_number = 31;
% Substructure = 0;
% threshold_std = 2.25;
% 
% event_start    = '2017-12-01 14:12:40.000';
% event_end      = '2017-12-01 14:13:40.000';
% 
% leading_start  = '2017-12-01 14:12:55.440';
% leading_end    = '2017-12-01 14:12:55.940';
% 		
% trailing_start = '2017-12-01 14:13:14.355';
% trailing_end   = '2017-12-01 14:13:15.355';

%% % %Event 30
% Event_number = 30;
% Substructure = 0;
% threshold_std = 3.5;
% 
% event_start    = '2017-12-01 14:11:20.000';
% event_end      = '2017-12-01 14:12:20.000';
% 
% leading_start  = '2017-12-01 14:11:34.275';
% leading_end    = '2017-12-01 14:11:38.275';
% 		
% trailing_start = '2017-12-01 14:12:01.643';
% trailing_end   = '2017-12-01 14:12:04.643';

%% % %Event 29 
% Event_number = 29;
% Substructure = 1;
% threshold_std = 2.75;
% 
% event_start    = '2017-11-20 08:06:30.000';
% event_end      = '2017-11-20 08:09:00.000';
% 
% leading_start  = '2017-11-20 08:07:06.555';
% leading_end    = '2017-11-20 08:07:07.555';
% 	
% trailing_start = '2017-11-20 08:08:14.931';
% trailing_end   = '2017-11-20 08:08:15.931';

%% % %Event 28
% Event_number = 28;
% Substructure = 0;
% threshold_std = 2.75;
% 
% event_start    = '2017-11-04 22:45:40.000';
% event_end      = '2017-11-04 22:46:45.000';
% 
% leading_start  = '2017-11-04 22:45:56.138';
% leading_end    = '2017-11-04 22:45:58.138';
% 
% trailing_start = '2017-11-04 22:46:24.466';
% trailing_end   = '2017-11-04 22:46:26.467';

%% % %Event 27
% Event_number = 27;
% Substructure = 0;
% threshold_std = 2.25;
% 
% event_start    = '2017-10-21 21:18:55.000';
% event_end      = '2017-10-21 21:20:00.000';
% 
% leading_start  = '2017-10-21 21:19:28.038';
% leading_end    = '2017-10-21 21:19:29.038';
% 
% trailing_start = '2017-10-21 21:19:32.882';
% trailing_end   = '2017-10-21 21:19:34.882';

%% % %Event 0 %too turbulent before and after/ foreshock bubble with no clear leading.
% Event_number = 0;
% Substructure = 0;
% threshold_std = 3.75;
%
% event_start    = '2017-01-14 06:49:15.000';
% event_end      = '2017-01-14 06:50:05.000';
%
% leading_start  = '2017-11-16 12:13:22.457';
% leading_end    = '2017-11-16 12:13:23.457';
%
% trailing_start = '2017-11-16 12:13:37.090';
% trailing_end   = '2017-11-16 12:13:37.590';

%% % %Event 26
% Event_number = 26;
% Substructure = 0; %slight substrcutures that divides interior into two regions
% threshold_std = 1.50;
% 
% event_start    = '2017-11-19 07:43:20.000';
% event_end      = '2017-11-19 07:44:08.000';
% 
% leading_start  = '2017-11-19 07:43:27.580';
% leading_end    = '2017-11-19 07:43:28.580';
% 
% trailing_start = '2017-11-19 07:43:58.041';
% trailing_end   = '2017-11-19 07:43:59.041';

%% % %Event 25 %Transverse Speed too high!
% Event_number = 25;
% Substructure = 1;
% threshold_std = 2.75;
% 
% event_start    = '2016-12-08 09:59:50.000';
% event_end      = '2016-12-08 10:01:05.000';
% 
% leading_start  = '2016-12-08 10:00:03.785';
% leading_end    = '2016-12-08 10:00:04.785';
% 
% trailing_start = '2016-12-08 10:00:40.059';
% trailing_end   = '2016-12-08 10:00:41.559';

%% % %Event 24
% Event_number = 24;
% Substructure = 1;
% threshold_std = 3.25;
% 
% event_start    = '2015-12-28 04:27:04.000';
% event_end      = '2015-12-28 04:28:34.000';
% 
% leading_start  = '2015-12-28 04:27:26.374';
% leading_end    = '2015-12-28 04:27:26.874';
% 
% trailing_start = '2015-12-28 04:27:58.413';
% trailing_end   = '2015-12-28 04:27:59.413';

%% % %Event 23 %SHFA
% Event_number = 23;
% Substructure = 0;
% threshold_std = 2.75;
% 
% event_start    = '2017-11-16 12:13:00.000';
% event_end      = '2017-11-16 12:14:00.000';
% 
% leading_start  = '2017-11-16 12:13:22.457';
% leading_end    = '2017-11-16 12:13:23.457';
% 
% 
% trailing_start = '2017-11-16 12:13:37.090';
% trailing_end   = '2017-11-16 12:13:37.590';

%% % %Event 20b/22
% Event_number = 22;
% Substructure = 0;
% threshold_std = 1;
% 
% event_start    = '2017-12-29 19:11:21.000';
% event_end      = '2017-12-29 19:11:50.000';
% 
% leading_start  = '2017-12-29 19:11:24.880';
% leading_end    = '2017-12-29 19:11:25.880';
% 
% trailing_start = '2017-12-29 19:11:31.974';
% trailing_end   = '2017-12-29 19:11:33.974';

%% % %Event 20a/21
% Event_number = 21;
% Substructure = 0;
% threshold_std = 0.1;
% 
% event_start    = '2017-12-29 19:11:00.000';
% event_end      = '2017-12-29 19:11:24.000';
% 
% leading_start  = '2017-12-29 19:11:12.576';
% leading_end    = '2017-12-29 19:11:13.076';
% 
% trailing_start = '2017-12-29 19:11:18.498';
% trailing_end   = '2017-12-29 19:11:19.498';

%% % %Event 20 
% Event_number = 20;
% Substructure = 0;
% threshold_std = 2.00;
% 
% event_start    = '2017-11-10 18:59:40.000';
% event_end      = '2017-11-10 19:01:26.000';
% 
% leading_start  = '2017-11-10 19:00:00.720';
% leading_end    = '2017-11-10 19:00:03.720';
% 
% trailing_start = '2017-11-10 19:01:06.635';
% trailing_end   = '2017-11-10 19:01:09.635';

%% % %Event 19 %foreshock bubble
% Event_number = 19;
% Substructure = 0;
% threshold_std = 2.75;
% 
% event_start    = '2018-04-03 13:57:25.000';
% event_end      = '2018-04-03 13:58:35.000';
% 
% leading_start  = '2018-04-03 13:57:32.095';
% leading_end    = '2018-04-03 13:57:34.095';
% 
% trailing_start = '2018-04-03 13:58:03.228';
% trailing_end   = '2018-04-03 13:58:05.228';

%% % %Event 18
% Event_number = 18;
% Substructure = 0;
% threshold_std = 2.5;
% 
% event_start    = '2017-12-29 17:20:35.000';
% event_end      = '2017-12-29 17:22:10.000';
% 
% leading_start  = '2017-12-29 17:21:01.570';
% leading_end    = '2017-12-29 17:21:02.570';
% 
% trailing_start = '2017-12-29 17:21:34.461';
% trailing_end   = '2017-12-29 17:21:34.961';

%% % %Event 17
% Event_number = 17;
% Substructure = 1;
% threshold_std = 2.25;
% 
% event_start    = '2018-04-27 19:47:30.000';
% event_end      = '2018-04-27 19:48:22.000';
% 
% leading_start  = '2018-04-27 19:47:45.500';
% leading_end    = '2018-04-27 19:47:47.500';
% 
% trailing_start = '2018-04-27 19:47:58.750';
% trailing_end   = '2018-04-27 19:47:59.750';

%% % %Event 16
% Event_number = 16;
% Substructure = 0;
% threshold_std = 3.75;
% 
% event_start = '2018-04-01 01:10:34.000';
% event_end = '2018-04-01 01:11:22.000';
% 
% leading_start  = '2018-04-01 01:10:50.500';
% leading_end    = '2018-04-01 01:10:54.500';
% 
% trailing_start = '2018-04-01 01:10:58.000';
% trailing_end   = '2018-04-01 01:11:00.000';

%% % %Event 15
% Event_number = 15;
% Substructure = 0;
% threshold_std = 3.25;
% 
% event_start    = '2017-12-18 12:55:30.000';
% event_end      = '2017-12-18 12:57:15.000';
% 
% leading_start  = '2017-12-18 12:55:58.753';
% leading_end    = '2017-12-18 12:56:01.753';
% 
% trailing_start = '2017-12-18 12:56:40.636';
% trailing_end   = '2017-12-18 12:56:42.636';

%% % %Event 14
% Event_number = 14;
% Substructure = 1;
% threshold_std = 2.75;
% 
% event_start = '2018-03-26 10:00:40.000';
% event_end = '2018-03-26 10:01:20.000';
% 
% leading_start  = '2018-03-26 10:00:57.850';
% leading_end    = '2018-03-26 10:01:01.850';
% 
% trailing_start = '2018-03-26 10:01:07.257';
% trailing_end   = '2018-03-26 10:01:09.257';

%% % %Event 13
% Event_number = 13;
% Substructure = 1;
% threshold_std = 2.75;
% 
% event_start = '2018-03-18 22:00:04.000';
% event_end = '2018-03-18 22:00:52.000';
% 
% leading_start  = '2018-03-18 22:00:26.000';
% leading_end    = '2018-03-18 22:00:27.000';
% 
% trailing_start = '2018-03-18 22:00:32.500';
% trailing_end   = '2018-03-18 22:00:33.500';

%% % %Event 12
% Event_number = 12;
% Substructure = 1;
% threshold_std = 2.75;
% 
% event_start = '2018-01-29 18:34:30.000';
% event_end = '2018-01-29 18:35:30.000';
% 
% leading_start  = '2018-01-29 18:34:57.000';
% leading_end    = '2018-01-29 18:34:57.750';
% 
% trailing_start = '2018-01-29 18:35:09.000';
% trailing_end   = '2018-01-29 18:35:12.000';

%% % %Event 11
% Event_number = 11;
% Substructure = 0;
% threshold_std = 5;
% 
% event_start = '2017-12-29 18:10:55.000';
% event_end = '2017-12-29 18:12:25.000';
% 
% leading_start  = '2017-12-29 18:11:16.750';
% leading_end    = '2017-12-29 18:11:18.500';
% 
% trailing_start = '2017-12-29 18:11:52.750';
% trailing_end   = '2017-12-29 18:11:54.750';

%% % %Event 10
% Event_number = 10; %mec srvy is broken, use brst.
% Substructure = 1;  %waves before event/after
% threshold_std = 2.75;
% event_start    = '2018-01-26 08:23:05.000';
% event_end      = '2018-01-26 08:23:35.000';
% 
% leading_start  = '2018-01-26 08:23:15.500';
% leading_end    = '2018-01-26 08:23:17.250';
% 
% trailing_start = '2018-01-26 08:23:25.250';
% trailing_end   = '2018-01-26 08:23:26.250';

%% % %Event 9
% Event_number = 9; %lots of wave activity
% Substructure = 1;
% threshold_std = 2.75;
% event_start    = '2018-01-14 23:06:50.000';
% event_end      = '2018-01-14 23:08:05.000';
% 
% leading_start  = '2018-01-14 23:07:02.000';
% leading_end    = '2018-01-14 23:07:02.500';
% 
% trailing_start = '2018-01-14 23:07:48.250';
% trailing_end   = '2018-01-14 23:07:51.000';

%% % %Event 8
% Event_number = 8;
% Substructure = 0;
% threshold_std = 2.75;
% 
% event_start    = '2018-01-12 01:50:21.000';
% event_end      = '2018-01-12 01:52:30.000';
% 
% leading_start  = '2018-01-12 01:50:37.250';
% leading_end    = '2018-01-12 01:50:38.000';
% 
% trailing_start = '2018-01-12 01:51:47.375';
% trailing_end   = '2018-01-12 01:51:48.250';

%% % %Event 7
% Event_number = 7;
% Substructure = 1;
% threshold_std = 2.5;
% event_start    = '2018-01-09 08:34:27.000';
% event_end      = '2018-01-09 08:35:02.000';
% 
% leading_start  = '2018-01-09 08:34:34.750';
% leading_end    = '2018-01-09 08:34:35.250';
% 
% trailing_start = '2018-01-09 08:34:49.500';
% trailing_end   = '2018-01-09 08:34:51.000';

%% % %Event 6
% Event_number = 6;
% Substructure = 1;
% threshold_std=3;
% event_start      = '2017-12-18 11:59:58.000';
% event_end        = '2017-12-18 12:01:20.000';
% 
% leading_start    = '2017-12-18 12:00:21.500';
% leading_end      = '2017-12-18 12:00:23.500';
% 
% trailing_start   = '2017-12-18 12:01:00.375';
% trailing_end     = '2017-12-18 12:01:01.875';

%% % %Event 5
%changing the event_start by 5 seconds causes event distance by 20 times...
%drastiically changes cs normal. by 17 degrees between 45 and 50.
% Event_number = 5;
% Substructure = 0; %srvy mec
% threshold_std = 2.75;
% 
% event_start = '2017-11-10 17:25:45.000';
% event_end = '2017-11-10 17:27:57.000';
% 
% leading_start = '2017-11-10 17:26:17.500';
% leading_end = '2017-11-10 17:26:19.000';
% 
% trailing_start = '2017-11-10 17:27:38.750';
% trailing_end = '2017-11-10 17:27:40.000';

%% % %Event 4
% Event_number = 4;
% Substructure = 0;
% threshold_std = 2.75;
% 
% event_start = '2018-04-07 17:52:40.000';
% event_end = '2018-04-07 17:53:15.000';
% 
% leading_start = '2018-04-07 17:52:51.000';
% leading_end = '2018-04-07 17:52:52.000';
% 
% trailing_start = '2018-04-07 17:53:04.000';
% trailing_end = '2018-04-07 17:53:05.000';

%% % %Event 3
% Event_number = 3; %Use Mec Srvy
% Substructure = 0;
% threshold_std = 2.75;
% 
% event_start = '2015-12-28 05:26:48.000';
% event_end = '2015-12-28 05:27:36.000';
% 
% leading_start = '2015-12-28 05:26:57.750';
% leading_end = '2015-12-28 05:26:59.000';
% 
% trailing_start = '2015-12-28 05:27:25.750';
% trailing_end = '2015-12-28 05:27:27.500';

%% % %Event 2
% Event_number = 2; %compression of wave activity after event, averages therefore would be inaccurate
% Substructure = 0;
% threshold_std = 2.125;
% 
% event_start = '2018-03-01 01:03:44.000';
% event_end = '2018-03-01 01:04:19.000';
% 
% leading_start = '2018-03-01 01:03:54.600';
% leading_end = '2018-03-01 01:03:55.100';
% 
% trailing_start = '2018-03-01 01:04:10.500';
% trailing_end = '2018-03-01 01:04:11.250';

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%Directory Management%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if event_analysis ~= 2
    cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis'
    eventFolderName = strcat(event_start(1:10),'/',event_start(12:13),'-',event_start(15:16));
    if ~exist(eventFolderName,'dir')
        mkdir(eventFolderName)
    end
    cd(eventFolderName)
    
elseif event_analysis == 2
    cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analyzed Events'
    eventFolderName = strcat('Event Number/',num2str(Event_number));
    if exist(eventFolderName,'dir')
        rmdir(eventFolderName,'s')
        mkdir(eventFolderName)
    else
        mkdir(eventFolderName)
    end
    cd(eventFolderName)
end



%% %Data Retrieval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Data Retrieval%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data_type = 'brst';
%Load MEC Data
[mms1_mec_timedata_raw, mms1_mec_rdata_raw] = load_mec(event_start,1,'srvy');
[mms2_mec_timedata_raw, mms2_mec_rdata_raw] = load_mec(event_start,2,'srvy');
[mms3_mec_timedata_raw, mms3_mec_rdata_raw] = load_mec(event_start,3,'srvy');
[mms4_mec_timedata_raw, mms4_mec_rdata_raw] = load_mec(event_start,4,'srvy');

%load FGM data
[mms1_fgm_timedata_raw, mms1_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,1,data_type);
[mms2_fgm_timedata_raw, mms2_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,2,data_type);
[mms3_fgm_timedata_raw, mms3_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,3,data_type);
[mms4_fgm_timedata_raw, mms4_fgm_bdata_raw, ~, ~] = load_fgm(event_start,event_end,4,data_type);

[mms1_fgm_timedata_srvy, mms1_fgm_bdata_srvy] = load_fgm(event_start,event_end,1,'srvy');

%Load FPI_e
[fpi_e_timedata,fpi_e_ndata,fpi_e_vdata,fpi_e_tparadata,fpi_e_tperpdata,...
    fpi_e_edata,fpi_e_espectdata] = load_fpi(event_start,event_end,1,data_type,'e');
%Load FPI_i
[fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_tparadata,fpi_i_tperpdata,...
    fpi_i_edata,fpi_i_espectdata] = load_fpi(event_start,event_end,1,data_type,'i');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% %Plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Plot Plots%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if threshold_std == 0
    %Current Sheet Normal Calculation Plot, manually.
    [n_cs,B_pre,B_post] = manualCurrentSheet(event_start,event_end,leading_leftmost_date,leading_rightmost_date,...
        trailing_leftmost_date,trailing_rightmost_date,mms1_fgm_timedata_srvy,mms1_fgm_bdata_srvy);
else
    %Calculate Current Sheet
    [n_cs,B_pre,B_post] = calculateCurrentSheet(event_start,event_end,mms1_fgm_timedata_srvy,mms1_fgm_bdata_srvy,threshold_std);
end

%Summary Plot
[bowshock_n,currentsheet_n,r_sc] = plot_mms_observation_summary_mod(event_start,event_end,...
    mms1_mec_timedata_raw,mms1_mec_rdata_raw,mms1_fgm_timedata_raw, mms1_fgm_bdata_raw,...
    fpi_e_timedata,fpi_e_ndata,fpi_e_vdata,fpi_e_tparadata,fpi_e_tperpdata,fpi_e_edata,fpi_e_espectdata,...
    fpi_i_timedata,fpi_i_ndata,fpi_i_vdata,fpi_i_tparadata,fpi_i_tperpdata,fpi_i_edata,fpi_i_espectdata,...
    n_cs)


%Either 0, do sliding windows, or 1, do MVA,MDD,Timing, not both., 1.5 is
%for MVA sliding window only.
if event_analysis == 0 || event_analysis == 0.5 || event_analysis == 1.5
    %Sliding Window MVA
    plot_MVABslidingwindow(event_start,event_end,mms1_fgm_timedata_srvy, mms1_fgm_bdata_srvy,currentsheet_n)
    
    if event_analysis == 1.5
    else
        %Sliding Window Timing
        [time_ranges] =  plot_Timingslidingwindow(event_start,event_end,...
            mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
            mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
            mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
            mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
            mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
            mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
            mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
            mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
            timing_window,data_type);
        
%         [time_ranges] =  plot_Timingslidingwindow(event_start,event_end,...
%             mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
%             mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
%             mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
%             mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
%             mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
%             mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
%             mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
%             mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
%             0.25*256,data_type);
%         
%         [time_ranges] =  plot_Timingslidingwindow(event_start,event_end,...
%             mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
%             mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
%             mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
%             mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
%             mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
%             mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
%             mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
%             mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
%             0.5*256,data_type);
%         
%             [time_ranges] =  plot_Timingslidingwindow(event_start,event_end,...
%             mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
%             mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
%             mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
%             mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
%             mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
%             mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
%             mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
%             mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
%             1*256,data_type);
%         
%         
%             [time_ranges] =  plot_Timingslidingwindow(event_start,event_end,...
%             mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
%             mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
%             mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
%             mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
%             mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
%             mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
%             mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
%             mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
%             1.5*256,data_type);
%         
%         
%             [time_ranges] =  plot_Timingslidingwindow(event_start,event_end,...
%             mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
%             mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
%             mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
%             mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
%             mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
%             mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
%             mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
%             mms4_mec_timedata_raw,mms4_mec_rdata_raw,...
%             2*256,data_type);
    end
end


if event_analysis == 1 || event_analysis == 0.5 || event_analysis == 1.5 || event_analysis == 2
    %Inner Leading Boundary
    [MVA_n1,timing_n1,timing_v1,MVAB_timing_angle1] = plot_complete_boundary_analysis(event_start,event_end,leading_start,leading_end,...
        mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
        mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
        mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
        mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
        mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
        mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
        mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
        mms4_mec_timedata_raw,mms4_mec_rdata_raw)
    
    %Inner Trailing Boundary
    [MVA_n2,timing_n2,timing_v2,MVAB_timing_angle2] = plot_complete_boundary_analysis(event_start,event_end,trailing_start,trailing_end,...
        mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,...
        mms2_fgm_timedata_raw,mms2_fgm_bdata_raw,...
        mms3_fgm_timedata_raw,mms3_fgm_bdata_raw,...
        mms4_fgm_timedata_raw,mms4_fgm_bdata_raw,...
        mms1_mec_timedata_raw,mms1_mec_rdata_raw,...
        mms2_mec_timedata_raw,mms2_mec_rdata_raw,...
        mms3_mec_timedata_raw,mms3_mec_rdata_raw,...
        mms4_mec_timedata_raw,mms4_mec_rdata_raw)
    
end
%% %Event Parameters
if event_analysis == 1.5 || event_analysis == 2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%Event Parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    duration = findDuration(leading_end,trailing_start);
    %currentsheet_n=currentsheet_n.*[-1,1,1]; %-X GSE for CS normal,
    currentsheet_n=currentsheet_n; %#ok<ASGSL> %+X GSE for CS normal,according to Schwartz et al. 2018
    %Find average values before and after event
    %[B_pre,B_post] = pre_post(mms1_fgm_timedata_raw,mms1_fgm_bdata_raw(:,1:3),event_start,event_end);
    %[B_pre,B_post] = pre_post(mms1_fgm_timedata_raw,mms1_fgm_bdata_raw,event_start,event_end);
    [V_pre,V_post] = pre_post(fpi_i_timedata,fpi_i_vdata,event_start,event_end);
    E_pre = -cross(V_pre,B_pre(1:3));
    E_post= (-cross(V_post,B_post(1:3)));
    
    %Shock Angles before and after event
    Theta_B_pre_n   = acosd(dot(B_pre,bowshock_n)/norm(B_pre))
    Theta_B_post_n  = acosd(dot(B_post,bowshock_n)/norm(B_post))
    Theta_BpreBpost = acosd(dot(B_pre,B_post)/(norm(B_pre)*norm(B_post)));
    
    
    
    %Current Sheet Speeds, Use V_sw before Event
    V_pre_ncs  = dot(V_pre,currentsheet_n);
    V_post_ncs = dot(V_post,currentsheet_n);
    V_ncs      = V_pre_ncs
    
    %Electric field Direction to CS
    E_pre_dot_cs_n  = dot(E_pre,currentsheet_n) %positive is towards
    E_post_dot_cs_n = (-1)*dot(E_post,currentsheet_n) %positive is towards
    
    %Angle between CS and BS
    Theta_csbs = acosd(dot(bowshock_n,currentsheet_n))
    
    %magnetic shear angle, if less than 30, probably a SHFA, no
    %discontinuity.
    Shear_angle = angle(B_pre,B_post)
    
    %Transversal Speed
    V_tr = V_ncs/((sind(Theta_csbs))^2)*(currentsheet_n-bowshock_n*cosd(Theta_csbs))
    mag_V_tr = norm(V_tr)
    
    %angle with solar wind velocity
    Theta_csswpre=acosd(dot(V_pre,currentsheet_n)/norm(V_pre));
    Theta_bsswpre=acosd(dot(V_pre,bowshock_n)/norm(V_pre));
    Theta_csswpost=acosd(dot(V_post,currentsheet_n)/norm(V_post));
    Theta_bsswpost=acosd(dot(V_post,bowshock_n)/norm(V_post));
    
    %Pre Vg
    V_g_pre = 2 * abs(dot(V_pre,bowshock_n) * sind(Theta_B_pre_n));
    V_tr_pre_over_V_g = abs(mag_V_tr/V_g_pre)
    
    %V_tr_pre_over_V_g_2 = abs(cosd(Theta_csswpre))/(2*cosd(Theta_bsswpre)*sind(Theta_B_pre_n)*sind(Theta_csbs))
    
    %Post Vg
    V_g_post = 2*abs(dot(V_post,bowshock_n) * sind(Theta_B_post_n));
    V_tr_post_over_V_g = abs(mag_V_tr/V_g_post)
    
    %V_tr_post_over_V_g_2 = abs(cosd(Theta_csswpost))/(2*cosd(Theta_bsswpost)*sind(Theta_B_post_n)*sind(Theta_csbs))
    
    
    %size
    Size = mag_V_tr*duration/6371.2
    
    Vn_1 = timing_n1.*timing_v1
    Vn_2 = timing_n2.*(timing_v2)
    
    Delta_V =(Vn_2-Vn_1)
    %HFA_exp_V = abs(dot(Delta_V,N_cs))
    %Delta_V=(392-325);
    HFA_exp_V=abs(dot(Delta_V,currentsheet_n))
    
    
    Event_age = Size/HFA_exp_V*6371.2
    Distance_traveled = Event_age*mag_V_tr/6371.2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% %Write to Excel Spreadsheet; Event Parameters Database
    cd '~/Library/Mobile Documents/com~apple~CloudDocs/Research/Analysis'
    databaseFilename = 'EventParameterDatabase.xlsx';
    
    % if Substructure == 0
    %    sub = 'no';
    % else
    %     sub = 'yes';
    % end
    
    T = table(Event_number,Substructure,convertCharsToStrings(event_start),convertCharsToStrings(event_end),r_sc,currentsheet_n, bowshock_n,Theta_csbs,Shear_angle,...
        convertCharsToStrings(leading_start),convertCharsToStrings(leading_end),MVA_n1,timing_n1,timing_v1,MVAB_timing_angle1,...
        convertCharsToStrings(trailing_start),convertCharsToStrings(trailing_end),MVA_n2,timing_n2,timing_v2,MVAB_timing_angle2,...
        B_pre,B_post,Theta_BpreBpost,V_pre,V_post,Theta_B_pre_n,Theta_B_post_n,E_pre_dot_cs_n,E_post_dot_cs_n,V_ncs,V_tr,mag_V_tr,...
        V_g_pre,V_tr_pre_over_V_g,V_g_post,V_tr_post_over_V_g,...
        Vn_1,Vn_2, Delta_V, Size, HFA_exp_V,...
        duration,Event_age,Distance_traveled);
    
    if Event_number == 0
        %Don't write to Excel, this is a practice run.
    else
        %Write to Excel on the next line, which corresponds to Event number.
        EventRow = strcat('A',num2str(Event_number));
        writetable(T,databaseFilename,'Sheet',1,'WriteVariableNames',false,'Range',EventRow) %Master sheet of all events
        
        
        %Substructure yes/no sheet division
        if Substructure == 1 %yes substructure
            writetable(T,databaseFilename,'Sheet',2,'WriteVariableNames',false,'Range',EventRow)
        elseif Substructure == 0 %no substructure
            writetable(T,databaseFilename,'Sheet',3,'WriteVariableNames',false,'Range',EventRow)
        end
        
    end
end

