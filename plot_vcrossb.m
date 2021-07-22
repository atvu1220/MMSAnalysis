function [] = plot_vcrossb(date_start,date_end,fgm_timedata,fgm_bdata,fpi_timedata,fpi_vdata,num_plots,plot_order,flag)
    %plots the plasma parameter with the given start date, end date, datas,
    %number of plots and plot order
    
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    tstart = datenum(date_start,formatIn);
    tend = datenum(date_end,formatIn);
    
    %Find the start and end limits of the event in the data
    start_index = find(fpi_timedata > tstart, 1);
    end_index = find(fpi_timedata > tend, 1);
    
    %interpolate data since fgm brst has more data points than fpi
    fgm_bdata_interp=zeros(length(fpi_timedata),3);
    
    for i=1:3
        [x,index] = unique(fgm_timedata);
        fgm_bdata_interp(:,i) = interp1(x,fgm_bdata(index,i),fpi_timedata);
    end
    
    
    
    %crop data
    fpi_timedata = fpi_timedata(start_index:end_index,1);
    
    for i=1:3
        v(:,i)=(fpi_vdata(start_index:end_index,i));
        b(:,i) = fgm_bdata_interp(start_index:end_index,i);
    end
    
    v_cross_b = -cross(v,b,2);
    
    subplot(num_plots,1,plot_order);
    
    n_cs = tdnormal(date_start,date_end,fgm_timedata,fgm_bdata,'event');
    v_cross_b_dot_n_cs= dot(v_cross_b,repmat(n_cs,length(v_cross_b),1),2);
%     v_cross_b_dot_n_cs = acosd(v_cross_b_dot_n_cs/norm(v_cross_b_dot_n_cs));
    
    %Plotting the electric field in TD normal coordinates or all 4;
    %If the leading edge Edotn_td is negative, the electric field points
    %towards the discontintuiy.
    %if the trailing edge Edotn_td is positive, the electric field points
    %towards the discontinuity
    
    if strcmp(flag,'n_cs')
        plot(fpi_timedata,v_cross_b_dot_n_cs,'LineWidth',1)
        line([fpi_timedata(1),fpi_timedata(end)],[0,0],'Color','k','LineStyle','--')
        legend({'E\cdotn_{cs}'},'FontSize',10)
    else
        plot(fpi_timedata,[v_cross_b,v_cross_b_dot_n_cs],'LineWidth',1)
        legend({'E_x','E_y','E_z','E\cdotn_{cs}'},'FontSize',10)
    end
    
    
    
    legend('boxoff')
    legend('Location','eastoutside')
    colormap('winter');
    datetick
    xlim([fpi_timedata(1) fpi_timedata(end)])
    ylabel({'E';'[nT\cdotkm/s]'},'FontSize', 14)
    set(gca, 'XTickLabel', [],'XMinorTick','on','YMinorTick','on','linewidth',1.25)
    
    
    
    
    
    
    
    
end

