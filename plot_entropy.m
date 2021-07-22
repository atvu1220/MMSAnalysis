function [] = plot_entropy(date_start,date_end,time,n,temp_para,temp_perp,specie,num_plots,plot_order)
    %plots the general plasma parameter with the given start date, end date, datas,
    %number of plots and plot order
    
    formatIn='yyyy-mm-dd HH:MM:SS.FFF';
    tstart = datenum(date_start,formatIn);
    tend = datenum(date_end,formatIn);
    
    
    %Find the start and end limits of the event in the data
    start_index = find(time >= tstart, 1);
    end_index = find(time >= tend, 1);
    
    %crop data
    time = time(start_index:end_index,1);
    n = n(start_index:end_index,:);
    temp_para = temp_para(start_index:end_index,:);
    temp_perp = temp_perp(start_index:end_index,:);
    
    %Entropy calculation
    eVtoJ = 1.60218e-19;
    boltzmannJ = 1.38064852e-23;
    boltzmanneV = 8.61733e-5;
    if strcmp(specie,'i')
        mass = 1.6726219e-27;
        legendString = '\eta_i';
    elseif strcmp(specie,'e')
        mass = 9.10938356e-31;
        legendString = '\eta_e';
    end
    
    temperature = ((temp_para).^2 + (temp_perp).^2).^(1/2);
    %entropy = mass.*n.*log(n.*boltzmanneV.*temperature./ (mass.*n).^(5/3));
    entropy = (temperature .* boltzmanneV./(n).^(2/3)); %Units of K eV / cm^-2
    
    
    subplot(num_plots,1,plot_order);
    
    
    plot(time,entropy(:,:),'LineWidth',1)
    legend(legendString,'FontSize',14)
    legend('boxoff')
    legend('Location','eastoutside')
    colormap('winter');
    datetick
    xlim([time(1) time(end)])
    ylabel({'Entropy';'[eV/cm^{-2}]'},'FontSize', 14)
    set(gca,'Yscale','log')
    set(gca, 'XTickLabel', [],'XMinorTick','on','YMinorTick','on','linewidth',1.25)
    
    
end

