function get_G_dist(requested_voltages)
    close all
    clearvars -except requested_voltages

    exp_data_dir='../../../DBs/memdiode_measurements/exp_data/M2';
    meas_files=dir(fullfile(exp_data_dir,'*.txt'));

    for i=2:size(meas_files,1)
        if ~isempty(strfind(meas_files(i,1).name,'_B')) && isempty(strfind(meas_files(i,1).name,'_B1.txt'))
            data=importdata(fullfile(exp_data_dir,meas_files(i,1).name));
            if exist('It_data_set_I','var')
                It_data_set_I=horzcat(It_data_set_I,data(:,2));
                It_data_set_t=horzcat(It_data_set_t,data(:,3));
            else
                It_data_set_I=data(:,2);
                It_data_set_t=data(:,3);                
            end

            if exist('Vt_data_set_V','var')
                Vt_data_set_V=horzcat(Vt_data_set_V,data(:,1));
                Vt_data_set_t=horzcat(Vt_data_set_t,data(:,3));
            else
                Vt_data_set_V=data(:,1);
                Vt_data_set_t=data(:,3);                
            end
        end

        if ~isempty(strfind(meas_files(i,1).name,'_R')) && isempty(strfind(meas_files(i,1).name,'_R1.txt'))
            data=importdata(fullfile(exp_data_dir,meas_files(i,1).name));
            if exist('It_data_reset_I','var')
                It_data_reset_I=horzcat(It_data_reset_I,data(:,2));
                It_data_reset_t=horzcat(It_data_reset_t,data(:,3));
            else
                It_data_reset_I=data(:,2);
                It_data_reset_t=data(:,3);                
            end

            if exist('Vt_data_reset_V','var')
                Vt_data_reset_V=horzcat(Vt_data_reset_V,data(:,1));
                Vt_data_reset_t=horzcat(Vt_data_reset_t,data(:,3));
            else
                Vt_data_reset_V=data(:,1);
                Vt_data_reset_t=data(:,3);                
            end
        end
        %fprintf('%d\n',i);
    end

    Rt_data_set=(Vt_data_set_V./It_data_set_I);
    Gt_data_set=(It_data_set_I./Vt_data_set_V);
    mean_It_data_set=[mean(It_data_set_t,2) mean(It_data_set_I,2)];
    mean_Vt_data_set=[mean(Vt_data_set_t,2) mean(Vt_data_set_V,2)];
    mean_Rt_data_set=[mean_Vt_data_set(:,2) mean_Vt_data_set(:,2)./mean_It_data_set(:,2)];  
    mean_Gt_data_set=[mean_Vt_data_set(:,2) mean_It_data_set(:,2)./mean_Vt_data_set(:,2)]; 
    std_It_data_set=[std(It_data_set_t,0,2) std(It_data_set_I,0,2)];
    std_Vt_data_set=[std(Vt_data_set_t) std(Vt_data_set_V)];
    std_Rt_data_set=[std(Rt_data_set')];
    std_Gt_data_set=[std(Gt_data_set')];

    mean_It_data_reset=[mean(It_data_reset_t,2) mean(It_data_reset_I,2)];
    mean_Vt_data_reset=[mean(Vt_data_reset_t,2) mean(Vt_data_reset_V,2)];

    mean_It_data_all=vertcat(mean_It_data_set,mean_It_data_reset+mean_It_data_set(end,:));
    mean_Vt_data_all=vertcat(mean_Vt_data_set,mean_Vt_data_reset+mean_Vt_data_set(end,:));

    dlmwrite('../../../models/memdiode/model_fit/mean_It_data_all.txt',mean_It_data_all,'\t');
    dlmwrite('../../../models/memdiode/model_fit/mean_Vt_data_all.txt',mean_Vt_data_all,'\t');
    
    %% Figure containing I-V for SET and RESET lin and log
    figure();
    subplot(1,2,1)
    for i=1:size(It_data_set_I,2)
        plot(Vt_data_set_V(:,i),It_data_set_I(:,i),'Color',[0.5 0.5 0.5],'LineWidth',0.75)
        if i==1
            hold on
        end
    end
    for i=1:size(It_data_reset_I,2)
        plot(Vt_data_reset_V(:,i),It_data_reset_I(:,i),'Color',[0.5 0.5 0.5],'LineWidth',0.75)
        if i==1
            hold on
        end
    end
    plot(mean_Vt_data_set(:,2),mean_It_data_set(:,2),'MarkerSize',4,'Color','red','LineWidth',1.25);
    plot(mean_Vt_data_reset(:,2),mean_It_data_reset(:,2),'MarkerSize',4,'Color','red','LineWidth',1.25);

    xlabel('Voltage [V]')
    ylabel('Mean Current [A]')    
    axis([-2 1.2 0 12e-3]);
    
    subplot(1,2,2)
    for i=1:size(It_data_set_I,2)
        semilogy(Vt_data_set_V(:,i),It_data_set_I(:,i),'Color',[0.5 0.5 0.5],'LineWidth',0.75)
        if i==1
            hold on
        end
    end
    for i=1:size(It_data_reset_I,2)
        semilogy(Vt_data_reset_V(:,i),abs(It_data_reset_I(:,i)),'Color',[0.5 0.5 0.5],'LineWidth',0.75)
        if i==1
            hold on
        end
    end
    %plot(mean_Vt_data_set(:,2),mean_It_data_set(:,2),'MarkerSize',4,'Color','red','LineWidth',1.25);
    %plot(mean_Vt_data_reset(:,2),mean_It_data_reset(:,2),'MarkerSize',4,'Color','red','LineWidth',1.25);

    xlabel('Voltage [V]')
    ylabel('Mean Current [A]')      
    axis([-2 1.2 2e-7 2e-2]);
    
    %% Figure containing G, R and I for SET
    figure();
    subplot(2,2,1)
    for i=1:size(It_data_set_I,2)
        plot(Vt_data_set_V(:,i),It_data_set_I(:,i),'Color',[0.5 0.5 0.5],'LineWidth',0.75)
        if i==1
            hold on
        end
    end
    errorbar(mean_Vt_data_set(:,2),mean_It_data_set(:,2),std_It_data_set(:,2),'o','MarkerSize',4,'Color','red','LineWidth',1.25);
    xlabel('Voltage [V]')
    ylabel('Mean Current [A]')

    subplot(2,2,2)
    for i=1:size(Rt_data_set,2)
        plot(Vt_data_set_V(:,i),Rt_data_set(:,i),'Color',[0.5 0.5 0.5],'LineWidth',0.75)
        if i==1
            hold on
        end
    end
    errorbar(mean_Rt_data_set(:,1),mean_Rt_data_set(:,2),std_Rt_data_set,'o','MarkerSize',4,'Color','red','LineWidth',1.25);
    xlabel('Voltage [V]')
    ylabel('Mean Resistance [\Omega]')

    subplot(2,2,3)
    for i=1:size(Gt_data_set,2)
        plot(Vt_data_set_V(:,i),Gt_data_set(:,i),'Color',[0.5 0.5 0.5],'LineWidth',0.75)
        if i==1
            hold on
        end
    end
    errorbar(mean_Gt_data_set(:,1),mean_Gt_data_set(:,2),std_Gt_data_set,'o','MarkerSize',4,'Color','red','LineWidth',1.25);
    xlabel('Voltage [V]')
    ylabel('Mean Conductance [S]')


    for i=1:length(requested_voltages)
        aux_v=abs(Vt_data_set_V(:,1)-requested_voltages(i));
        [min_aux_v_val,min_aux_v_index]=min(aux_v);
        if min_aux_v_val>0.0025
            sprintf('The requested voltage value does not exist in the Voltage sweep. The closet one will be considered (%.2g)',Vt_data_set_V(min_aux_v_index));
        end
        if ~exist('voltage_idx','var')
            voltage_idx=min_aux_v_index;
        else
            voltage_idx=vertcat(voltage_idx,min_aux_v_index);
        end
    end
    for i=1:length(voltage_idx)
        figure();
        histogram(Gt_data_set(voltage_idx(i),:));
        %if i==1
        %axis([min((min(Gt_data_set(:,1))))*0.8 max((max(Gt_data_set(2:end-1,1))))*1.2 -inf inf]);
        xlabel('Conductance [S]')
        %hold on
        %end
        fit_distro=fitdist(Gt_data_set(voltage_idx(i),:)','Normal');
        annotation('textbox',...
        [0.15 0.65 0.22 0.23],...
        'String',sprintf('Voltage: %.2g, %s distribution, \\mu=%.2g, \\sigma=%.2g', Vt_data_set_V(voltage_idx(i)), fit_distro.DistributionName, fit_distro.mu, fit_distro.sigma),...
        'FontSize',11,...
        'FontName','Arial',...
        'LineStyle','-',...
        'EdgeColor',[0 0 0],...
        'LineWidth',1.2,...
        'BackgroundColor',[1 1 1],...
        'Color',[0 0 0]);
    end
    %range2plot=linspace(min((min(Gt_data_set(:,1))))*0.8,max((max(Gt_data_set(2:end-1,1))))*1.2,1000);
    %y=pdf(fit_distro,range2plot);
    %plot(range2plot,y,'LineWidth',2)

    %diff_Rt=diff(mean_Vt_data_set(:,2))./diff(mean_It_data_set(:,2));

    % figure();
    % plot(mean_Vt_data_set(1:end-1,2),diff_Rt);
end