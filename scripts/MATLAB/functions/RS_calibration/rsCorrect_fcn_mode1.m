function [outputArg1,outputArg2] = rsCorrect_fcn_mode1(inputArg1,inputArg2)
    cont_it=1;
    iteration=1;
    iterations=[1];
    mean_rel_I1=1;
    mean_rel_I2=1;
    %W_init_pos_corr=W_init_pos;
    %W_init_neg_corr=W_init_neg;

    op_map_fresh=simRRAM_NW_dualPart(matrix_size,'simNetlist','yes',...
                          'Rs',0.001,...
                          'R_shunt',R_shunt,...
                          'Rcs',0.001 ,...
                          'simTime',25e-6,...
                          'runMode','-big',...
                          'simWR',write,...
                          'simEngine',simEngine,...
                          'genPlot','no',...
                          'minStep',1e-6,...
                          'model_ver',model_ver,...
                          'RUNLVL',3,...
                          'vProg',1.25,...
                          'rVoltage',read_voltage(rVoltage_i),...
                          'nProc',8,...
                          'useWData',h_data,...
                          'deg_folder',deg_folder,...
                          'mappingInfo',map_array,...
                          'W_matrix_pos',W_matrix_pos,...%W_init',W_init,...
                          'W_matrix_neg',W_matrix_neg,...
                          'W_init_pos',W_init_pos,...%W_init',W_init,...
                          'W_init_neg',W_init_neg,...
                          'memdiodeModel',memdiode_model,...
                          'dualSide',dualSide_connect,...
                          'inputMatrix',ones(matrix_size(1),1),...%,Vin_mean',...
                          'outputFMT','binary',...
                          'nPart',partitions_N,...
                          'analysis','o_point',...
                          'simDir',sim_folder);

    I_pos=op_map_fresh.I_pos;                      
    I_neg=op_map_fresh.I_neg;
    Vpn_pos=op_map_fresh.Vpn_pos;
    Vpn_neg=op_map_fresh.Vpn_neg;

    save(fullfile('..','mat_files',strcat('NW_',num2str(matrix_size(1)),'by',num2str(matrix_size(2)),'_closed_loop_ideal_op.mat')),'I_pos','I_neg','Vpn_pos','Vpn_neg');

    while abs(mean_rel_I1)>0.1 || abs(mean_rel_I2)>0.1
        if iteration>25
            scale=scale-0.1;
            iteration=1;
            iterations=[1];
            mean_rel_I1=1;
            mean_rel_I2=1;
            clear min_rel_I1_vect min_rel_I2_vect mean_rel_I1_vect mean_rel_I2_vect
        end
        op_map=simRRAM_NW_dualPart(matrix_size,'simNetlist','yes',...
                                      'Rs',series_resistance(Rs_i),...
                                      'R_shunt',R_shunt,...
                                      'Rcs',R_cs(R_cs_i) ,...
                                      'simTime',25e-6,...
                                      'runMode','-big',...
                                      'simWR',write,...
                                      'simEngine',simEngine,...
                                      'genPlot','no',...
                                      'minStep',1e-6,...
                                      'model_ver',model_ver,...
                                      'RUNLVL',3,...
                                      'vProg',1.25,...
                                      'rVoltage',read_voltage(rVoltage_i),...
                                      'nProc',8,...
                                      'useWData',h_data,...
                                      'deg_folder',deg_folder,...
                                      'W_matrix_pos',W_matrix_pos,...%W_init',W_init,...
                                      'W_matrix_neg',W_matrix_neg,...
                                      'W_init_pos',W_init_pos_corr,...%W_init',W_init,...
                                      'W_init_neg',W_init_neg_corr,...
                                      'memdiodeModel',memdiode_model,...
                                      'dualSide',dualSide_connect,...
                                      'mappingInfo',map_array,...
                                      'inputMatrix',ones(matrix_size(1),1),...%Vin_mean',...
                                      'outputFMT','binary',...
                                      'nPart',partitions_N,...
                                      'correctsRs',correct_rs,...
                                      'analysis','o_point');

        str_to_display=op_map.cmd_history;

    %                             rel_Vpn1=(op_map.Vpn_pos)./(Vpn_pos);                          
    %                             rel_Vpn2=(op_map.Vpn_neg)./(Vpn_neg);                          
    %                             rel_I1=(op_map.I_pos)./(scale*I_pos);                          
    %                             rel_I2=(op_map.I_neg)./(scale*I_neg);      

        rel_Vpn1=(op_map.Vpn_pos-Vpn_pos)./(Vpn_pos);                          
        rel_Vpn2=(op_map.Vpn_neg-Vpn_neg)./(Vpn_neg);                          
        rel_I1=(op_map.I_pos-scale*I_pos)./(scale*I_pos);                          
        rel_I2=(op_map.I_neg-scale*I_neg)./(scale*I_neg);      

        %W_init_pos=W_init;
        %W_init_pos(W_init_pos<0)=0;
        corr_rel_I1=rel_I1(:,:,1);
    %                             corr_rel_I1(corr_rel_I1<-0.99)=-0.99;
    %                             corr_rel_I1(corr_rel_I1>0.9)=0.9;
        W_init_pos_corr=10.^(log10(W_init_pos_corr)+log10(1-(corr_rel_I1)));
        W_init_pos_corr(W_init_pos_corr>1)=1;
        W_init_pos_corr(W_init_pos_corr<0.01)=0.01;

        %W_init_neg=W_init;
        %W_init_neg(W_init_neg>0)=0;
        corr_rel_I2=rel_I2(:,:,1);
    %                             corr_rel_I2(corr_rel_I2<-0.99)=-0.99;
    %                             corr_rel_I2(corr_rel_I2>0.9)=0.9;
        W_init_neg_corr=10.^(log10(W_init_neg_corr)+log10(1-(corr_rel_I2)));
        W_init_neg_corr(W_init_neg_corr>1)=1;
        W_init_neg_corr(W_init_neg_corr<0.01)=0.01;
        %W_init=W_init_pos+W_init_neg;

        [W_pos_rows, W_pos_cols] = ind2sub(matrix_size,find(W_init_pos>0));
        [W_neg_rows, W_neg_cols] = ind2sub(matrix_size,find(W_init_neg>0));

        max_rel_I1=max(max(rel_I1(:,:,1)));
        max_rel_I2=max(max(rel_I2(:,:,1)));

        mean_rel_I1=mean(mean(rel_I1(:,:,1)));
        mean_rel_I2=mean(mean(rel_I2(:,:,1)));    

        min_rel_I1=min(min(rel_I1(:,:,1)));
        min_rel_I2=min(min(rel_I2(:,:,1)));

        if exist('min_rel_I1_vect','var')
            min_rel_I1_vect=horzcat(min_rel_I1_vect, min_rel_I1);
        else
            min_rel_I1_vect=min_rel_I1;
        end
        if exist('mean_rel_I1_vect','var')
            mean_rel_I1_vect=horzcat(mean_rel_I1_vect, mean_rel_I1);
        else
            mean_rel_I1_vect=mean_rel_I1;
        end
        if exist('min_rel_I2_vect','var')
            min_rel_I2_vect=horzcat(min_rel_I2_vect, min_rel_I2);
        else
            min_rel_I2_vect=min_rel_I2;
        end
        if exist('mean_rel_I2_vect','var')
            mean_rel_I2_vect=horzcat(mean_rel_I2_vect, mean_rel_I2);
        else
            mean_rel_I2_vect=mean_rel_I2;
        end

        max_rel_I1_nZ=max(max(rel_I1(W_pos_rows, W_pos_cols,1)));
        max_rel_I2_nZ=max(max(rel_I2(W_neg_rows, W_neg_cols,1)));

        mean_rel_I1_nZ=mean(mean(rel_I1(W_pos_rows, W_pos_cols,1)));
        mean_rel_I2_nZ=mean(mean(rel_I2(W_neg_rows, W_neg_cols,1)));      

        min_rel_I1_nZ=min(min(rel_I1(W_pos_rows, W_pos_cols,1)));
        min_rel_I2_nZ=min(min(rel_I2(W_neg_rows, W_neg_cols,1)));

    %                         err_I1_max=max([abs(max_rel_I1)*100;
    %                         err_I1_min=max([abs(max_rel_I1-1) abs(min_rel_I1-1)])*100;
    % 
    %                         err_I2=max([abs(max_rel_I2-1) abs(min_rel_I2-1)])*100;

        save(fullfile(op_map.folder,'op_map.mat'),'op_map');

        if iteration==1
            h=figure('Name',sprintf('Relative Error Map %s, %s, Rs=%.2f, Rcs=%.2f, N=%d, V_read=%f, N_partitions=%d', connections, correct_folder, series_resistance(Rs_i), R_cs(R_cs_i), number_of_images(N_i), read_voltage(rVoltage_i), partitions_N));
        else
            figure(h)
        end
        subplot(2,2,1)
        imagesc(rel_Vpn1(:,:,1));
        colorbar
        caxis([-1.05 1.05])
        title('Vpn_{pos}')

        subplot(2,2,2)
        imagesc(rel_Vpn2(:,:,1));
        colorbar
        caxis([-1.05 1.05])
        title('Vpn_{neg}')

        subplot(2,2,3)
        imagesc(rel_I1(:,:,1));
        colorbar
        caxis([-1.05 1.05])
        title('I_{pos}')

        subplot(2,2,4)
        imagesc(rel_I2(:,:,1));
        colorbar
        caxis([-1.05 1.05]) 
        title('I_{neg}')

        if iteration==1
            h_error=figure('Name',sprintf('Error %s, %s, Rs=%.2f, Rcs=%.2f, N=%d, V_read=%f, N_partitions=%d', connections, correct_folder, series_resistance(Rs_i), R_cs(R_cs_i), number_of_images(N_i), read_voltage(rVoltage_i), partitions_N));
        else
            figure(h_error)
        end
        plot(iterations, min_rel_I1_vect)
        axis([-inf inf -1.1 1.1]);
        xlabel('Iterations')
        hold on
        plot(iterations, mean_rel_I1_vect)
        plot(iterations, min_rel_I2_vect)
        plot(iterations, mean_rel_I2_vect)
        hold off
        legend({'min. re. I_{pos}','Avg. rel. I_{pos}','min. rel. I_{neg}','Avg. rel. I_{neg}'},'Location','Northwest');
    %                         cont_it=input('Stop iteration?[y/N]: ','s');
    %                         if isempty(cont_it)==0
    %                             if(cont_it=='y' || cont_it=='Y')
    %                             cont_it=0;
    %                             end
    %                         else
    %                             cont_it=1;
    %                         end     
        iteration=iteration+1;
        iterations=horzcat(iterations,iteration);
    end
    clear min_rel_I1_vect min_rel_I2_vect mean_rel_I1_vect mean_rel_I2_vect
end

