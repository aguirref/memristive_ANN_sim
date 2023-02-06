function processed_data=generate_plots(RAW_DATA,gen_Ht, gen_It, gen_In_t,gen_sys_crl,results_folder,varargin)

    simEngine=RAW_DATA.simEngine;
    omit_size_search=0;
    %% Optional values assignment
    if mod(nargin-6,2)
        display('A odd number of arguments have been introduced. Please revise the function call');
    else
        for i=1:nargin-6
            if strcmp(varargin{i},'networkSize')
                size_A=varargin{i+1}(1);
                size_B=varargin{i+1}(2);
                omit_size_search=1;
            end
        end
    end
    
    if omit_size_search==0
        size_A=0;
        size_B=0;
        i=1;
        while i<=size(RAW_DATA.variable_name_list,2)
            ii=1;        
            while ii<=size(RAW_DATA.variable_name_list,2)
                for iii=1:size(RAW_DATA.variable_name_list,2)
                    variable_name=RAW_DATA.variable_name_list{iii};
                    string_search=sprintf('V(mmed_network:memda%d-%db%d-%d:h)', i, ii, ii, i);
                    if strcmp(string_search,variable_name)
                        if i>size_A
                            size_A=i;
                        end
                        if ii>size_B
                            size_B=ii;
                        end
                        break
                    end
                end
                ii=ii+1;
            end
            i=i+1;
        end
    end
    
    size_matrix=[size_A size_B];
    num_neurons_out=size_matrix(2);
    root_dir='D:\UAB_RS_neural_networks';
    results_file='netlist10_results';    
    time=RAW_DATA.time_vect;
    
    %% Get H(t) data for each Memdiodo
    data_cols=1:1:size(RAW_DATA.variable_name_list,2);
    for i=1:size_A
        for ii=1:size_B
            for iii=1:size(data_cols,2)
                %Search the H cols
                variable_name=RAW_DATA.variable_name_list{data_cols(iii)};
                if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                    string_search=sprintf('v(xmmed_network.xmemda%d-%db%d-%d.h)', i, ii, ii, i);
                else
                    string_search=sprintf('v(mmed_network:memda%d-%db%d-%d:h)', i, ii, ii, i);
                end
                if strcmpi(string_search,variable_name)
                    if exist('variable_cols_H','var')
                        variable_cols_H=horzcat(variable_cols_H,data_cols(iii));
                    else
                        variable_cols_H=data_cols(iii);
                    end
                    %delete_H=iii;
                    if exist('delete_columns','var')==0
                        delete_columns=iii;
                    else
                        delete_columns=vertcat(delete_columns,iii);
                    end
                end
                %search I cols
                if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                    string_search=sprintf('i(xmmed_network.rcsa%d-%db%d-%d_1)', i, ii, ii, i);
                else
                    string_search=sprintf('i(mmed_network:rcsa%d-%db%d-%d_1)', i, ii, ii, i);
                end
                if strcmpi(string_search,variable_name) 
                    if exist('variable_cols_I','var')
                        variable_cols_I=horzcat(variable_cols_I,data_cols(iii));
                    else
                        variable_cols_I=data_cols(iii);
                    end
                    %delete_I=iii;
                    if exist('delete_columns','var')==0
                        delete_columns=iii;
                    else
                        delete_columns=vertcat(delete_columns,iii);
                    end
                end                
                %search In cols
                if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                    string_search=sprintf('i(rneuronb%d)', i);
                else
                    string_search=sprintf('i(rneuronb%d)', i);
                end
                if strcmpi(string_search,variable_name) 
                    if exist('variable_cols_In','var')
                        variable_cols_In=horzcat(variable_cols_In,data_cols(iii));
                    else
                        variable_cols_In=data_cols(iii);
                    end
                    %delete_In=iii;
                    if exist('delete_columns','var')==0
                        delete_columns=iii;
                    else
                        delete_columns=vertcat(delete_columns,iii);
                    end                    
                end                
                %search Vp cols
                if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                    string_search=sprintf('v(xmmed_network.a%d-%d_p)', i, ii);
                else
                    string_search=sprintf('v(mmed_network:a%d-%d_p)', i, ii);
                end
                if strcmpi(string_search,variable_name) 
                    if exist('variable_cols_Vp','var')
                        variable_cols_Vp=horzcat(variable_cols_Vp,data_cols(iii));
                    else
                        variable_cols_Vp=data_cols(iii);
                    end
                    %delete_Vp=iii;
                    if exist('delete_columns','var')==0
                        delete_columns=iii;
                    else
                        delete_columns=vertcat(delete_columns,iii);
                    end                   
                end                
                %search Vn cols
                if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                    string_search=sprintf('v(xmmed_network.b%d-%d_n)', ii, i);
                else
                    string_search=sprintf('v(mmed_network:b%d-%d_n)', ii, i);
                end
                if strcmpi(string_search,variable_name) 
                    if exist('variable_cols_Vn','var')
                        variable_cols_Vn=horzcat(variable_cols_Vn,data_cols(iii));
                    else
                        variable_cols_Vn=data_cols(iii);
                    end
                    %delete_Vn=iii;
                    if exist('delete_columns','var')==0
                        delete_columns=iii;
                    else
                        delete_columns=vertcat(delete_columns,iii);
                    end                    
                end     
%                 if strcmp('v(a1_pad_1)',variable_name)
%                     A1_PAD_1=RAW_DATA.variable_mat(iii,:,:);
%                 end
%                 if strcmp('v(b1_pad_1)',variable_name)
%                     B1_PAD_1=RAW_DATA.variable_mat(iii,:,:);
%                 end
%                 if strcmp('v(a2_pad_1)',variable_name)
%                     A2_PAD_1=RAW_DATA.variable_mat(iii,:,:);
%                 end
%                 if strcmp('v(b2_pad_1)',variable_name)
%                     B2_PAD_1=RAW_DATA.variable_mat(iii,:,:);
%                 end                
%                 if strcmp('v(vconta1-1)',variable_name)
%                     vconta1_1=RAW_DATA.variable_mat(iii,:,:);
%                 end
%                 if strcmp('v(vcontb1-1)',variable_name)
%                     vcontb1_1=RAW_DATA.variable_mat(iii,:,:);
%                 end
%                 if strcmp('v(vconta2-1)',variable_name)
%                     vconta2_1=RAW_DATA.variable_mat(iii,:,:);
%                 end
%                 if strcmp('v(vcontb2-1)',variable_name)
%                     vcontb2_1=RAW_DATA.variable_mat(iii,:,:);
%                 end 
%                 if strcmp('v(vconta1-1_out)',variable_name)
%                     vconta1_1_out=RAW_DATA.variable_mat(iii,:,:);
%                 end
%                 if strcmp('v(vcontb1-1_out)',variable_name)
%                     vcontb1_1_out=RAW_DATA.variable_mat(iii,:,:);
%                 end
%                 if strcmp('v(vconta2-1_out)',variable_name)
%                     vconta2_1_out=RAW_DATA.variable_mat(iii,:,:);
%                 end
%                 if strcmp('v(vcontb2-1_out)',variable_name)
%                     vcontb2_1_out=RAW_DATA.variable_mat(iii,:,:);
%                 end                 
%                 if strcmp('v(vpulsedsignal)',variable_name)
%                     vpulsedsignal=RAW_DATA.variable_mat(iii,:,:);
%                 end
%                 if strcmp('v(write_signal_a1_1)',variable_name)
%                     write_signal_a1_1=RAW_DATA.variable_mat(iii,:,:);
%                 end
%                 if strcmp('v(write_signal_b1_1)',variable_name)
%                     write_signal_b1_1=RAW_DATA.variable_mat(iii,:,:);
%                 end
%                 if strcmp('v(wr_user_n)',variable_name)
%                     wr_user_n=RAW_DATA.variable_mat(iii,:,:);
%                 end                
%                 if strcmp('v(eventsa)',variable_name)
%                     eventsa=RAW_DATA.variable_mat(iii,:,:);
%                 end                
%                 if strcmp('v(eventsb)',variable_name)
%                     eventsb=RAW_DATA.variable_mat(iii,:,:);
%                 end                
%                 if strcmp('v(eventsb_in_buf)',variable_name)
%                     eventsb_in_buf=RAW_DATA.variable_mat(iii,:,:);
%                 end    
                if strcmpi('v(write_ok!)',variable_name) && ~exist('write_ok','var')
                    write_ok=RAW_DATA.variable_mat(iii,:,:);
                    processed_data.WR_ok=write_ok;
                end  
                
                if strcmpi('v(rw)',variable_name) && ~exist('RW','var')
                    RW=RAW_DATA.variable_mat(iii,:,:);
                    processed_data.RW=RW;
                end   
                if strcmpi('v(comparator)',variable_name) && ~exist('comparator','var')
                    comparator=RAW_DATA.variable_mat(iii,:,:);
                    processed_data.comparator=comparator;
                end
%                 if strcmp('v(v_ref)',variable_name)
%                     v_ref=RAW_DATA.variable_mat(iii,:,:);
%                 end       
%                 if strcmp('v(clear)',variable_name)
%                     clear_signal=RAW_DATA.variable_mat(iii,:,:);
%                 end 
%                 if strcmp('v(xcountera.n004)',variable_name)
%                     xcountera_n004=RAW_DATA.variable_mat(iii,:,:);
%                 end       
%                 if strcmp('v(xcountera.output_2)',variable_name)
%                     xcountera_output_2=RAW_DATA.variable_mat(iii,:,:);
%                 end  
%                 if strcmp('v(xcounterb.n004)',variable_name)
%                     xcounterb_n004=RAW_DATA.variable_mat(iii,:,:);
%                 end       
%                 if strcmp('v(xcounterb.output_2)',variable_name)
%                     xcounterb_output_2=RAW_DATA.variable_mat(iii,:,:);
%                 end  
%                 if strcmp('v(event!)',variable_name)
%                     event=RAW_DATA.variable_mat(iii,:,:);
%                 end                  
                if strcmpi('v(sys_clk)',variable_name) && ~exist('sys_clk','var')
                    sys_clk=RAW_DATA.variable_mat(iii,:,:);
                    processed_data.sys_clk=sys_clk;
                end  
                if strcmpi('v(enable_write)',variable_name) && ~exist('enable_write','var')
                    enable_write=RAW_DATA.variable_mat(iii,:,:);
                    processed_data.enable_write=enable_write;
                end
                if strcmpi('v(eventsa)',variable_name) && ~exist('eventsa','var')
                    eventsa=RAW_DATA.variable_mat(iii,:,:);
                    processed_data.enable_write=eventsa;
                end
                if strcmpi('v(event!)',variable_name) && ~exist('event','var')
                    event=RAW_DATA.variable_mat(iii,:,:);
                    processed_data.event=event;
                end                
%                 if strcmp('v(sr_rw_sync)',variable_name)
%                     sr_rw_sync=RAW_DATA.variable_mat(iii,:,:);
%                 end
%                 if strcmp('v(event_out_q)',variable_name)
%                     event_out_q=RAW_DATA.variable_mat(iii,:,:);
%                 end
%                 if strcmp('v(event_out_q_c)',variable_name)
%                     event_out_q_C=RAW_DATA.variable_mat(iii,:,:);
%                 end
%                 if strcmp('v(clear_event)',variable_name)
%                     clear_event=RAW_DATA.variable_mat(iii,:,:);
%                 end
%                 if strcmp('v(event_out!)',variable_name)
%                     event_out=RAW_DATA.variable_mat(iii,:,:);
%                 end 
%                 if strcmp('v(event_in!)',variable_name)
%                     event_in=RAW_DATA.variable_mat(iii,:,:);
%                 end  
%                 if strcmp('v(write_pulse_ok!)',variable_name)
%                     write_pulse_ok=RAW_DATA.variable_mat(iii,:,:);
%                 end  
%                 if strcmp('v(write_pulse_ok_n!)',variable_name)
%                     write_pulse_ok_n=RAW_DATA.variable_mat(iii,:,:);
%                 end  
%                 if strcmp('v(write_okn!)',variable_name)
%                     write_okn=RAW_DATA.variable_mat(iii,:,:);
%                 end                    
            end
            delete_columns=sortrows(delete_columns,'descend');
            for iiii=1:size(delete_columns,1)
                data_cols(delete_columns(iiii))=[];
            end
            clear delete_columns
        end
    end
    
    DATA_H=reshape(RAW_DATA.variable_mat(variable_cols_H,:,:),[size_B, size_A, size(RAW_DATA.variable_mat,2), size(RAW_DATA.variable_mat,3)]);
    if size(DATA_H,4)>1
        for i=2:size(DATA_H,4)
            ERROR_H(:,:,i-1)=(DATA_H(:,:,end,i)-DATA_H(:,:,end,i-1))./DATA_H(:,:,end,1)*100;
        end
        processed_data.ERROR_H=ERROR_H;
    end
    processed_data.DATA_H=DATA_H;
    clear variable_cols
    
    %% Get I(t) data for each Memdiodo
    DATA_I=reshape(RAW_DATA.variable_mat(variable_cols_I,:,:),[size_B, size_A, size(RAW_DATA.variable_mat,2), size(RAW_DATA.variable_mat,3)]);
    if size(DATA_I,4)>1
        for i=2:size(DATA_I,4)
            ERROR_I(:,:,i-1)=(DATA_I(:,:,end,i)-DATA_I(:,:,end,i-1))./DATA_I(:,:,end,1)*100;
        end    
        processed_data.ERROR_I=ERROR_I;
    end
    processed_data.DATA_I=DATA_I;
    clear variable_cols

        
    
    %% Get I_neuron(t) data for each Memdiodo
    DATA_In=RAW_DATA.variable_mat(variable_cols_In,:,:);
    if size(DATA_In,3)>1
        for i=2:size(DATA_In,3)
            ERROR_In(:,i-1)=(DATA_In(:,end,i)-DATA_In(:,end,i-1))./DATA_In(:,end,1)*100;
        end 
        processed_data.ERROR_In=ERROR_In;
    end
    processed_data.DATA_In=DATA_In;
    clear variable_cols
    %% Get Vp(t) data for each Memdiodo
    DATA_Vp=reshape(RAW_DATA.variable_mat(variable_cols_Vp,:,:),[size_B, size_A, size(RAW_DATA.variable_mat,2), size(RAW_DATA.variable_mat,3)]);
    if size(DATA_Vp,4)>1
        for i=2:size(DATA_H,4)
            ERROR_Vp(:,:,i-1)=(DATA_Vp(:,:,end,i)-DATA_Vp(:,:,end,i-1))./DATA_Vp(:,:,end,1)*100;
        end
        processed_data.ERROR_Vp=ERROR_Vp;            
    end
    processed_data.DATA_Vp=DATA_Vp;    
    clear variable_cols
    
    %% Get Vn(t) data for each Memdiodo
    DATA_Vn=reshape(RAW_DATA.variable_mat(variable_cols_Vn,:,:),[size_B, size_A, size(RAW_DATA.variable_mat,2), size(RAW_DATA.variable_mat,3)]);
    if size(DATA_Vn,4)>1
        for i=2:size(DATA_H,4)
            ERROR_Vn(:,:,i-1)=(DATA_Vn(:,:,end,i)-DATA_Vn(:,:,end,i-1))./DATA_Vn(:,:,end,1)*100;
        end
        processed_data.ERROR_Vn=ERROR_Vn;            
    end
    processed_data.DATA_Vn=DATA_Vn;
    clear variable_cols

    %% plot generation of H(t) for each Rs value 
    if gen_Ht==1
            image_folder=fullfile(results_folder,'figures',strcat('plot_H_vs_pos_3D_',num2str(size_matrix(1)),'_by_',num2str(size_matrix(2))));

            if exist(image_folder,'dir')==0
            mkdir(image_folder);
        end
        p=numSubplots(size(DATA_H,4));
        figure();
        for ii=1:size(DATA_H,4)

            subplot(p(1),p(2),ii)

            bar3(squeeze(DATA_H(:,:,end,ii)));
            title(sprintf('H(t), time= %.6g',time(end)))
            xlabel('Layer A')
            ylabel('Layer B')
            zlabel('Memdiode H(t)')   
            if size(DATA_H,4)>1
                axis([0 size_matrix(1)+1 0 size_matrix(2)+1 0 max(max(max(max(DATA_H(:,:,:,end)))))*1.2]); 
            else
                axis([0 size_matrix(1)+1 0 size_matrix(2)+1 0 max(max(DATA_H(:,:,end)))*1.2]); 
            end
        end

        print (gcf,'-dpng',fullfile(image_folder,strcat('imagen_end.png')))
        savefig(fullfile(image_folder,strcat('imagen_end.fig')));

        %Plot generation of the normalized error in H(t)
        if size(DATA_H,4)>1
            p=numSubplots(size(ERROR_H,3));
            figure();
            for i=1:size(ERROR_H,3)

                subplot(p(1),p(2),i)

                bar3(abs(squeeze(ERROR_H(:,:,i))));
                title(sprintf('abs(Normalized ERROR of H(t_{end}))'))
                xlabel('Layer A')
                ylabel('Layer B')           
                max_z=max(max(max(abs(ERROR_H))))*1.2;
                if max_z==0
                    max_z=1e-3;
                end
                axis([0 size_matrix(1)+1 0 size_matrix(2)+1 0 max_z]);        
            end
        end
    end

    %% plot generation of I(t) for each Rs value 
    if gen_It==1
        image_folder=fullfile(results_folder,'figures',strcat('plot_I_vs_pos_3D_',num2str(size_matrix(1)),'_by_',num2str(size_matrix(2))));

        if exist(image_folder,'dir')==0
            mkdir(image_folder);
        end    
        p=numSubplots(size(DATA_I,4));
        figure();
        for ii=1:size(DATA_I,4)

            subplot(p(1),p(2),ii)

            bar3(squeeze(DATA_I(:,:,end,ii)));
            title(sprintf('I(t), time= %.6g',time(end)))
            xlabel('Layer A')
            ylabel('Layer B')
            zlabel('Memdiode current')
            if size(DATA_I,4)>1
                axis([0 size_matrix(1)+1 0 size_matrix(2)+1 0 max(max(max(max(DATA_I(:,:,:,end)))))*1.2]); 
            else
                axis([0 size_matrix(1)+1 0 size_matrix(2)+1 0 max(max(DATA_I(:,:,end)))*1.2]); 
            end
        end

        print (gcf,'-dpng',fullfile(image_folder,strcat('imagen_end.png')))
        savefig(fullfile(image_folder,strcat('imagen_end.fig')));

        %Plot generation of the normalized error in I(t)
        if size(DATA_I,4)>1
            p=numSubplots(size(ERROR_I,3));
            figure();
            for i=1:size(ERROR_I,3)

                subplot(p(1),p(2),i)

                bar3(abs(squeeze(ERROR_I(:,:,i))));
                title(sprintf('abs(Normalized ERROR of I(t_{end}))'))
                xlabel('Layer A')
                ylabel('Layer B')            
                max_z=max(max(max(abs(ERROR_I))))*1.2;
                if max_z==0
                    max_z=1e-3;
                end
                axis([0 size_matrix(1)+1 0 size_matrix(2)+1 0 max_z]);   
            end
        end
    end

    %% plot generation of In(t) for each Rs value 
    if gen_In_t==1
        image_folder=fullfile(results_folder,'figures',strcat('plot_In_vs_pos_',num2str(size_matrix(1)),'_by_',num2str(size_matrix(2))));

        if exist(image_folder,'dir')==0
            mkdir(image_folder);
        end    
        p=numSubplots(size(DATA_In,3));
        figure();
        for ii=1:size(DATA_In,3)

            subplot(p(1),p(2),ii)

            semilogy(squeeze(DATA_In(:,end,ii)));
            title(sprintf('I_{neuron}(t), time= %.6g',time(ii)))
            %ylabel('Layer A')
            xlabel('Layer B')
            ylabel('Output Neuron current')
            %axis([0 size_matrix(1)+1 0 size_matrix(2)+1 0 5e-5]);  
            if size(DATA_In,4)>1
                axis([0 size_matrix(2)+1 10^floor(log10(min(min(min(DATA_In(:,end,:)))))) 10^ceil(log10(max(max(max(DATA_In(:,end,:))))))]);
            else
                axis([0 size_matrix(2)+1 10^floor(log10(min(min(DATA_In(:,end))))) 10^ceil(log10(max(max(DATA_In(:,end)))))]);        
            end            
        end

        print (gcf,'-dpng',fullfile(image_folder,strcat('imagen_end.png')))
        savefig(fullfile(image_folder,strcat('imagen_end.fig')));
        
        %Plot generation of the normalized error in I(t)
        if size(DATA_In,3)>1
            p=numSubplots(size(ERROR_In,2));
            figure();
            for i=1:size(ERROR_In,2)

                subplot(p(1),p(2),i)

                bar(abs(squeeze(ERROR_In(:,i))));
                title(sprintf('abs(Normalized ERROR of In(t_{end}))'))
                xlabel('Layer A')
                ylabel('Layer B')
                max_z=max(max(max(abs(ERROR_In))))*1.2;
                if max_z==0
                    max_z=1e-3;
                end
                axis([0 size_matrix(2)+1 0 max_z]);
            end
        end
    end
    %% plot generation of Vp-n(t) for each Rs value 
    if gen_Ht==1
        image_folder=fullfile(results_folder,'figures',strcat('plot_Vpn_vs_pos_3D_',num2str(size_matrix(1)),'_by_',num2str(size_matrix(2))));

        if exist(image_folder,'dir')==0
            mkdir(image_folder);
        end
        p=numSubplots(size(DATA_Vp,4));
        figure();
        for ii=1:size(DATA_Vp,4)

            subplot(p(1),p(2),ii)

            bar3(squeeze(DATA_Vp(:,:,end,ii)-DATA_Vn(:,:,end,ii)));
            title(sprintf('V_{p-n}(t), time= %.6g',time(end)))
            xlabel('Layer A')
            ylabel('Layer B')
            zlabel('Memdiode App. Voltage(t)') 
            if size(DATA_I,4)>1
                axis([0 size_matrix(1)+1 0 size_matrix(2)+1 0 max(max(max(max(DATA_Vp(:,:,:,end)-DATA_Vn(:,:,:,end)))))*1.2]); 
            else
                axis([0 size_matrix(1)+1 0 size_matrix(2)+1 0 max(max(DATA_Vp(:,:,end)-DATA_Vn(:,:,end)))*1.2]); 
            end            
        end
    
        print (gcf,'-dpng',fullfile(image_folder,strcat('imagen_end.png')))
        savefig(fullfile(image_folder,strcat('imagen_end.fig')));
        
        %Plot generation of the normalized error in H(t)
        if size(DATA_Vp,4)>1
            p=numSubplots(size(ERROR_Vp,3));
            figure();
            for i=1:size(ERROR_Vp,3)

                subplot(p(1),p(2),i)

                bar3(abs(squeeze(ERROR_Vp(:,:,i))));
                title(sprintf('abs(Normalized ERROR of Memdiode App. Voltage(t))'))
                xlabel('Layer A')
                ylabel('Layer B')
                max_z=max(max(max(abs(ERROR_Vp))))*1.2;
                if max_z==0
                    max_z=1e-3;
                end
                axis([0 size_matrix(1)+1 0 size_matrix(2)+1 0 max_z]);
            end
        end
    end
    if gen_sys_crl==1
        for step_param_iteration=1:size(comparator,3)
            image_folder=fullfile(results_folder,'figures',strcat('plot_crl_vars_vs_time_',num2str(size_matrix(1)),'_by_',num2str(size_matrix(2)),num2str(step_param_iteration)));

            if exist(image_folder,'dir')==0
                mkdir(image_folder);
            end

            figure()
            subplot(5,1,1)
            plot(RAW_DATA.time_vect,squeeze(comparator(1,:,step_param_iteration)));
            if min(squeeze(comparator(1,:,step_param_iteration)))==max(squeeze(comparator(1,:,step_param_iteration))) && min(squeeze(comparator(1,:,step_param_iteration)))==0
                axis([-inf inf -inf inf]);
            else
                axis([-inf inf min(squeeze(comparator(1,:,step_param_iteration)))*0.9 max(squeeze(comparator(1,:,step_param_iteration)))*1.1]);
            end
            ylabel('Comparator Voltage [V]');
            xlabel('Time [Sec.]');


            subplot(5,1,2)
            plot(RAW_DATA.time_vect,squeeze(enable_write(1,:,step_param_iteration)))
            if min(squeeze(enable_write(1,:,step_param_iteration)))==max(squeeze(enable_write(1,:,step_param_iteration))) && min(squeeze(enable_write(1,:,step_param_iteration)))==0
                axis([-inf inf -inf inf]);
            else
                axis([-inf inf min(squeeze(enable_write(1,:,step_param_iteration)))*0.9 max(squeeze(enable_write(1,:,step_param_iteration)))*1.1]);
            end            
            ylabel('en. Write Voltage [V]');
            xlabel('Time [Sec.]');


            subplot(5,1,3)
            plot(RAW_DATA.time_vect,squeeze(write_ok(1,:,step_param_iteration)))
            if min(squeeze(write_ok(1,:,step_param_iteration)))==max(squeeze(write_ok(1,:,step_param_iteration))) && min(squeeze(write_ok(1,:,step_param_iteration)))==0
                axis([-inf inf -inf inf]);
            else
                axis([-inf inf min(squeeze(write_ok(1,:,step_param_iteration)))*0.9 max(squeeze(write_ok(1,:,step_param_iteration)))*1.1]);
            end
            ylabel('write_ok Voltage [V]');
            xlabel('Time [Sec.]');


            subplot(5,1,4)
            plot(RAW_DATA.time_vect,squeeze(RW(1,:,step_param_iteration)))
            if min(squeeze(RW(1,:,step_param_iteration)))==max(squeeze(RW(1,:,step_param_iteration))) && min(squeeze(RW(1,:,step_param_iteration)))==0
                axis([-inf inf -inf inf]);
            else
                axis([-inf inf min(squeeze(RW(1,:,step_param_iteration)))*0.9 max(squeeze(RW(1,:,step_param_iteration)))*1.1]);
            end
            ylabel('RW Voltage [V]');
            xlabel('Time [Sec.]');


            subplot(5,1,5)
            plot(RAW_DATA.time_vect,squeeze(sys_clk(1,:,step_param_iteration)))
            if min(squeeze(sys_clk(1,:,step_param_iteration)))==max(squeeze(sys_clk(1,:,step_param_iteration))) && min(squeeze(sys_clk(1,:,step_param_iteration)))==0
                axis([-inf inf -inf inf]);
            else
                axis([-inf inf min(squeeze(sys_clk(1,:,step_param_iteration)))*0.9 max(squeeze(sys_clk(1,:,step_param_iteration)))*1.1]);
            end            
            ylabel('sys_clk Voltage [V]');
            xlabel('Time [Sec.]');

            print (gcf,'-dpng',fullfile(image_folder,strcat('imagen_sys_crl.png')))
            savefig(fullfile(image_folder,strcat('imagen_sys_crl.fig')));        
        end
    end
end

