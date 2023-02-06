function processed_data=generate_plots_dual(RAW_DATA, partitions_N, gen_Ht, gen_It, gen_In_t,gen_V, gen_sys_crl,results_folder,varargin)

    simEngine=RAW_DATA.simEngine;
    omit_size_search=0;
    pos_neg={'pos' 'neg'};
    
    %% Optional values assignment
    if mod(nargin-8,2)
        display('A odd number of arguments have been introduced. Please revise the function call');
    else
        for i=1:nargin-8
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
    for pos_neg_i=1:2
        for partitions_i=1:partitions_N
            for i=1:size_A/partitions_N
                for ii=1:size_B
                    for iii=1:size(data_cols,2)
                        %Search the H cols
                        variable_name=RAW_DATA.variable_name_list{data_cols(iii)};
                        if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                            string_search=sprintf('v(xmmed_network_P%d_%s.xmemda%d-%db%d-%d.h)', partitions_i, pos_neg{pos_neg_i}, i, ii, ii, i);
                        else
                            string_search=sprintf('v(mmed_network_P%d_%s:memda%d-%db%d-%d:h)', partitions_i, pos_neg{pos_neg_i}, i, ii, ii, i);
                        end
                        if strcmpi(string_search,variable_name)
                            if exist(sprintf('variable_cols_H_%s',pos_neg{pos_neg_i}),'var')
                                eval(sprintf('variable_cols_H_%s=horzcat(variable_cols_H_%s,data_cols(iii));',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
                            else
                                eval(sprintf('variable_cols_H_%s=data_cols(iii);',pos_neg{pos_neg_i}));
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
                            string_search=sprintf('i(xmmed_network_P%d_%s.rcsa%d-%db%d-%d_1)', partitions_i, pos_neg{pos_neg_i}, i, ii, ii, i);
                        else
                            string_search=sprintf('i(mmed_network_P%d_%s:rcsa%d-%db%d-%d_1)', partitions_i, pos_neg{pos_neg_i}, i, ii, ii, i);
                        end
                        if strcmpi(string_search,variable_name) 
                            if exist(sprintf('variable_cols_I_%s',pos_neg{pos_neg_i}),'var')
                                eval(sprintf('variable_cols_I_%s=horzcat(variable_cols_I_%s,data_cols(iii));',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
                            else
                                eval(sprintf('variable_cols_I_%s=data_cols(iii);', pos_neg{pos_neg_i}));
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
                            string_search=sprintf('i(rneuronb%d_P%d_%s)', i, partitions_i, pos_neg{pos_neg_i});
                        else
                            string_search=sprintf('i(rneuronb%d_P%d_%s)', i, partitions_i, pos_neg{pos_neg_i});
                        end
                        if strcmpi(string_search,variable_name) 
                            if exist(sprintf('variable_cols_In_%s',pos_neg{pos_neg_i}),'var')
                                eval(sprintf('variable_cols_In_%s=horzcat(variable_cols_In_%s,data_cols(iii));',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
                            else
                                eval(sprintf('variable_cols_In_%s=data_cols(iii);',pos_neg{pos_neg_i}));
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
                            string_search=sprintf('v(xmmed_network_P%d_%s.a%d-%d_p)', partitions_i, pos_neg{pos_neg_i}, i, ii);
                        else
                            string_search=sprintf('v(mmed_network_P%d_%s:a%d-%d_p)', partitions_i, pos_neg{pos_neg_i}, i, ii);
                        end
                        if strcmpi(string_search,variable_name) 
                            if exist(sprintf('variable_cols_Vp_%s', pos_neg{pos_neg_i}),'var')
                                eval(sprintf('variable_cols_Vp_%s=horzcat(variable_cols_Vp_%s,data_cols(iii));', pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
                            else
                                eval(sprintf('variable_cols_Vp_%s=data_cols(iii);',pos_neg{pos_neg_i}));
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
                            string_search=sprintf('v(xmmed_network_P%d_%s.b%d-%d_n)', partitions_i, pos_neg{pos_neg_i}, ii, i);
                        else
                            string_search=sprintf('v(mmed_network_P%d_%s:b%d-%d_n)', partitions_i, pos_neg{pos_neg_i}, ii, i);
                        end
                        if strcmpi(string_search,variable_name) 
                            if exist(sprintf('variable_cols_Vn_%s', pos_neg{pos_neg_i}),'var')
                                eval(sprintf('variable_cols_Vn_%s=horzcat(variable_cols_Vn_%s,data_cols(iii));',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
                            else
                                eval(sprintf('variable_cols_Vn_%s=data_cols(iii);',pos_neg{pos_neg_i}));
                            end
                            %delete_Vn=iii;
                            if exist('delete_columns','var')==0
                                delete_columns=iii;
                            else
                                delete_columns=vertcat(delete_columns,iii);
                            end                    
                        end     
  
                        if strcmpi('v(write_ok!)',variable_name) && ~exist('write_ok','var')
                            write_ok=RAW_DATA.variable_mat(data_cols(iii),:,:);
                            processed_data.WR_ok=write_ok;
                        end  

                        if strcmpi('v(rw)',variable_name) && ~exist('RW','var')
                            RW=RAW_DATA.variable_mat(data_cols(iii),:,:);
                            processed_data.RW=RW;
                        end   
                        if strcmpi(sprintf('v(comparator_P%d_%s)', partitions_i, pos_neg{pos_neg_i}),variable_name) && ~exist(sprintf('comparator_P%d_%s', partitions_i, pos_neg{pos_neg_i}),'var')
                            eval(sprintf('comparator_P%d_%s=RAW_DATA.variable_mat(data_cols(iii),:,:);', partitions_i, pos_neg{pos_neg_i}));
                            eval(sprintf('processed_data.comparator_P%d_%s=comparator_P%d_%s;',partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i}));
                        end
                 
                        if strcmpi('v(sys_clk)',variable_name) && ~exist('sys_clk','var')
                            sys_clk=RAW_DATA.variable_mat(data_cols(iii),:,:);
                            processed_data.sys_clk=sys_clk;
                        end  
                        if strcmpi(sprintf('v(enable_write_P%d_%s)',partitions_i, pos_neg{pos_neg_i}),variable_name) && ~exist(sprintf('enable_write_P%d_%s', partitions_i, pos_neg{pos_neg_i}),'var')
                            eval(sprintf('enable_write_P%d_%s=RAW_DATA.variable_mat(data_cols(iii),:,:);', partitions_i, pos_neg{pos_neg_i}));
                            eval(sprintf('processed_data.enable_write_P%d_%s=enable_write_P%d_%s;', partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i}));
                        end
                        if strcmpi('v(eventsa)',variable_name) && ~exist('eventsa','var')
                            eventsa=RAW_DATA.variable_mat(data_cols(iii),:,:);
                            processed_data.enable_write=eventsa;
                        end
                        if strcmpi(sprintf('v(event_P%d_%s!)', partitions_i, pos_neg{pos_neg_i}),variable_name) && ~exist(sprintf('event_P%d_%s', partitions_i, pos_neg{pos_neg_i}),'var')
                            eval(sprintf('event_P%d_%s=RAW_DATA.variable_mat(data_cols(iii),:,:);', partitions_i, pos_neg{pos_neg_i}));
                            eval(sprintf('processed_data.event_P%d_%s=event_P%d_%s;',partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i}));
                        end                
                  
                    end
                    delete_columns=sortrows(delete_columns,'descend');
                    for iiii=1:size(delete_columns,1)
                        data_cols(delete_columns(iiii))=[];
                    end
                    clear delete_columns
                end
            end
        end
    end
    
    for pos_neg_i=1:2
        eval(sprintf('DATA_H_%s=reshape(RAW_DATA.variable_mat(variable_cols_H_%s,:,:),[size_B, size_A, size(RAW_DATA.variable_mat,2), size(RAW_DATA.variable_mat,3)]);',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
        if size(DATA_H_pos,4)>1
            for i=2:size(DATA_H_pos,4)
                eval(sprintf('ERROR_H_%s(:,:,i-1)=(DATA_H_%s(:,:,end,i)-DATA_H_%s(:,:,end,i-1))./DATA_H_%s(:,:,end,1)*100;',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}, pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
            end
            eval(sprintf('processed_data.ERROR_H_%s=ERROR_H_%s;',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
        end
        eval(sprintf('processed_data.DATA_H_%s=DATA_H_%s;',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
        clear variable_cols

        %% Get I(t) data for each Memdiodo
        eval(sprintf('DATA_I_%s=reshape(RAW_DATA.variable_mat(variable_cols_I_%s,:,:),[size_B, size_A, size(RAW_DATA.variable_mat,2), size(RAW_DATA.variable_mat,3)]);', pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
        if size(DATA_I_pos,4)>1
            for i=2:size(DATA_I_pos,4)
                eval(sprintf('ERROR_I_%s(:,:,i-1)=(DATA_I_%s(:,:,end,i)-DATA_I_%s(:,:,end,i-1))./DATA_I_%s(:,:,end,1)*100;',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}, pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
            end    
            eval(sprintf('processed_data.ERROR_I_%s=ERROR_I_%s;',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
        end
        eval(sprintf('processed_data.DATA_I_%s=DATA_I_%s;',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
        clear variable_cols

        %% Get I_neuron(t) data for each Memdiodo
        eval(sprintf('DATA_In_%s=RAW_DATA.variable_mat(variable_cols_In_%s,:,:);',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
        if size(DATA_In_pos,3)>1
            for i=2:size(DATA_In_pos,3)
                eval(sprintf('ERROR_In_%s(:,i-1)=(DATA_In_%s(:,end,i)-DATA_In_%s(:,end,i-1))./DATA_In_%s(:,end,1)*100;',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}, pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
            end 
            eval(sprintf('processed_data.ERROR_In_%s=ERROR_In_%s;',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
        end
        eval(sprintf('processed_data.DATA_In_%s=DATA_In_%s;',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
        clear variable_cols
        
        %% Get Vp(t) data for each Memdiodo
        eval(sprintf('DATA_Vp_%s=reshape(RAW_DATA.variable_mat(variable_cols_Vp_%s,:,:),[size_B, size_A, size(RAW_DATA.variable_mat,2), size(RAW_DATA.variable_mat,3)]);',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
        if size(DATA_Vp_pos,4)>1
            for i=2:size(DATA_H_pos,4)
                eval(sprintf('ERROR_Vp_%s(:,:,i-1)=(DATA_Vp_%s(:,:,end,i)-DATA_Vp_%s(:,:,end,i-1))./DATA_Vp_%s(:,:,end,1)*100;',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}, pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
            end
            eval(sprintf('processed_data.ERROR_Vp_%s=ERROR_Vp_%s;',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));            
        end
        eval(sprintf('processed_data.DATA_Vp_%s=DATA_Vp_%s;',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));    
        clear variable_cols

        %% Get Vn(t) data for each Memdiodo
        eval(sprintf('DATA_Vn_%s=reshape(RAW_DATA.variable_mat(variable_cols_Vn_%s,:,:),[size_B, size_A, size(RAW_DATA.variable_mat,2), size(RAW_DATA.variable_mat,3)]);',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
        if size(DATA_Vn_pos,4)>1
            for i=2:size(DATA_H_pos,4)
                eval(sprintf('ERROR_Vn_%s(:,:,i-1)=(DATA_Vn_%s(:,:,end,i)-DATA_Vn_%s(:,:,end,i-1))./DATA_Vn_%s(:,:,end,1)*100;',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}, pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
            end
            eval(sprintf('processed_data.ERROR_Vn_%s=ERROR_Vn_%s;',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));            
        end
        eval(sprintf('processed_data.DATA_Vn_%s=DATA_Vn_%s;',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));    
        clear variable_cols
    end
    
    for pos_neg_i=1:2
    %% plot generation of H(t) for each Rs value 
        if gen_Ht==1
                image_folder=fullfile(results_folder,'figures',strcat('plot_H_vs_pos_3D_',num2str(size_matrix(1)),'_by_',num2str(size_matrix(2))));

                if exist(image_folder,'dir')==0
                mkdir(image_folder);
            end
            p=numSubplots(size(DATA_H_pos,4));
            figure();
            for ii=1:size(DATA_H_pos,4)

                subplot(p(1),p(2),ii)

                eval(sprintf('bar3(squeeze(DATA_H_%s(:,:,end,ii)));', pos_neg{pos_neg_i}));
                title(sprintf('H_%s(t), time= %.6g', pos_neg{pos_neg_i}, time(end)))
                xlabel('Layer A')
                ylabel('Layer B')
                zlabel('Memdiode H(t)')   
                if size(DATA_H_pos,4)>1
                    eval(sprintf('axis([0 size_matrix(1)+1 0 size_matrix(2)+1 0 max(max(max(max(DATA_H_%s(:,:,:,end)))))*1.2]);',pos_neg{pos_neg_i}));
                else
                    eval(sprintf('axis([0 size_matrix(1)+1 0 size_matrix(2)+1 0 max(max(DATA_H_%s(:,:,end)))*1.2]);',pos_neg{pos_neg_i}));
                end
            end

            print (gcf,'-dpng',fullfile(image_folder,sprintf('imagen_H_%s_end.png', pos_neg{pos_neg_i})))
            savefig(fullfile(image_folder,sprintf('imagen_H_%s_end.fig', pos_neg{pos_neg_i})));

            %Plot generation of the normalized error in H(t)
            if size(DATA_H_pos,4)>1
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
            p=numSubplots(size(DATA_I_pos,4));
            figure();
            for ii=1:size(DATA_I_pos,4)

                subplot(p(1),p(2),ii)

                eval(sprintf('bar3(squeeze(DATA_I_%s(:,:,end,ii)));', pos_neg{pos_neg_i}));
                title(sprintf('I_%s(t), time= %.6g', pos_neg{pos_neg_i}, time(end)));
                xlabel('Layer A')
                ylabel('Layer B')
                zlabel('Memdiode current')
                if size(DATA_I_pos,4)>1
                    eval(sprintf('axis([0 size_matrix(1)+1 0 size_matrix(2)+1 0 max(max(max(max(DATA_I_%s(:,:,:,end)))))*1.2]);',pos_neg{pos_neg_i}));
                else
                    eval(sprintf('axis([0 size_matrix(1)+1 0 size_matrix(2)+1 0 max(max(DATA_I_%s(:,:,end)))*1.2]);',pos_neg{pos_neg_i}));
                end
            end

            print (gcf,'-dpng',fullfile(image_folder,sprintf('imagen_I_%s_end.png', pos_neg{pos_neg_i})))
            savefig(fullfile(image_folder,sprintf('imagen_I_%s_end.fig', pos_neg{pos_neg_i})));

            %Plot generation of the normalized error in I(t)
            if size(DATA_I_pos,4)>1
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
            p=numSubplots(size(DATA_In_pos,3));
            figure();
            for ii=1:size(DATA_In_pos,3)

                subplot(p(1),p(2),ii)

                eval(sprintf('semilogy(squeeze(DATA_In_%s(:,end,ii)));',pos_neg{pos_neg_i}));
                title(sprintf('I_{neuron}%d(t), time= %.6g', pos_neg{pos_neg_i}, time(ii)))
                %ylabel('Layer A')
                xlabel('Layer B')
                ylabel('Output Neuron current')
                %axis([0 size_matrix(1)+1 0 size_matrix(2)+1 0 5e-5]);  
                if size(DATA_In_pos,4)>1
                    eval(sprintf('axis([0 size_matrix(2)+1 10^floor(log10(min(min(min(DATA_In_%s(:,end,:)))))) 10^ceil(log10(max(max(max(DATA_In_%s(:,end,:))))))]);',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));
                else
                    eval(sprintf('axis([0 size_matrix(2)+1 10^floor(log10(min(min(DATA_In_%s(:,end))))) 10^ceil(log10(max(max(DATA_In_%s(:,end)))))]);',pos_neg{pos_neg_i}, pos_neg{pos_neg_i}));       
                end            
            end

            print (gcf,'-dpng',fullfile(image_folder,sprintf('imagen_In_%s_end.png',pos_neg{pos_neg_i})))
            savefig(fullfile(image_folder,sprintf('imagen_In_%s_end.fig', pos_neg{pos_neg_i})));

            %Plot generation of the normalized error in I(t)
            if size(DATA_In_pos,3)>1
                p=numSubplots(size(ERROR_In_pos,2));
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
        if gen_V==1
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
            for partitions_i=1:partitions_N
                image_folder=fullfile(results_folder,'figures',strcat('plot_crl_vars_vs_time_',num2str(size_matrix(1)),'_by_',num2str(size_matrix(2))));

                if exist(image_folder,'dir')==0
                    mkdir(image_folder);
                end

                figure()
                subplot(5,1,1)
                eval(sprintf('plot(RAW_DATA.time_vect,squeeze(comparator_P%d_%s(1,:)));',partitions_i, pos_neg{pos_neg_i}));
                if eval(sprintf('min(squeeze(comparator_P%d_%s(1,:)))==max(squeeze(comparator_P%d_%s(1,:))) && min(squeeze(comparator_P%d_%s(1,:)))==0',partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i}))
                    axis([-inf inf -inf inf]);
                else
                    eval(sprintf('axis([-inf inf min(squeeze(comparator_P%d_%s(1,:)))*0.9 max(squeeze(comparator_P%d_%s(1,:)))*1.1]);',partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i}));
                end
                ylabel(sprintf('comparator_P%d_%s Voltage [V]', partitions_i, pos_neg{pos_neg_i}));
                xlabel('Time [Sec.]');


                subplot(5,1,2)
                eval(sprintf('plot(RAW_DATA.time_vect,squeeze(enable_write_P%d_%s(1,:)));',partitions_i, pos_neg{pos_neg_i}));
                if eval(sprintf('min(squeeze(enable_write_P%d_%s(1,:)))==max(squeeze(enable_write_P%d_%s(1,:))) && min(squeeze(enable_write_P%d_%s(1,:)))==0',partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i}))
                    axis([-inf inf -inf inf]);
                else
                    eval(sprintf('axis([-inf inf min(squeeze(enable_write_P%d_%s(1,:)))*0.9 max(squeeze(enable_write_P%d_%s(1,:)))*1.1]);',partitions_i, pos_neg{pos_neg_i}, partitions_i, pos_neg{pos_neg_i}));
                end     
                ylabel(sprintf('enable_write_P%d_%s Voltage [V]',partitions_i, pos_neg{pos_neg_i}));
                xlabel('Time [Sec.]');


                subplot(5,1,3)
                plot(RAW_DATA.time_vect,squeeze(write_ok(1,:)))
                if min(squeeze(write_ok(1,:)))==max(squeeze(write_ok(1,:))) && min(squeeze(write_ok(1,:)))==0
                    axis([-inf inf -inf inf]);
                else
                    axis([-inf inf min(squeeze(write_ok(1,:)))*0.9 max(squeeze(write_ok(1,:)))*1.1]);
                end
                ylabel('write_ok Voltage [V]');
                xlabel('Time [Sec.]');


                subplot(5,1,4)
                plot(RAW_DATA.time_vect,squeeze(RW(1,:)))
                if min(squeeze(RW(1,:)))==max(squeeze(RW(1,:))) && min(squeeze(RW(1,:)))==0
                    axis([-inf inf -inf inf]);
                else
                    axis([-inf inf min(squeeze(RW(1,:)))*0.9 max(squeeze(RW(1,:)))*1.1]);
                end
                ylabel('RW Voltage [V]');
                xlabel('Time [Sec.]');


                subplot(5,1,5)
                plot(RAW_DATA.time_vect,squeeze(sys_clk(1,:)))
                if min(squeeze(sys_clk(1,:)))==max(squeeze(sys_clk(1,:))) && min(squeeze(sys_clk(1,:)))==0
                    axis([-inf inf -inf inf]);
                else
                    axis([-inf inf min(squeeze(sys_clk(1,:)))*0.9 max(squeeze(sys_clk(1,:)))*1.1]);
                end            
                ylabel('sys_clk Voltage [V]');
                xlabel('Time [Sec.]');

                print (gcf,'-dpng',fullfile(image_folder,sprintf('imagen_sys_crl_%d.png', pos_neg{pos_neg_i})))
                savefig(fullfile(image_folder,sprintf('imagen_sys_crl_%d.fig', pos_neg{pos_neg_i})));        
            end
        end
    end
end

