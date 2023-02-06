function processed_data = read_hData_MLP(finesim_results, partitions_N, neural_layers, results_folder, h_vPN, varargin)

    simEngine='FineSim';
    omit_size_search=0;
    pos_neg={'pos' 'neg'};
    
    if strcmpi(h_vPN,'vPN_DATA')
        hdr_name='vPN';
    else 
        hdr_name='h';
    end

    %% Optional values assignment
    if mod(nargin-5,2)
        display('A odd number of arguments have been introduced. Please revise the function call');
    else
        for i=1:nargin-5
            if strcmp(varargin{i},'networkSize')
                sizes=varargin{i+1};
                omit_size_search=1;
            end
        end
    end
 
    if ispc
        current_dir=split(finesim_results,"\");
        current_dir=current_dir(1:end-1);
        result_folder=cell2mat(join(current_dir,"\"));
        result_folder=strrep(root,'\','\\');
    elseif isunix
        current_dir=split(finesim_results,"/");
        current_dir=current_dir(1:end-1);
        result_folder=cell2mat(join(current_dir,"/"));
        %root_str=strrep(root,'\','\\');       
    end
    
    NW_files_id='MLP_';
    for layer_str_i=1:length(neural_layers)-1
        NW_files_id=strcat(NW_files_id,num2str(sizes(layer_str_i)),'by',num2str(sizes(layer_str_i+1)),'_');
    end
    NW_files_id=strcat(NW_files_id,'closed_loop');  
    
    sim_data=importdata(finesim_results);

    RAW_DATA.variable_mat=sim_data.data(:,1:end);
    RAW_DATA.variable_name_list=sim_data.colheaders(1:end);
    
    %% Get H(t) data for each Memdiodo
    if exist(strcat(NW_files_id,'_ideal_h.mat'),'file')
        load(strcat(NW_files_id,'_ideal_h.mat'));
        create_ops_file=0;
    else
        create_ops_file=1;
        search_ops=1;
    end
    data_cols=1:1:size(RAW_DATA.variable_name_list,2);
    for CPA_i=1:length(neural_layers)-1
        for pos_neg_i=1:2
            for partitions_i=1:partitions_N(CPA_i)
                for partitions_ii=1:partitions_N(CPA_i+1)
                    for i=1:sizes(CPA_i)/partitions_N(CPA_i)
                        for ii=1:sizes(CPA_i+1)/partitions_N(CPA_i+1)
                            found_h=0;

                            if create_ops_file==1
                                search_ops=1;
                            else
                                search_ops=0;
                            end

                            if create_ops_file==0
                                if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                                    if strcmpi(h_vPN,'vPN_DATA')
                                        string_search=sprintf('V(xmemd_network_P%d_%s.a%d-b%d_pn)', partitions_i, pos_neg{pos_neg_i},i,ii);                                        
                                    else 
                                        string_search=sprintf('V(xmemd_network_%s-%s_P%d-%d_%s.xmemdr%d-%dc%d-%d.h)', neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i},i,ii,ii,i);                                        
                                    end
                                else
                                    if strcmpi(h_vPN,'vPN_DATA')
                                        string_search=sprintf('V(xmemd_network_P%d_%s.a%d-b%d_pn)', partitions_i, pos_neg{pos_neg_i},i,ii);
                                    else 
                                        string_search=sprintf('V(xmemd_network_P%d_%s.a%d-b%d.h)', partitions_i, pos_neg{pos_neg_i}, i, ii);
                                    end
                                end

                                if ~eval(sprintf('strcmpi(RAW_DATA.variable_name_list{%s_%s%s_%s_loc(i+(partitions_i-1)*sizes(CPA_i)/partitions_N(CPA_i),ii+(partitions_ii-1)*sizes(CPA_i+1)/partitions_N(CPA_i+1),1)},string_search)',hdr_name, neural_layers{CPA_i}, neural_layers{CPA_i+1}, pos_neg{pos_neg_i}))
                                    search_ops=1;
                                else
                                    eval(sprintf('%s_%s%s_%s(i+(partitions_i-1)*sizes(CPA_i)/partitions_N(CPA_i),ii+(partitions_ii-1)*sizes(CPA_i+1)/partitions_N(CPA_i+1),1)=RAW_DATA.variable_mat(:,%s_%s%s_%s(i+(partitions_i-1)*sizes(CPA_i)/partitions_N,ii+(partitions_ii-1)*sizes(CPA_i+1)/partitions_N(CPA_i+1),2));', hdr_name, neural_layers{CPA_i}, neural_layers{CPA_i+1}, pos_neg{pos_neg_i}, hdr_name, neural_layers{CPA_i}, neural_layers{CPA_i+1}, pos_neg{pos_neg_i}));
                                end                
                            end

                            if search_ops==1
                                for iii=1:size(data_cols,2)              
                                    variable_name=RAW_DATA.variable_name_list{data_cols(iii)};
                                    %search Vp cols
                                    if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                                        if strcmpi(h_vPN,'vPN_DATA')
                                            string_search=sprintf('V(xmemd_network_%s-%s_P%d-%d_%s.r%d-c%d_pn)', neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i}, i, ii);
                                        else 
                                        string_search=sprintf('V(xmemd_network_%s-%s_P%d-%d_%s.xmemdr%d-%dc%d-%d.h)', neural_layers{CPA_i}, neural_layers{CPA_i+1}, partitions_i, partitions_ii, pos_neg{pos_neg_i},i,ii,ii,i);                                        
                                        end
                                    end
                                    if strcmpi(string_search,variable_name) 
                                        eval(sprintf('%s_%s%s_%s(i+(partitions_i-1)*sizes(CPA_i)/partitions_N(CPA_i),ii+(partitions_ii-1)*sizes(CPA_i+1)/partitions_N(CPA_i+1),1)=RAW_DATA.variable_mat(:,data_cols(iii));', hdr_name, neural_layers{CPA_i}, neural_layers{CPA_i+1}, pos_neg{pos_neg_i}));
                                        eval(sprintf('%s_%s%s_%s_loc(i+(partitions_i-1)*sizes(CPA_i)/partitions_N(CPA_i),ii+(partitions_ii-1)*sizes(CPA_i+1)/partitions_N(CPA_i+1),1)=data_cols(iii);',hdr_name, neural_layers{CPA_i}, neural_layers{CPA_i+1}, pos_neg{pos_neg_i}));
                                        found_h=1;
                                        if exist('delete_columns','var')==0
                                            delete_columns=iii;
                                        else
                                            delete_columns=vertcat(delete_columns,iii);
                                        end                   
                                    end                

                                    if found_h==1
                                        found_h=0;
                                        break
                                    end
                                end
%                                 if exist('delete_columns','var')
                                    delete_columns=sortrows(delete_columns,'descend');
%                                 else
%                                     prueba=1;
%                                 end
                                for iiii=1:size(delete_columns,1)
                                    data_cols(delete_columns(iiii))=[];
                                end
                                clear delete_columns
                            end
                        end
                    end
                end
            end
            eval(sprintf('processed_data.%s_%s{CPA_i,1}=%s_%s%s_%s;', hdr_name, pos_neg{pos_neg_i}, hdr_name, neural_layers{CPA_i}, neural_layers{CPA_i+1}, pos_neg{pos_neg_i}));
            eval(sprintf('processed_data.%s_%s_loc{CPA_i,1}=%s_%s%s_%s_loc;', hdr_name, pos_neg{pos_neg_i}, hdr_name, neural_layers{CPA_i}, neural_layers{CPA_i+1}, pos_neg{pos_neg_i}));
        end
    end
end

