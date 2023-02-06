function processed_data = read_opData(finesim_results,results_folder,varargin)

    simEngine='FineSim';
    omit_size_search=0;
    %% Optional values assignment
    if mod(nargin-2,2)
        display('A odd number of arguments have been introduced. Please revise the function call');
    else
        for i=1:nargin-2
            if strcmp(varargin{i},'networkSize')
                size_A=varargin{i+1}(1);
                size_B=varargin{i+1}(2);
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
    
    size_matrix=[size_A size_B];
    num_neurons_out=size_matrix(2);  
    
    sim_data=importdata(finesim_results);

    RAW_DATA.variable_mat=sim_data.data(:,1:end);
    RAW_DATA.variable_name_list=sim_data.colheaders(1:end);
    
    %% Get H(t) data for each Memdiodo
    if exist(strcat('NW_',num2str(size_matrix(1)),'by',num2str(size_matrix(2)),'_closed_loop_ideal_op.mat'),'file')
        load(strcat('NW_',num2str(size_matrix(1)),'by',num2str(size_matrix(2)),'_closed_loop_ideal_op.mat'));
        create_ops_file=0;
    else
        create_ops_file=1;
        search_ops=1;
    end
    data_cols=1:1:size(RAW_DATA.variable_name_list,2);
    for pos_neg=1:2
        for i=1:size_A
            for ii=1:size_B
                found_I=0;
                found_Vpn=0;
                
                if create_ops_file==1
                    search_ops=1;
                else
                    search_ops=0;
                end
                
                if create_ops_file==0                
                    if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                        string_search=sprintf('i(xmmed_network_%d.rcsa%d-%db%d-%d_1)', pos_neg, i, ii, ii, i);
                    else
                        string_search=sprintf('i(mmed_network_%d:rcsa%d-%db%d-%d_1)', pos_neg, i, ii, ii, i);
                    end

                    if ~eval(sprintf('strcmpi(RAW_DATA.variable_name_list{I%d(i,ii,2)},string_search)',pos_neg))
                        search_ops=1;
                    else
                        eval(sprintf('I%d(i,ii,1)=RAW_DATA.variable_mat(:,I%d(i,ii,2));', pos_neg, pos_neg));
                    end

                    if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                        string_search=sprintf('V(xmmed_network_%d.a%d-b%d_pn)', pos_neg, i, ii);
                    else
                        string_search=sprintf('v(xmmed_network_%d.a%d-b%d_pn)', pos_neg, i, ii);
                    end

                    if ~eval(sprintf('strcmpi(RAW_DATA.variable_name_list{Vpn%d(i,ii,2)},string_search)',pos_neg))
                        search_ops=1;
                    else
                        eval(sprintf('Vpn%d(i,ii,1)=RAW_DATA.variable_mat(:,Vpn%d(i,ii,2));', pos_neg, pos_neg));
                    end                
                end
                
                if search_ops==1
                    for iii=1:size(data_cols,2)
    %                     %Search the H cols
                         variable_name=RAW_DATA.variable_name_list{data_cols(iii)};

                        %search I cols
                        if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                            string_search=sprintf('i(xmmed_network_%d.rcsa%d-%db%d-%d_1)', pos_neg, i, ii, ii, i);
                        else
                            string_search=sprintf('i(mmed_network_%d:rcsa%d-%db%d-%d_1)', pos_neg, i, ii, ii, i);
                        end
                        if strcmpi(string_search,variable_name) 
                            eval(sprintf('I%d(i,ii,1)=RAW_DATA.variable_mat(:,data_cols(iii));', pos_neg));
                            eval(sprintf('I%d(i,ii,2)=data_cols(iii);', pos_neg));                        
                            found_I=1;
                            %delete_I=iii;
                            if exist('delete_columns','var')==0
                                delete_columns=iii;
                            else
                                delete_columns=vertcat(delete_columns,iii);
                            end
                        end                

    %                     %search Vp cols
                        if strcmp(simEngine,'HSPICE') || strcmp(simEngine,'FineSim')
                            string_search=sprintf('V(xmmed_network_%d.a%d-b%d_pn)', pos_neg, i, ii);
                        else
                            string_search=sprintf('v(xmmed_network_%d.a%d-b%d_pn)', pos_neg, i, ii);
                        end
                        if strcmpi(string_search,variable_name) 
                            eval(sprintf('Vpn%d(i,ii,1)=RAW_DATA.variable_mat(:,data_cols(iii));', pos_neg));
                            eval(sprintf('Vpn%d(i,ii,2)=data_cols(iii);', pos_neg));
                            found_Vpn=1;
                            %delete_Vp=iii;
                            if exist('delete_columns','var')==0
                                delete_columns=iii;
                            else
                                delete_columns=vertcat(delete_columns,iii);
                            end                   
                        end                

                        if found_I==1 && found_Vpn==1
                            found_I=0;
                            found_Vpn=0;
                            break
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
        eval(sprintf('processed_data.I%d=I%d;', pos_neg, pos_neg));
        eval(sprintf('processed_data.Vpn%d=Vpn%d;', pos_neg, pos_neg));
    end
end

