function [RESULTS_TABLE] = save2table(cfg_variables,results_directory)
    
    save_file=0;
    table_entry = struct2table(cfg_variables,'AsArray',true);
    
    if exist(fullfile(results_directory,'RESULTS_TABLE.mat'))
        load(fullfile(results_directory,'RESULTS_TABLE.mat'));
    end
    
    if ~exist('RESULTS_TABLE','var')
        RESULTS_TABLE=table_entry;
        save_file=1;
    else
        SUB_TABLE=RESULTS_TABLE;
        for m=1:length(RESULTS_TABLE.Properties.VariableNames)
            variable=table2array(table_entry(1,m));
            SUB_TABLE=find_in_table(SUB_TABLE,table_entry.Properties.VariableNames{1,m},variable);
            if isempty(SUB_TABLE)
                %if length(table_entry.Properties.VariableNames)==length(RESULTS_TABLE.Properties.VariableNames) && isempty(setdiff(fields(cfg_variables),fields(table2struct(RESULTS_TABLE(1,:))))) && isempty(setdiff(fields(table2struct(RESULTS_TABLE(1,:))),fields(cfg_variables))) && length(cfg_variables.series_resistance)==1
                    RESULTS_TABLE=[RESULTS_TABLE;table_entry];
                %else
                %    display(cfg_variables.results_folder)
                %end
                save_file=1;
                break;
            end
        end
    end
    if save_file==1
        save_file=0;
        save(fullfile(results_directory,'RESULTS_TABLE.mat'),'RESULTS_TABLE');
    end
end

