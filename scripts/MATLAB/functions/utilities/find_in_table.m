function [SUB_TABLE] = find_in_table(SUB_TABLE,field,variable)
    if iscell(variable)
        variable=variable{1,1};
    end
    if ~isnumeric(variable) & ~ismissing(variable) & ~isscalar(variable) & ~iscell(variable)
        SUB_TABLE=eval(sprintf('SUB_TABLE(strcmpi(SUB_TABLE.%s,variable),:)',field));         
    elseif isnumeric(variable) & ~ismissing(variable) & isscalar(variable) & ~iscell(variable)
        SUB_TABLE=eval(sprintf('SUB_TABLE(SUB_TABLE.%s==variable,:)',field));   
    elseif ismissing(variable) & ~iscell(variable)
        %TF = ismissing(SUB_TABLE,{'' '.' 'NA' NaN -99});
        eval(sprintf('TF = ismissing(SUB_TABLE.%s);',field));
        SUB_TABLE = SUB_TABLE(any(TF,2),:);
    elseif isnumeric(variable) & ~isscalar(variable) & ~iscell(variable)
        eval(sprintf('SUB_TABLE=SUB_TABLE(find(cellfun(@(c) isequaln(c, variable), SUB_TABLE.%s)),:);',field));
    elseif iscell(variable)
        eval(sprintf('variable_col=SUB_TABLE.%s;',field));
        for n=1:size(variable_col,1)
            aux_var=setdiff(variable_col{n,1},variable);
            if isempty(aux_var)
                if ~exist('aux_selected_rows','var')
                    aux_selected_rows=n;
                else
                    aux_selected_rows=[aux_selected_rows,n];
                end
            end
        end
        SUB_TABLE=SUB_TABLE(aux_selected_rows,:);
        clear aux_selected_rows
    end
end

