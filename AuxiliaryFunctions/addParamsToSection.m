function section = addParamsToSection(section, params)
import mlreportgen.report.*
import mlreportgen.dom.*
% 创建表格并填入参数数据
paramNames = fieldnames(params);
numParams = length(paramNames);
paramTableData = cell(numParams, 2);
for i = 1:numParams
    paramName = paramNames{i};
    paramValue = params.(paramName);
    paramTableData{i, 1} = paramName;
    if isnumeric(paramValue) || islogical(paramValue)
        paramTableData{i, 2} = num2str(paramValue);
    else
        paramTableData{i, 2} = char(paramValue);
    end
end
% 创建MATLAB表格
T_params = table(paramTableData(:,1), paramTableData(:,2), ...
    'VariableNames', {'Parameter', 'Value'});
% 将表格添加到报告章节
tbl = Table(T_params);
tbl.Style = {Width('100%'), Border('solid'), ...
    ColSep('solid'), RowSep('solid')};
tbl.TableEntriesVAlign = 'middle';
tbl.TableEntriesHAlign = 'center';
tbl.ColSep = 'single';
tbl.RowSep = 'single';
tbl.Border = 'single';

% 设置标题和添加表格到章节
tblTitle = Text('Parameters Used');
tblTitle.Bold = false;
tblTitle.FontSize = '12';
add(section, tblTitle);
add(section, tbl);
return
end