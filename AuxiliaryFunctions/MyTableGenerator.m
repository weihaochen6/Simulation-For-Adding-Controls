classdef MyTableGenerator<handle
    properties(Access = public)
        Data = []
        rowNames = []
        columnNames = []
        DecimalNum = 4
    end

    %% Construct Methods
    methods(Access = public)
        function obj = MyTableGenerator(inputStruct)
            % Allow create empty instance
            if nargin == 0
                return;
            end

            % Validate that the input argument is a structure
            if ~isstruct(inputStruct)
                error([class(obj) ':Constructor:InputNotStruct'], 'Input must be a structure.');
            end

            % Proceed if input is a structure
            propNames = fieldnames(inputStruct);
            for i = 1:length(propNames)
                propName = propNames{i};
                if isprop(obj, propName)
                    obj.(propName) = inputStruct.(propName);
                end
            end
        end
    end

    %% Main Public Methods
    methods(Access = public)
        function getReportTable(obj)
            import mlreportgen.dom.*
            import mlreportgen.report.*

            if(isempty(obj.Data))
                error('The Data is not imported.');
            end
            [DataRowNum,DataColumn] = size(obj.Data);

            if isempty(obj.rowNames)
                % 生成行名
                obj.rowNames = cell(DataRowNum, 1);  % 初始化行名的 cell 数组
                for i = 1:DataRowNum
                    obj.rowNames{i} = sprintf('Row%d', i);  % 格式化生成行名，如 Row1, Row2, ...
                end
            end

            if isempty(obj.columnNames)
                % 生成列名
                obj.columnNames = cell(1, DataColumn);  % 初始化列名的 cell 数组
                for j = 1:DataColumn
                    obj.columnNames{j} = sprintf('Col%d', j);  % 格式化生成列名，如 Col1, Col2, ...
                end
            end

            % 准备包含行名和列名的数据
            fullData = [cellstr(obj.columnNames); num2cell(data)];
            fullData = [cellstr([{' '}; obj.rowNames']) fullData];

            % 创建表格
            table = Table(fullData);
            table.Border = 'solid';
            table.TableEntriesStyle = {FontFamily('Times Newman'), Width('1in'), Color('black')};
            table.Style = {Width('100%'), Border('solid')};
            table.TableEntriesHAlign = 'center';
            table.TableEntriesVAlign = 'middle';
            % 构造格式化字符串
            numberFormatStr = sprintf('%%0.%df', obj.DecimalNum);

            % 将格式化字符串用于设置表格的数字格式
            table.Style = [table.Style {NumberFormat(numberFormatStr)}];
        end
    end
end

