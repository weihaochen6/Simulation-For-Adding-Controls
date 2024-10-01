function sect = addTableToSection(sect, tbl, titleText, decimalPlaces)
import mlreportgen.dom.*;
import mlreportgen.report.*;

% 创建表格标题
tblTitle = Paragraph(titleText);
tblTitle.Style = {Bold(false), FontSize('14pt')};
add(sect, tblTitle);

% 使用 Table(tbl) 自动转换 MATLAB table 为 Report Generator Table
rptTable = Table(tbl);

% 构造格式化字符串
numberFormatStr = sprintf('%%0.%df', decimalPlaces);


% 应用样式
rptTable.Style = {Width('100%'), Border('solid')};
rptTable.TableEntriesHAlign = 'center';
rptTable.TableEntriesVAlign = 'middle';
% 将格式化字符串用于设置表格的数字格式
rptTable.Style = [rptTable.Style {NumberFormat(numberFormatStr)}];
add(sect, rptTable);
end
