#FlowJo Plotter

##Overview
This script plots FlowJo-Exported Data to a chart or Excel table.
Expected input: Table exported from FlowJo, with sample names renamed as Group_Sample in the first column.

##Options

**Chart Width**: Width of output chart(s) in inches
**Chart Height**: Height of output chart(s) in inches
**[Group][Sample] Separator**: The character or string by which the Group and Sample in the first column are separated. Default is _
**Error Bars**: Chart and Table error can be set to sample standard deviation or S.E.M.
**Load File**: Prompts user for the input table
**Plot Data**: Creats charts from data loaded in file. One chart will be created per column in the input table.
**Save to Excel**: Creates an Excel file with the input data pre-grouped for convenient plotting
**Exit**: Self-explanatory

**Normalize to Selected**: Data can be normalized to a sample, group, or column. Up to one of each can be chosen. Error will be propagated.
**Clear** Clears the selected normalizing sample/group/column
