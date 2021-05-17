# FlowJo Plotter

## Overview
This script plots FlowJo-Exported Data to a chart or conveniently formatted table.

Expected input: Table exported from FlowJo, with sample names renamed as Group_Sample in the first column.

## Installation

### Windows

1: Download and install Anaconda from [here](https://docs.anaconda.com/anaconda/install/)

2: Download the .py and .bat file from this GitHub Repo, or clone the whole repository.

3: Edit the .bat file in a text editor to replace my installation paths with yours for Anaconda\Scripts, Anaconda, and FlowJo Plotter. Save the bat file.

4: Double click on the bat file to launch.

### Mac

1: Go to Best Buy

2: Buy a Windows PC

3: See **Windows** section

OR

Find someone with a Mac who can help you with installing Anaconda and launching Python scripts

## Options

**Chart Width**: Width of output chart(s) in inches

**Chart Height**: Height of output chart(s) in inches

**[Group][Sample] Separator**: The character or string by which the Group and Sample in the first column are separated. Default is _

**Error Bars**: Chart and Table error can be set to sample standard deviation or S.E.M.

**Legend Position**: Controls where the legend will be (standard MatPlotLib positions), or Manual. X and Y adjustments will only be used when position is Manual

**Load File**: Prompts user for the input table exported from FlowJo

**Plot Data**: Creats charts from data loaded in file. One chart will be created per column in the input table.

**Save to Excel**: Creates an Excel file with the input data pre-grouped for convenient plotting, or for you to verify the math by-hand.

**Save to CSV**: Creates an CSV file with the input data pre-grouped for convenient plotting. Pretty much same as **Save to Excel** but without the default formatting.

**Exit**: Self-explanatory, hopefully.

**Normalize to Selected**: Data can be normalized to a sample, group, or column. Up to one of each can be chosen. Error will be propagated.
**Clear** Clears the selected normalizing sample/group/column

## Example Usages

Background Excel file is the FlowJo-outputted table. Only modifications are renaming the samples & column headers.
Bottom Excel file is the new one generated by FlowJo Plotter's **Save to Excel** function.

1: Reporting receptor activation - raw MFI values
![Screenshot](https://raw.github.com/dpiraner/FlowJoPlotte/main/Docs/Img/Example1.PNG)

2: Reporting receptor activation - fold change relative to vanilla cells
![Screenshot](https://raw.github.com/dpiraner/FlowJoPlotte/main/Docs/Img/Example2.PNG)

3: Reporting receptor activation - fold change relative to Reporter-Only T-cells
![Screenshot](https://raw.github.com/dpiraner/FlowJoPlotte/main/Docs/Img/Example3.PNG)

4: Killing assay - raw K562 counts
![Screenshot](https://raw.github.com/dpiraner/FlowJoPlotte/main/Docs/Img/Example4.PNG)

5: Killing assay - K562 counts normalized to vanilla counts for each type of T-cell
![Screenshot](https://raw.github.com/dpiraner/FlowJoPlotte/main/Docs/Img/Example5.PNG)

6: Killing assay - K562 counts normalized to counting beads, and to control wells with untransduced T-cells.
![Screenshot](https://raw.github.com/dpiraner/FlowJoPlotte/main/Docs/Img/Example6.PNG)
