#https://www.geeksforgeeks.org/file-explorer-in-python-using-tkinter/
  
# import all components
# from the tkinter library
import tkinter as tk
  
# import filedialog module
from tkinter import filedialog
  
import pandas as pd
import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
import matplotlib
matplotlib.use("TkAgg")
import pickle 
import os
from functools import partial
import copy

#create settings
settings = {
    'chartWidth': 7,
    'chartHeight': 5,
    'sepChar': '_',
    'errType': "Standard Deviation",
    'initDir': "",
    'legendPos': 'Best',
    'legendPosAdjustX': '',
    'legendPosAdjustY': ''
    }

#create global variables for UI display
settingsfile = 'settings.pk'
errorTypes = ["Standard Deviation", "Standard Error"]
legendPositions = ['Best', 'Upper Right', 'Upper Left', 'Lower Left', 'Lower Right', 'Right', 'Center Left', 'Center Right', 'Lower Center', 'Upper Center', 'Center', 'Manual']       
currentFile = ''
data = {}
samples = [] #sample names (e.g. circuit1, circuit2)
categories = [] # sample groups (e.g. K562, K562 CD19)
headers = [] #FlowJo statistics; e.g. MFI, % positive

#create global variables for passing UI settings because tkinter is obnoxious
global entry_Width
global entry_Height
global entry_sepChar
global legendPosition
global entry_LegendAdjustX
global entry_LegendAdjustY
global list_Samples
global list_Groups
global list_Headers

#functions for moving listbox items

def moveup(listbox, name, *args):
   pos=1
   try:
       listbox.idxs = listbox.curselection()
       if not listbox.idxs:
           return
       for pos in listbox.idxs:
           if pos==0:
               continue
           text=listbox.get(pos)
           listbox.delete(pos)
           listbox.insert(pos-1, text)
           listbox.pop(pos)
           listbox.insert(pos-1, text)
           
   except:
       pass

   if name=='samples':
       global samples
       samples=list(listbox.get(0, tk.END)) #explicit cast to list because listbox.get() returns tuple because god knows why
       global list_Samples
       list_Samples.selection_set(pos-1)
   elif name=='groups':
       global categories
       categories=list(listbox.get(0, tk.END))
       global list_Groups
       list_Groups.selection_set(pos-1)
   elif name=='headers':
       global headers
       headers = list(listbox.get(0, tk.END))
       headers.insert(0, 'samples')
       global list_Headers
       list_Headers.selection_set(pos-1)

def movedown(listbox, name, *args):
   pos=1
   try:
       listbox.idxs = listbox.curselection()
       if not listbox.idxs:
           return
       for pos in listbox.idxs:
           text=listbox.get(pos)
           listbox.delete(pos)
           listbox.insert(pos+1, text)
           listbox.pop(pos)
           listbox.insert(pos+1, text)
           
   except:
       pass

   if name=='samples':
       global samples
       samples=list(listbox.get(0, tk.END))
       global list_Samples
       list_Samples.selection_set(pos+1)
   elif name=='groups':
       global categories
       categories=list(listbox.get(0, tk.END))
       global list_Groups
       list_Groups.selection_set(pos+1)
   elif name=='headers':
       global headers
       headers=list(listbox.get(0, tk.END))
       headers.insert(0, 'samples')
       global list_Headers
       list_Headers.selection_set(pos+1)


# Loading and plotting functions
def strToNum(s):
    try:
        return float(s)
    except:
        return 0
    
def getMaxColumnLength(col):
    max =  0
    for entry in col:
        l = len(str(entry))
        if (l > max):
            max = l
    return max

def computeError(numArr):
    eType = str(errorBarType.get())
    if (eType == "Standard Deviation"):
        return np.std(numArr, ddof=1)  #ddof=1 makes it compute sample, rather than population, stDev
    elif (eType == "Standard Error"):
        return sp.sem(numArr, ddof=1)

def parseData():
    
    df = pd.read_excel (currentFile)
    
    global data
    global samples
    global categories
    global headers
    
    #clear previous
    data = {}
    samples = []
    categories = []
    headers = []
    
    #prepare data structure
    headers = list(df.columns)
    categories, samples = getCategoriesAndSamples(df)
    
    for header in headers[1:]:
        data[header] = {}
        for sample in samples:
            data[header][sample] = {}
            for category in categories:
                data[header][sample][category] = {}
                data[header][sample][category]['values'] = []
    
    #populate data
    for index, row in df.iterrows():
        for i, value in enumerate(row):
            if i == 0:
                category, sample = getCategorySampleFromStr(value)
            else:
                header = headers[i]
                
                if isinstance(value, str) and ' %' in value:
                    value = value.replace(' %', '')
                    value = float(value)
                
                if math.isnan((value)):
                    continue

                else:
                    data[header][sample][category]['values'].append(value)
                
    #compute average and error
    for sample in samples:
        for header in headers[1:]:
            for category in categories:
                values = data[header][sample][category]['values']
                data[header][sample][category]['average'] = np.mean(values)
                data[header][sample][category]['error'] = computeError(values)

def divideAndPropagate(numerator, numErr, denominator, denomErr):
    quotient = numerator / denominator
    
    fac1 = (numErr / numerator)**2
    fac2 = (denomErr / denominator)**2
    err = quotient * (fac1 + fac2)**0.5
    return (quotient, err)
    

def normalizeBySample(inData, name):
    for header in inData:
        normalizations = copy.deepcopy(inData[header][name]) ## deep copy to avoid changing itself
        for sample in inData[header]:
            for category in inData[header][sample]:
                currentNorm = normalizations[category]
                (newMean, newErr)  = divideAndPropagate(inData[header][sample][category]['average'], inData[header][sample][category]['error'], currentNorm['average'], currentNorm['error'])
                inData[header][sample][category]['average'] = newMean
                inData[header][sample][category]['error'] = newErr     
    
    
def normalizeByGroup(inData, name):
    for header in inData:
        for sample in inData[header]:
            currentNorm = copy.deepcopy(inData[header][sample][name]) ## deep copy to avoid changing itself
            for category in inData[header][sample]:
                (newMean, newErr)  = divideAndPropagate(inData[header][sample][category]['average'], inData[header][sample][category]['error'], currentNorm['average'], currentNorm['error'])
                inData[header][sample][category]['average'] = newMean
                inData[header][sample][category]['error'] = newErr     
                
def normalizebyHeader(inData, name):
    normalizations = copy.deepcopy(inData[name])
    for header in inData:
        for sample in inData[header]:
            for category in inData[header][sample]:
                currentNorm = normalizations[sample][category] ## deep copy to avoid changing itself
                (newMean, newErr)  = divideAndPropagate(inData[header][sample][category]['average'], inData[header][sample][category]['error'], currentNorm['average'], currentNorm['error'])
                inData[header][sample][category]['average'] = newMean
                inData[header][sample][category]['error'] = newErr    

def normalizationRequired():
    normalization_Sample = normalizing_Sample.get()
    normalization_Group = normalizing_Group.get()
    normalization_Header = normalizing_Header.get()
    
    if normalization_Sample != "":
        return True
    if normalization_Group != "":
        return True
    if normalization_Header != "":
        return True
    return False

def normalizeData():
    
    output = copy.deepcopy(data)
    
    normalization_Sample = normalizing_Sample.get()
    normalization_Group = normalizing_Group.get()
    normalization_Header = normalizing_Header.get()
    
    if normalization_Sample != "":
        normalizeBySample(output, normalization_Sample)
    if normalization_Group != "":
        normalizeByGroup(output, normalization_Group)
    if normalization_Header != "":
        normalizebyHeader(output, normalization_Header)
        
    return output

def plotData():
    
    normalized = normalizeData() #will return the original input if no normalization required
    
    isNormalized = normalizationRequired()
    
    ind = np.arange(len(samples))  
    numGroups = len(categories)
    gap = .8 / len(categories)
    
    startPos = gap / numGroups * -1
    

    ### the official Roybal Lab colors...
    colors = [
        (127/255,127/255,127/255,255/255), #grey
        (0/255, 153/255, 153/255, 255/255), #teal
        (4/255, 7/255, 7/255, 255/255), #black
        (214/255, 22/255, 116/255, 255/255), #Axel Pink
        (21/255, 126/255, 170/255, 255/255), #Axel Dark Blue
        (253/255, 200/255, 121/255, 255/255) #Rayowis Pale Orange
        ]
    
    plotW = float(entry_Width.get())
    plotH = float(entry_Height.get())
    
    for header in headers[1:]:
        
        plt.figure(figsize = (plotW, plotH))
        
        for i, category in enumerate(categories):
            dataForCat = []
            err = []
            hasReplicates = False
            for sample in samples:
                dataForCat.append(normalized[header][sample][category]['average'])
                err.append(normalized[header][sample][category]['error'])
                if len(normalized[header][sample][category]['values']) > 1:
                    hasReplicates = True

            if hasReplicates == True:
                plt.bar(ind + startPos + (i+1) * gap, dataForCat, yerr=err, capsize=gap*10, width = gap, label = category, color = colors[i])
            else:
                plt.bar(ind + startPos + (i+1) * gap, dataForCat, width = gap, label = category, color = colors[i]) 
        
        plt.xticks([r + gap for r in range(len(samples))], samples, rotation = 45, fontsize=12, fontname='Arial', ha='right')
        if not isNormalized:
            plt.ylabel(header, fontsize=12, fontname='Arial')
        else:
            yLabel = header + " (Normalized)"
            plt.ylabel(yLabel, fontsize=12, fontname='Arial')
        
        
        
        legendLocation = legendPosition.get().lower()
        if not legendLocation == 'manual':
            plt.legend(loc=legendLocation, frameon=False, fancybox=False, shadow=False)
        else:
            adjX = strToNum(entry_LegendAdjustX.get())  / 10
            adjY = -strToNum(entry_LegendAdjustY.get()) / 10
            plt.legend(loc='lower center', bbox_to_anchor=(adjX, adjY), fancybox=False, shadow=False, frameon=False, ncol=len(categories))
            plt.subplots_adjust(bottom=0.25)
            
        plt.tight_layout()
        plt.show()
        
    defaultDir = os.path.dirname(currentFile)    
    matplotlib.rcParams['savefig.directory'] = defaultDir #set the default figure save directory to the data input directory

def getMaxValueCount(header, category):
    max = 0
    for sample in samples:
        if len(data[header][sample][category]['values']) > max:
            max = len(data[header][sample][category]['values'])
    return max

def GetLongestSampleNameLength():
    max = 0
    global samples
    for sample in samples:
        if len(sample) > max:
            max = len(sample)
    return max

def prepareData():
    isNormalized = normalizationRequired()
    normalizedData = {}
    if isNormalized:
        normalizedData = normalizeData()
        
    maxColumn = 0
    dfs = []
    for header in headers[1:]:
        #new excel sheet
        rows = []
        columnHeaders = []
        
        maxSubColumn = 0
        
        #write column headers
        for i, category in enumerate(categories):
            for j in range(getMaxValueCount(header, category)):
                columnHeaders.append(category)
                maxSubColumn += 1
        for i, category in enumerate(categories):
                columnHeaders.append('Average: ' + category)
                maxSubColumn += 1
        for i, category in enumerate(categories):
                columnHeaders.append('Error: ' + category)
                maxSubColumn += 1
        if isNormalized:
            for i, category in enumerate(categories):
                columnHeaders.append('Normalized Average: ' + category)
                maxSubColumn += 1
            for i, category in enumerate(categories):
                columnHeaders.append('Normalized Error: ' + category) 
                maxSubColumn += 1
        
        if maxSubColumn > maxColumn:
            maxColumn = maxSubColumn
        
        #write data
        for sample in samples:
            row = []
            #fill in raw data
            for i, category in enumerate(categories):
                row = row + data[header][sample][category]['values']
                #if there are less values for this sample than the max sample, pad this array with empty cells
                for j in range(getMaxValueCount(header, category) - len(data[header][sample][category]['values'])):
                    row = row + ['']
                
            #fill in means
            for i, category in enumerate(categories):
                row = row + [data[header][sample][category]['average']]
            #fill in errors
            for i, category in enumerate(categories):
                row = row + [data[header][sample][category]['error']]
                
            if isNormalized:
                #fill in normalized means
                for i, category in enumerate(categories):
                    row = row + [normalizedData[header][sample][category]['average']]
                #fill in normalized errors
                for i, category in enumerate(categories):
                    row = row + [normalizedData[header][sample][category]['error']]
            
            rows.append(row)
        df = pd.DataFrame(rows, index=samples, columns=columnHeaders)
        dfs.append(df)
    return dfs

def saveExcel():
    defaultDir = os.path.dirname(currentFile)
    types = [('Excel Files (*.xlsx)', '*.xlsx')]
    outputFile = filedialog.asksaveasfilename(initialdir = defaultDir, filetypes=types, defaultextension=types)
    if outputFile == '':
        return
    
    dfs = prepareData()
    
    writer = pd.ExcelWriter(outputFile, engine='xlsxwriter')
    wb = writer.book    
    cell_Format_Header = wb.add_format({'bold': True, 'align': 'left'})
    
    for index, header in enumerate(headers[1:]):
        df = dfs[index]
        df.to_excel(writer, sheet_name=header, startrow=1, header=False) #turn off auto-header to alow custom format
        
        #write headers with custom format https://xlsxwriter.readthedocs.io/example_pandas_header_format.html
        ws = writer.sheets[header]
        writer.sheets[header].set_column(0,0,GetLongestSampleNameLength() * 1.1) #set first column width according to max sample name length
        for col_num, value in enumerate(df.columns.values):  
            ws.write(0, col_num + 1, value, cell_Format_Header)
            #set column width according to header width
            writer.sheets[header].set_column(col_num+1, col_num+1, len(value) * 1.1)
        
        #Make row headers also left-justified
        s_idx = 0
        for row_name, row in df.iterrows(): 
            ws = writer.sheets[header]
            ws.write(s_idx + 1, 0, row_name, cell_Format_Header)
            s_idx += 1

    writer.save()

def saveCSV():
    defaultDir = os.path.dirname(currentFile)
    outputDir = filedialog.askdirectory(initialdir = defaultDir)
    if outputDir == '':
        return
    
    dfs = prepareData()
    for index, header in enumerate(headers[1:]):
        df = dfs[index]
        outputPath = os.path.join(outputDir, header) + '.csv'
        
        df.to_csv(outputPath)
     
def unique(input):
    
    output = []
    for elem in input:
        if not elem in output:
            output.append(elem)
    return output

def getCategorySampleFromStr(str):
    separator = entry_sepChar.get()
    splt = str.split(separator)
    cat = splt[0]
    sample = separator.join(splt[1:])
    return cat, sample

def getCategoriesAndSamples(df):
    
    categories = []
    samples = []
    
    firstCol = df.iloc[:,0].values
    for name in firstCol:
        cat, sample = getCategorySampleFromStr(name)
        categories.append(cat)
        samples.append(sample)
    
    categories = unique(categories)
    samples = unique(samples)
    
    categories.sort(key=len)
    
    return unique(categories), unique(samples)

    
def chooseFile(initDir):
    chosenFile = filedialog.askopenfilename(initialdir = initDir,
                                          title = "Select a File",
                                          filetypes = (("Excel files",
                                                        "*.xls*"),
                                                       ("all files",
                                                        "*.*")))
    
    if chosenFile == '':
        return
    if not (os.path.exists(chosenFile)):
        print(chosenFile + " is not a valid file path")
        return
    
    
    
    global currentFile 
    currentFile = chosenFile
    
def loadData(initDir, list_Samples, list_Groups, list_Headers):
    chooseFile(initDir)
    if currentFile == '':
        return
    parseData()
    
    #update UI lists
    list_Samples.delete(0, tk.END)
    global samples
    for sample in samples:
        list_Samples.insert(tk.END, sample)
        
    list_Groups.delete(0, tk.END)
    global categories
    for cat in categories:
        list_Groups.insert(tk.END, cat)
    
    list_Headers.delete(0, tk.END)
    global headers
    for h in headers[1:]:
        list_Headers.insert(tk.END, h)
        
def selectNormalizationSample(selectedItem):
    
    global normalizing_Sample
    normalizing_Sample.set(selectedItem.get(tk.ACTIVE))
        
def selectNormalizationGroup(selectedItem):
    
    global normalizing_Group
    normalizing_Group.set(selectedItem.get(tk.ACTIVE))
        
def selectNormalizationHeader(selectedItem):
    
    global normalizing_Header
    normalizing_Header.set(selectedItem.get(tk.ACTIVE))
    
def clearNormalizingSample():
    global normalizing_Sample
    normalizing_Sample.set('')

def clearNormalizingGroup():
    global normalizing_Group
    normalizing_Group.set('')
    
def clearNormalizingHeader():
    global normalizing_Header
    normalizing_Header.set('')

#Save settings on quit
def on_closing():  
    
    global entry_Width
    global entry_Height
    global entry_sepChar
    global legendPosition
    
    try:
        settings['chartWidth'] = int(entry_Width.get())
        settings['chartHeight'] = int(entry_Height.get())
        settings['sepChar'] = entry_sepChar.get()
        settings['errType'] = str(errorBarType.get())
        settings['initDir'] = os.path.dirname(currentFile)
        settings['legendPos'] = legendPosition.get()
        settings['legendPosAdjustX'] = entry_LegendAdjustX.get()
        settings['legendPosAdjustY'] = entry_LegendAdjustY.get()

        with open(settingsfile, "wb") as f:
            pickle.dump(settings, f, pickle.HIGHEST_PROTOCOL)

        window.destroy()
    except ValueError as e:
        print("Invalid value:", e)

#Load settings on open
try:
    with open(settingsfile, 'rb') as f:
        try:
            settings = pickle.load(f)
        except: 
            pass
except:
    pass
                                                                                          
# Create the root window
window = tk.Tk()

canvas1 = tk.Canvas(window, width = 450, height = 500,  relief = 'raised', bg='white')
canvas1.pack()
  
# Set window title
window.title('File Explorer')
  
window.wm_title('FlowJo Table Plotting For Lazy Academics')  

# Create Lists to Display
label_LS = tk.Label(window, text = "Loaded Samples ", bg="white")
list_Samples = tk.Listbox(window, listvariable=samples, selectmode=tk.SINGLE, width=20, height=10)

moveSampleUp = partial(moveup, list_Samples, 'samples')
button_MoveSampleUp = tk.Button(window,
                        text = u"\u1403", fg='white', bg='black',
                        command = moveSampleUp)
moveSampleDown = partial(movedown, list_Samples, 'samples')
button_MoveSampleDown = tk.Button(window,
                        text = u"\u1401", fg='white', bg='black',
                        command = moveSampleDown)

label_LG = tk.Label(window, text = "Loaded Groups ", bg="white")
list_Groups = tk.Listbox(window, listvariable=categories, selectmode=tk.SINGLE, width=20, height=10)

moveGroupUp = partial(moveup, list_Groups, 'groups')
button_MoveGroupUp = tk.Button(window,
                        text = u"\u1403", fg='white', bg='black',
                        command = moveGroupUp)
moveGroupDown = partial(movedown, list_Groups, 'groups')
button_MoveGroupDown = tk.Button(window,
                        text = u"\u1401", fg='white', bg='black',
                        command = moveGroupDown)

label_LH = tk.Label(window, text = "Loaded Headers ", bg="white")
list_Headers = tk.Listbox(window, listvariable=headers, selectmode=tk.SINGLE, width=20, height=10)

moveHeaderUp = partial(moveup, list_Headers, 'headers')
button_MoveHeaderUp = tk.Button(window,
                        text = u"\u1403", fg='white', bg='black',
                        command = moveHeaderUp)
moveHeaderDown = partial(movedown, list_Headers, 'headers')
button_MoveHeaderDown = tk.Button(window,
                        text = u"\u1401", fg='white', bg='black',
                        command = moveHeaderDown)

label_SS = tk.Label(window, text = "Normalize to Sample: ", bg="white")
label_SG = tk.Label(window, text = "Normalize to Group: ", bg="white")
label_SH = tk.Label(window, text = "Normalize to Header: ", bg="white")

#Define data loading functions
loadFromLastDir = partial(loadData, settings['initDir'], list_Samples, list_Groups, list_Headers) 
button_getFile = tk.Button(window,
                        text = "Load File", fg='white', bg='black',
                        command = loadFromLastDir)

button_plotFile = tk.Button(window,
                        text = "Plot Data", fg='white', bg='black',
                        command = plotData)

button_writeExcel = tk.Button(window,
                        text = "Save to Excel", fg='white', bg='black',
                        command = saveExcel)

button_writeCSV = tk.Button(window,
                        text = "Save to CSV", fg='white', bg='black',
                        command = saveCSV)
  
button_exit = tk.Button(window,
                     text = "Exit", fg='white', bg='black',
                     command = on_closing)

#Define other UI elements
label_Width = tk.Label(window, text = "Chart Width: ", bg="white")
entry_Width = tk.Entry(window, width=5)
entry_Width.insert(-1, settings['chartWidth'])

label_Height = tk.Label(window, text = "Chart Height: ", bg="white")
entry_Height = tk.Entry(window, width=5)
entry_Height.insert(-1, settings['chartHeight'])

label_SepChar = tk.Label(window, text = "[Group][Sample] Separator: ", bg="white")
entry_sepChar = tk.Entry(window, width=5)
entry_sepChar.insert(-1, settings['sepChar'])

label_Err = tk.Label(window, text = "Error Bars: ", bg="white")
errorBarType = tk.StringVar(window)
errorBarType.set(settings['errType'])
drop_Err = tk.OptionMenu(window, errorBarType, *errorTypes)
drop_Err.config(bg = "black", fg="white") 
drop_Err.pack()

label_LegendPos = tk.Label(window, text = "Legend Position: ", bg="white")
legendPosition = tk.StringVar(window)
legendPosition.set(settings['legendPos'])
drop_LegendPos = tk.OptionMenu(window, legendPosition, *legendPositions)
drop_LegendPos.config(bg = "black", fg="white") 
drop_LegendPos.pack()

label_LegendAdjustX = tk.Label(window, text = "X: ", bg="white")
entry_LegendAdjustX = tk.Entry(window, width=5)
entry_LegendAdjustX.insert(-1, settings['legendPosAdjustX'])
label_LegendAdjustY = tk.Label(window, text = "Y: ", bg="white")
entry_LegendAdjustY = tk.Entry(window, width=5)
entry_LegendAdjustY.insert(-1, settings['legendPosAdjustY'])

normalizing_Sample = tk.StringVar(window)
normalizing_Group = tk.StringVar(window)
normalizing_Header = tk.StringVar(window)

label_nS = tk.Label(window, textvariable = normalizing_Sample, bg="white")
label_nG = tk.Label(window, textvariable = normalizing_Group, bg="white")
label_nH = tk.Label(window, textvariable = normalizing_Header, bg="white")

setNormToSelSample = partial(selectNormalizationSample, list_Samples)
button_normToSelSample = tk.Button(window,
                        text = "Normalize to Selected", fg='white', bg='black',
                        command = setNormToSelSample)

button_ClearNormToSelSample = tk.Button(window,
                        text = "Clear", fg='white', bg='black',
                        command = clearNormalizingSample)

setNormToSelGroup = partial(selectNormalizationGroup, list_Groups)
button_normToSelGroup = tk.Button(window,
                        text = "Normalize to Selected", fg='white', bg='black',
                        command = setNormToSelGroup)

button_ClearNormToSelGroup = tk.Button(window,
                        text = "Clear", fg='white', bg='black',
                        command = clearNormalizingGroup)

setNormToSelHeader = partial(selectNormalizationHeader, list_Headers)
button_normToSelHeader = tk.Button(window,
                        text = "Normalize to Selected", fg='white', bg='black',
                        command = setNormToSelHeader)

button_ClearNormToSelHeader = tk.Button(window,
                        text = "Clear", fg='white', bg='black',
                        command = clearNormalizingHeader)

#Draw UI elements  
canvas1.create_window(50, 25, window=label_Width)
canvas1.create_window(110, 25, window=entry_Width)
canvas1.create_window(250, 25, window=label_Height)
canvas1.create_window(310, 25, window=entry_Height)

canvas1.create_window(85, 50, window=label_SepChar)
canvas1.create_window(175, 50, window=entry_sepChar)
canvas1.create_window(240, 50, window=label_Err)
canvas1.create_window(350, 50, window=drop_Err)

canvas1.create_window(60, 80, window=label_LegendPos)
canvas1.create_window(175, 80, window=drop_LegendPos)
canvas1.create_window(275, 80, window=label_LegendAdjustX)
canvas1.create_window(305, 80, window=entry_LegendAdjustX)
canvas1.create_window(375, 80, window=label_LegendAdjustY)
canvas1.create_window(405, 80, window=entry_LegendAdjustY)

canvas1.create_window(50, 115, window=button_getFile)
canvas1.create_window(125, 115, window=button_plotFile)
canvas1.create_window(210, 115, window=button_writeExcel)
canvas1.create_window(305, 115, window=button_writeCSV)
canvas1.create_window(380, 115, window=button_exit)

canvas1.create_window(75, 150, window=label_LS)
canvas1.create_window(75, 250, window=list_Samples)
canvas1.create_window(128, 150, window=button_MoveSampleUp)
canvas1.create_window(20, 150, window=button_MoveSampleDown)

canvas1.create_window(225, 150, window=label_LG)
canvas1.create_window(225, 250, window=list_Groups)
canvas1.create_window(275, 150, window=button_MoveGroupUp)
canvas1.create_window(173, 150, window=button_MoveGroupDown)

canvas1.create_window(375, 150, window=label_LH)
canvas1.create_window(375, 250, window=list_Headers)
canvas1.create_window(428, 150, window=button_MoveHeaderUp)
canvas1.create_window(320, 150, window=button_MoveHeaderDown)

canvas1.create_window(75, 400, window=button_normToSelSample)
canvas1.create_window(225, 400, window=button_normToSelGroup)
canvas1.create_window(375, 400, window=button_normToSelHeader)

canvas1.create_window(75, 425, window=button_ClearNormToSelSample)
canvas1.create_window(225, 425, window=button_ClearNormToSelGroup)
canvas1.create_window(375, 425, window=button_ClearNormToSelHeader)

canvas1.create_window(75, 450, window=label_SS)
canvas1.create_window(75, 475, window=label_nS)

canvas1.create_window(225, 450, window=label_SG)
canvas1.create_window(225, 475, window=label_nG)

canvas1.create_window(375, 450, window=label_SH)
canvas1.create_window(375, 475, window=label_nH)
     
window.protocol("WM_DELETE_WINDOW", on_closing)

# Let the window wait for any events
window.mainloop()