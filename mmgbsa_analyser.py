#!/usr/bin/env python
__author__ = 'Khaled Tumbi and Naeem Attari'
import numpy as np
import os
import argparse
def xtractor(fileName):

    '''
    THis function will extracts relevent energy values from FINAL_DECOMP_MMPBSA.dat file.
    Input details Details:
    fieName:
    lres, rres: Number of residues in ligand and receptor
    lGBT,rGBT,lPBT,rPBT: np zero arrays, Because I dont know how to use them. hahahahahaha.

    OutPut Details:
    :param fileName: list of filename[s]. All values from these files will be averaged. and then std. deviation will be\
     calculated.
    :param lres: Number of ligand residues. This is useful for Protein-Protien interaction calculations.
    :param rres: Number of receptor residues.
    :param lGBT: np zero arrays, Because I dont know how to use them. hahahahahaha.
    :param rGBT: .............do.......................................
    :param lPBT: ....................do....................
    :param rPBT: ...........do.....................
    :return: retune many arrays.
    '''
    import os
    import numpy as np
    from StringIO import StringIO
    lig = []
    pro = []
    for i in range (len(fileName)):
        fname = fileName[i]
        if os.path.isfile(fname) == True:
	    print (fname)
            fin = open (fname, 'r')
            fileData = fin.read()
            HeaderLine = ',,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean,Avg.,Std. Dev.,Std. Err. of Mean\r'
            fileDataArray = fileData.split('\n')
            indices = [i for i, x in enumerate(fileDataArray) if x == HeaderLine]
            print (indices)
            for i in range (1, len(indices)+1):
                SkipFooter = fileData.count('\n')+1-(indices[(i-1)]+lres+rres+9)
                #print ('This is SkipFooter', SkipFooter)
                #print ('This is SkipHeader', indices[i-1]+1)
                Data = np.genfromtxt(StringIO(fileData), dtype=float, delimiter=",",skip_header=indices[i-1]+1,skip_footer=SkipFooter+i+3, unpack=True, usecols = (17))
                #print (Data)
		d.append(Data)
		rGBT = np.delete(Data,np.s_[rres:])
		#print (rGBT)
                if i-1 == 0:
		    rGBT = np.delete(Data,np.s_[rres:])
		    pro.append(rGBT)
		    lGBT = np.delete(Data,np.s_[:rres])
		    lig.append(lGBT)
                    #print ('Entered Here')
                    #rGBT = np.vstack((rGBT,np.delete(Data,np.s_[rres:])))
                    #lGBT = np.vstack((lGBT,np.delete(Data,np.s_[:rres])))
                    residueNames = np.genfromtxt(StringIO(fileData), dtype=str, delimiter=",",skip_header=indices[i-1]+1,skip_footer=SkipFooter+i+3, unpack=True, usecols = (0))
                    lresName = np.delete(residueNames,np.s_[:rres])

                    rresName = np.delete(residueNames,np.s_[rres:])
		    #print (lresName,rresName )
                #if i-1 == 2:
                #    rPBT = np.vstack((rPBT,np.delete(Data,np.s_[rres:])))
                #   lPBT = np.vstack((lPBT,np.delete(Data,np.s_[:rres])))
        fin.close()
	
    #rGBTm = np.mean((np.delete((np.transpose(rGBT)),0,axis=1)),axis=1)
    #lGBTm = np.mean((np.delete((np.transpose(lGBT)),0,axis=1)),axis=1)
    	

    #rGBTstd = np.std((np.delete((np.transpose(rGBT)),0,axis=1)),axis=1)
    #lGBTstd = np.std((np.delete((np.transpose(lGBT)),0,axis=1)),axis=1)
    print (np.transpose(lig))
    lGBTm = (np.mean(np.transpose(lig),axis=1))
    rGBTm = (np.mean(np.transpose(pro),axis=1))
    lGBTstd = (np.std(np.transpose(lig),axis=1))
    rGBTstd = (np.std(np.transpose(pro),axis=1))

    return lGBTm,rGBTm,lresName, rresName,  lGBTstd, rGBTstd

def sorterlig(resName, meanGB, stdGB ):
    '''
    It sort stack and sort the provided array.
    :param resName:
    :param meanGB:
    :param stdGB:
    :param meanPB:
    :param stdPB:
    :return:
    '''
    import numpy as np
    data = np.vstack((resName,(meanGB.astype(np.float)),stdGB.astype(np.float)))
    #SortedData = data[:,np.argsort(data[1].astype(float))]
    return data

def sorter(resName, meanGB, stdGB ):
    '''
    It sort stack and sort the provided array.
    :param resName:
    :param meanGB:
    :param stdGB:
    :param meanPB:
    :param stdPB:
    :return:
    '''
    import numpy as np
    data = np.vstack((resName,(meanGB.astype(np.float)),stdGB.astype(np.float)))
    SortedData = data[:,np.argsort(data[1].astype(float))]
    return SortedData

def MMPBSAPlot(JobName, PrintLimit, SortedData,):
    '''

    :param JobName:
    :param PrintLimit:
    :param SortedData:
    :return:
    '''
    import numpy as np
    import matplotlib.pyplot as plt

    Residues = SortedData[0][:PrintLimit]
    #print (Residues)
    GB = SortedData[1][:PrintLimit].astype(np.float)
    #print (GB)
    #print (SortedData[2][:PrintLimit].astype(np.float))
    PB = SortedData[1][:PrintLimit].astype(np.float)

    N = len(GB)

    ind = np.arange(N)  # the x locations for the groups
    width = 0.35       # the width of the bars
    fig, ax = plt.subplots()
    rects1 = ax.bar(ind, GB, width, color='grey', yerr = SortedData[2][:PrintLimit].astype(np.float))

    rects2 = ax.bar(ind+width, PB, width, color='black', yerr = SortedData[2][:PrintLimit].astype(np.float))
    # add some text for labels, title and axes ticks
    ax.set_xlabel('Residues')
    ax.set_ylabel('Binding Free Energy Contribution kcal/mol')
    ax.set_title(JobName)
    ax.set_xticks(ind+width)
    ax.set_xticklabels( (SortedData[0]), rotation='vertical',size=7 )
    ax.legend( (rects1[0], rects2[0]), ('GB', 'PB'),loc=4 )


    plt.axhline(0, color='black')

    plt.savefig(JobName+'.pdf',format='pdf')
   # plt.show()



def nres(prmtop):
    '''

    :param prmtop: Reads the prmtop file to detect number of residues.
    :return: return number of residues.
    '''


    fin = open(prmtop,mode='r')
    for i in range(0,8):
        line = fin.readline()
        b = line.split()
        #print(line)
    c = int(b[1])
    return c

parser = argparse.ArgumentParser(description='This program is to analyse Per Residue Decomposition Resutls \
                                            from MMPBSA.py program of AMBER MD package. Author Tumbi Khaled.')
parser.add_argument('-j', '--jobName', help='Job Name', default='Results_MMPBSA_data_Analysis')
parser.add_argument('-lp', '--ligandPrmtop', help='Lignad prmtop file', required=True)
parser.add_argument('-rp', '--receptorPrmtop', help='Receptopn prmtop file', required=True)
parser.add_argument('-l', '--printLimit', help='Enter the number of residues of receptor to be plotted in results graph.')
parser.add_argument('-t', '--jobType', help='Define job type here. 1: Analyse only receptor. (Useful in Protein Ligand Interaction)'
                                            '2: Analyse only ligand. (Usful in Protein Protein interactions)'
                                            '3: Analsyse both receptor and ligand.', default=3)
args = parser.parse_args()

rres = nres(args.ligandPrmtop)
lres = nres(args.receptorPrmtop)
d = []

fname = []
for dirpath, dirnames, filenames in os.walk("."):
    for filename in [ f for f in filenames if f == "FINAL_DECOMP_MMPBSA.dat"]:
        fname.append(os.path.join(dirpath,filename))
lig_output = open(args.jobName+'_lig-Data.csv', mode='w')
pro_output = open(args.jobName+'_pro-Data.csv', mode='w')

print('Started your calculations. Now wait and watch.')
print (fname)
lGBTm,rGBTm,lresName, rresName,  lGBTstd, rGBTstd = xtractor (fname)

if args.jobType == str(1):
    SortedData = sorter(rresName,rGBTm,rGBTstd,rPBTm,rPBTstd)
    MMPBSAPlot(args.jobName, args.printLimit, SortedData)

if args.jobType == str(2):
    SortedData = sorterlig(lresName, lGBTm, lGBTstd, lPBTm,lPBTstd)
    MMPBSAPlot(args.jobName, lres, SortedData)

if args.jobType == str(3):
    SortedData = sorter(rresName,rGBTm,rGBTstd)
    np.savetxt(pro_output, SortedData, delimiter=',', fmt="%s")
    #print (SortedData)
    MMPBSAPlot(str(args.jobName)+'_(Receptor)', int(args.printLimit), SortedData)

    SortedData = sorterlig(lresName, lGBTm, lGBTstd)
    np.savetxt(lig_output, SortedData, delimiter=',', fmt="%s")
    MMPBSAPlot(str(args.jobName)+'_(Ligand)', lres, SortedData)


print('#####################################################################')
print('\n\n')
print('Results saved as pdf file. In the current directory.')
print('\n\n###################################################################')
