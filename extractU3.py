print 'Collecting data from .odb file'
#Type abaqus python extractU3.py of the script to launch it on command line.  
######################################################################################################
#Input parameters
input='File.odb' # the name of the database
disp='U3' # wanted displacement component of the node 
output_file='U3' # outputfile with extension _ added in the program (*.rpt file created)
######################################################################################################
from odbAccess import *
from abaqusConstants import *
from odbSection import *
steps=64
Nodes=[100901401, 100899438, 100893471, 100884500, 100871524, 100854543, 100834556, 100810564, 100784567, 100755564, 100723556, 100690543, 100654524, 100617500, 100579471, 100540438, 100501401, 100462360, 100423315, 100385268, 100348218, 100312165, 100279112, 100247056, 100218001, 100191946, 100167890, 100147837, 100130784, 100117734, 100108687, 100102642, 100100601, 100102564, 100108531, 100117502, 100130478, 100147459, 100167446, 100191438, 100217435, 100246438, 100276446, 100311459, 100347478, 100384502, 100422531, 100461564, 100500601, 100539642, 100578687, 100616734, 100653784, 100689837, 100722890, 100754946, 100784001, 100810056, 100834112, 100854165, 100871218, 100884268, 100893315, 100899360]
output=input
print output
odb = openOdb(output)  
step1 = odb.steps.values()[0] 
for k in range(len(Nodes)):
     output=output_file+'_'+str(k+1)+'.txt'
     nr=str(Nodes[k])
     modified='Node PART-1-1.'+nr
     region = step1.historyRegions[modified]
     u3Data = region.historyOutputs[disp].data
     dispFile = open(output,'w')
     for time,u3Disp in u3Data:
            dispFile.write('%10.4E\n' % (u3Disp))
     dispFile.close()
odb.close()    
print 'Finish'  
