#Pulls delta PSI values from 2 filtered MISO outputs.  Useful for correlating those values.  Takes a list of events to query.  Events must be in both MISO outputs.

#Usage: python getdeltaPSIs.py <eventfile> <FilteredMISOfile1> <FilteredMISOfile2> <EventType> <Outputfile>

import sys

eventsfh = open(sys.argv[1],'r')
N2Amisocomparisonfh = open(sys.argv[2], 'r')
CADmisocomparisonfh = open(sys.argv[3], 'r')
outfh = open(sys.argv[5], 'w')

eventtype = sys.argv[4]
events=[]
N2Amisocomparisons=[]
CADmisocomparisons=[]
N2Adeltapsis = []
CADdeltapsis = []
Combineddeltapsis = []

for line in eventsfh:
    line = line.strip()
    events.append(line)

for line in N2Amisocomparisonfh:
    line=line.strip()
    line=line.split('\t')
    N2Amisocomparisons.append(line)

for line in CADmisocomparisonfh:
    line = line.strip()
    line = line.split('\t')
    CADmisocomparisons.append(line)

    
#Pull deltaPSI values from N2Amisocomparison file where gene matches gene in query list

for event in events:
    for misocomparison in N2Amisocomparisons:
        if event != 'event_name':
            if event == misocomparison[0]:
                N2Adeltapsi = [misocomparison[0], misocomparison[7]]
                N2Adeltapsis.append(N2Adeltapsi)
                
#Pull deltaPSI values from CADmisocomparison file where gene matches gene in query list
                
for event in events:
    for misocomparison in CADmisocomparisons:
        if event != 'event_name':
            if event == misocomparison[0]:
                CADdeltapsi = [misocomparison[0], misocomparison[7]]
                CADdeltapsis.append(CADdeltapsi)

#Combine two lists.  Add Event type identifier (sys.argv[4]).

for N2Adeltapsi in N2Adeltapsis:
    for CADdeltapsi in CADdeltapsis:
        if N2Adeltapsi[0] == CADdeltapsi[0]:
            Combineddeltapsis.append([N2Adeltapsi[0], N2Adeltapsi[1], CADdeltapsi[1], eventtype])

for Combineddeltapsi in Combineddeltapsis:
    Combineddeltapsi = ('\t').join(Combineddeltapsi)
    outfh.write(Combineddeltapsi + '\n')
