# make a subset of UNEAK data to keep for tutorial
nMarker = 200 # number of markers to keep

filePrefix = "HapMapMsiSubset"
outPrefix = "HapMapMsiSubset1"

with open(filePrefix + ".fas.txt", mode = 'rt') as incon:
    with open(outPrefix + ".fas.txt", mode = 'wt') as outcon:
        for i in range(nMarker * 4):
            outcon.write(incon.readline())

with open(filePrefix + ".hmc.txt", mode = 'rt') as incon:
    with open(outPrefix + ".hmc.txt", mode = 'wt') as outcon:
        for i in range(nMarker + 1):
            outcon.write(incon.readline())
