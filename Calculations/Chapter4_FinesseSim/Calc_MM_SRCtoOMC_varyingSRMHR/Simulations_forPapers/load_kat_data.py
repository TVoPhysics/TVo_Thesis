import numpy

def load(data):
	freq = []
	darm = []
	
	for i in range(len(data)):
		freq.append(data[i][0])
		darm.append(data[i][1])
	
	return freq,darm
