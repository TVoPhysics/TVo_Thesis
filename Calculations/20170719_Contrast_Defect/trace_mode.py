import pykat

def from_cav_to_BS(kat,cavity):
	kat.maxtem = 2
	for cav in kat.getAll(pykat.commands.cavity):
		cav.enabled = False

	katcavstr='katcav=kat.'+cavity+'.enabled=True'
	exec(katcavstr)
	#       katcav.enabled = True

	out_get = kat.run(getTraceData = True)
	qx = out_get[1][0]['nPRBS'][0]._BeamParam__q
	qy = out_get[1][0]['nPRBS'][1]._BeamParam__q

	return qx, qy


def from_cav_to_node(kat,cavity,node):
	kat.maxtem = 2
	for cav in kat.getAll(pykat.commands.cavity):
		cav.enabled = False

	katcavstr='katcav=kat.'+cavity+'.enabled=True'
	exec(katcavstr)
	#       katcav.enabled = True


	out_get = kat.run(getTraceData = True)

	try:
		qx = out_get[1][0][str(node)][0]._BeamParam__q
		qy = out_get[1][0][str(node)][1]._BeamParam__q
	except:
		print("Could not extract node from output")

	return qx, qy

