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

def get_sqz_mode_BS(kat,q_in_x,q_iny):
    
    #Turn the sqz_q from simulation into a pykat object to easily calc parameters
    q_sqz_x_new = pykat.BeamParam(q=q_in_x)
    q_sqz_y_new = pykat.BeamParam(q=q_iny)
    
    q_sqz_in = "gauss cavSQZ sqz nsqz " +  str(q_sqz_x_new.w0) + " " + str(-q_sqz_x_new.z) +\
    " " + str(q_sqz_y_new.w0) + " " + str(-q_sqz_y_new.z)
        
    #copy kat code of the simulation
    kat_sqz = kat.deepcopy()
    kat_sqz.parseCommands(q_sqz_in)

    #### Get the sqz mode at the beamsplitter by turning off all the cavities so all the nodes refer to the sqzer
    for cav in kat_sqz.getAll(pykat.commands.cavity):
        cav.enabled = False
    
    out_get = kat_sqz.run(getTraceData = True)
    q_sqz_BS_x = out_get[1][0]['nPRBS'][0]._BeamParam__q
    q_sqz_BS_y = out_get[1][0]['nPRBS'][1]._BeamParam__q

    return q_sqz_BS_x, q_sqz_BS_y
