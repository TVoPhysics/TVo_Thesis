import pykat

def phase_space(q_list_x,q_list_y):
	wx = []
	wy = []
	sx = []
	sy = []
	for i in range(len(q_list_x)):
	    qx = pykat.BeamParam(q=q_list_x[i])
	    qy = pykat.BeamParam(q=q_list_y[i])

	    #cuvature in microdipoters
	    sx.append( -(1/qx.curvature())*1e6)
	    sy.append( -(1/qy.curvature())*1e6)
	    
	    #curvature in mm
	    wx.append( 1e3*qx.beamsize())
	    wy.append( 1e3*qy.beamsize())

	return wx,wy,sx,sy
