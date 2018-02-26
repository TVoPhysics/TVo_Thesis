import trace_mode
import pykat

def find_all(kat_findmode, q_sqz_x, q_sqz_y):
	#### Turn off all the cav commands except SRCY to extract the eigenmode at the BS
	[srcx,srcy] = trace_mode.from_cav_to_BS(kat_findmode,'cavSRY')

	#### Turn off all the cav commands except XARM to extract the eigenmode at the BS
	[ARMx,ARMy] = trace_mode.from_cav_to_BS(kat_findmode,'cavXARM')

	#### Turn off all the cav commands except OMC to extract the eigenmode at the BS
	[OMCx,OMCy] = trace_mode.from_cav_to_BS(kat_findmode,'cavOMC')

	#### Turn off all the cav commands except FC to extract the eigenmode at the BS
	[FCx,FCy] = trace_mode.from_cav_to_BS(kat_findmode,'cavFC')

	#### Get the sqz mode at the BS, this is more difficult because in Finesse, the sqzer is not a cavity.
	[q_sqz_BS_x,q_sqz_BS_y] = trace_mode.get_sqz_mode_BS(
			kat_findmode,q_sqz_x,q_sqz_y)

	strings = ['SRCx', 'SRCy',
			  'ARMx', 'ARMy',
			  'OMCx', 'OMCy',
			  'FCx', 'FCy',
			  'q_sqz_BS_x','q_sqz_BS_y']

	values = [srcx,srcy,
			 ARMx,ARMy,
			 OMCx,OMCy,
			 FCx,FCy,
			 q_sqz_BS_x,q_sqz_BS_y]

	strings_values = zip(strings, values)

	dictionary = {}
	for strings, values in strings_values:
		dictionary[strings] = values



	print("Overlap between Sqz and ARM = " + str(round(pykat.BeamParam.overlap(
					dictionary['q_sqz_BS_x'],dictionary['ARMx']) ,3)))

	print("Overlap between Sqz and FC = "  + str(round(pykat.BeamParam.overlap(
					dictionary['q_sqz_BS_x'],dictionary['FCx'])  ,3)))

	print("Overlap between Sqz and OMC = " + str(round(pykat.BeamParam.overlap(
					dictionary['q_sqz_BS_x'],dictionary['OMCx']) ,3)))

	print("Overlap between Sqz and SRC = " + str(round(pykat.BeamParam.overlap(
					dictionary['q_sqz_BS_x'],dictionary['SRCx']) ,3)))

	print("Overlap between ARM and SRC = " + str(round(pykat.BeamParam.overlap(
					dictionary['ARMx'],dictionary['SRCx']) ,3)))

	return dictionary
