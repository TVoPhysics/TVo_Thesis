#Commands to get QM-lim sensitivity
darm_commands = """
tf sus 1 0 p $mech_fres $mech_Q
const mech_fres 1  # Approx. resonance frequency
const mech_Q    1M # Guess for suspension Q factor
# Differentially modulate the strain in the arms
fsig darm  LXarm 1 0 1
fsig darm2 LYarm 1 180 1

qnoisedS NSR_with_RP 1 $fs nOMC_AROC_trans
qshotS NSR_without_RP 1 $fs nOMC_AROC_trans
pd1 pdAS $fs nOMC_AROC_trans

xaxis darm f log 5 5k 100
yaxis lin re:im
retrace off
"""

commands = """
tf sus 1 0 p $mech_fres $mech_Q
const mech_fres 1  # Approx. resonance frequency
const mech_Q    1M # Guess for suspension Q factor
# Differentially modulate the strain in the arms
fsig darm  LXarm 1 0 1
fsig darm2 LYarm 1 180 1
#qnoisedS NSR_with_RP 1 $fs nOMC_AROC_trans
#xaxis darm f log 1 5k 100
noxaxis
retrace off
"""
###Add amplitude detectors to get Higher-Order mordal content within the cavity at the carrier freq
amplitude_detectors = """
ad nSRMHRaTEM00 0 0 0 nSRMHRa
ad nSRMHRaTEM01 0 1 0 nSRMHRa
ad nSRMHRaTEM02 0 2 0 nSRMHRa
#ad nSRMHRaTEM03 0 3 0 nSRMHRa
#ad nSRMHRaTEM04 0 4 0 nSRMHRa
"""
###Add photodetectors to get gaussian beam parameters right outside of SRC after the anti-thermal lens correction
pds = """
bp SRCoutx x q nIBAin
bp SRCouty y q nIBAin

bp SRMYqx x q nSRMHRa
bp SRMYqy y q nSRMHRa

bp ITMXqx x q nITMX2
bp ITMXqy y q nITMX2

bp ITMYqx x q nITMY2
bp ITMYqy y q nITMY2

bp OMCqx x q nOMC_HROC_refl
bp OMCqy y q nOMC_HROC_refl

bp OFIqx x q nIBAin
bp OFIqy y q nIBAin

"""
quantum_detectors= """
qnoised NSR_shot_rad 1 $fs nOMC_AROC_trans

qd qdA 0 0 nOMC_AROC_trans
qd qdP 0 90 nOMC_AROC_trans

sd sd00 0 0 0 nOMC_AROC_trans
sd sd01 0 0 1 nOMC_AROC_trans
sd sd02 0 0 2 nOMC_AROC_trans

sd* sd00m 0 0 0 nOMC_AROC_trans
sd* sd01m 0 0 1 nOMC_AROC_trans
sd* sd02m 0 0 2 nOMC_AROC_trans

qnoised PDquantumnoise 1 $fs nOMC_AROC_trans
qshot PDshotnoise 1 $fs nOMC_AROC_trans


sd sd00OFI 0 0 0 nIMFC1
sd sd01OFI 0 0 1 nIMFC1
sd sd02OFI 0 0 2 nIMFC1

pd1 signal $fs nOMC_AROC_trans

ad OMCoutTEM00 0 0 0 nOMC_AROC_trans
ad OMCoutTEM01 0 1 0 nOMC_AROC_trans
ad OMCoutTEM02 0 2 0 nOMC_AROC_trans
"""

add_squeezing="""
sq sqz 0 10 0 nsqz
"""
