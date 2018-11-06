'''

[2] LIGO-D0902216-v8
[3] galaxy.ligo.caltech.edu/optics
[4] LLO alog 27208
[5] LIGO-E1200438-v1
[6] LIGO-E1200464-v1
[7] LIGO-E1200274-v3
[8] LIGO-E1200607-v1
[9] LIGO-E1200502-v2
[10] LIGO-E1200210-v1
[11] LIGO-E1200706-v1
[12] LIGO-E1300397-v1
[13] LIGO-E1200524-v1
[14] T1300954
'''

import pykat
import pykat.ifo.aligo as aligo
import pykat.exceptions as pkex
import numpy as np


def make_LHO(verbose=False):
    '''
    Script creating an Advanced LIGO Hanford model.
    
    Checked:
        SRM -> OMC Path
        SQZ -> SRM Path
        
    Still need to check all the IFO parameters and input 
    '''
    kat = aligo.make_kat("design_with_IMC_HAM2_FI_OMC", keepComments=True, preserveConstants=True)

    ###########################################################################
    # Constants
    ###########################################################################

    # Refractive indices
    n_silica = 1.44963098985906   # From other kat-files. Consistent with [7].

    ###########################################################################
    # Input
    ###########################################################################

    # Input power [W]
    Pin = 22.7  # Fix!
    # Modulation frequencies [Hz], from [1, 4].
    f1 = 9099055.0
    f2 = 45495275.0
    
    # Modulation depths [rad], from [1], derived from [4].
    midx1 = 0.13
    midx2 = 0.139985
    
    if (5 * f1 != f2):
        print(" ** Warning: modulation frequencies do not match: 5*f1!=f2")

    # Set values
    # -------------
    kat.L0.P = Pin
    kat.f1.value = f1
    kat.IFO.f1 = f1
    kat.mod1.midx = midx1
    kat.f2.value = f2
    kat.IFO.f2 = f2
    kat.mod2.midx = midx2
    
    ###########################################################################
    # IMC to PRM PATH
    ###########################################################################
        
    kat.sIM2_FI.L  = 260e-3   # T1700227
    kat.sFI_IM3.L  = 910e-3   # T1700227
    kat.sIM3_IM4.L = 1175e-3  # T1700227
    kat.sPRCin.L   = 413e-3   # T1700227
    
    ###########################################################################
    # PRC
    ###########################################################################

    # Lengths [m] from [1, 2].
    # Radii of curvature [m] from [1, 3].
    # Power reflections, transmissions and losses from [1, 3].
    # Mass from [3].
    # Thickness from [5].

    # PRM (PRM 04)
    # -------------
    # Radius of curvature
    RoC_prm =  10.948                        # [3] -1 sign because curvature is on side 2
    T_prm_hr = 0.031                         # [3]
    L_prm_hr = (8 + 0.5)*1e-6                # Scatter + absorbtion
    R_prm_hr = 1.0 - T_prm_hr - L_prm_hr     
    
    R_prm_ar = 0                       
    L_prm_ar = (35+4.5)*1e-6                 # [3]
    T_prm_ar = 1.0 - R_prm_ar - L_prm_ar     
    
    th_prm = 74.9e-3       
    M_prm = 2.890          

    # PR2 (PR2 04)
    # -------------
    # Radius of curvature
    RoC_pr2  = -4.543              # [3]
    T_pr2_hr = 225.0e-6            # [3] Average
    L_pr2_hr = (10.0 + 0.7)*1e-6   # [3] Scatter + absorbtion
    R_pr2_hr = 1.0 - T_pr2_hr - L_pr2_hr
    
    L_pr2_ar = 8e-6                # [3]
    M_pr2    = 2.894               # [3]
    

    # PR3 (PR3 01)
    # -------------
    # Radius of curvature
    RoC_pr3 = 36.02               # [3]
    T_pr3_hr = 3e-6               # [3]
    L_pr3_hr = (16.5 + 0.45)*1e-6 # [3]
    R_pr3_hr = 1.0 - T_pr3_hr - L_pr3_hr
    R_pr3_ar = 0                  # Not measured and unused. [1] using 20 ppm.
    M_pr3 = 12.145                # [3]

    # Spaces
    # -------------
    L_prm_pr2 = 16.6107      # PRM to PR2
    L_pr2_pr3 = 16.1648      # PR2 to PR3
    L_pr3_bs  = 19.5380      # PR3 to BS

    # Derived lpr length
    lpr = L_prm_pr2 + L_pr2_pr3 + L_pr3_bs

    # Set values to kat-object
    # -------------
    # Lengths
    kat.lp1.L = L_prm_pr2
    kat.lp2.L = L_pr2_pr3
    kat.lp3.L = L_pr3_bs
    # Substrate thicknesses
    kat.sPRMsub1.L = th_prm
    # RoCs
    kat.PRM.Rc = RoC_prm
    kat.PR2.Rc = RoC_pr2
    kat.PR3.Rc = RoC_pr3 
    # Refl, tran, and loss
    kat.PRM.R = R_prm_hr
    kat.PRM.T = T_prm_hr
    kat.PRM.L = L_prm_hr
    kat.PRMAR.R = 0    # Adding reflection to loss due to the wedge.
    kat.PRMAR.T = T_prm_ar
    kat.PRMAR.L = L_prm_ar + R_prm_ar  # See above.

    kat.PR2.R = R_pr2_hr
    kat.PR2.T = T_pr2_hr
    kat.PR2.L = L_pr2_hr

    kat.PR3.R = R_pr3_hr
    kat.PR3.T = T_pr3_hr
    kat.PR3.L = L_pr3_hr

    ###########################################################################
    # Beam splitter (BS 01)
    ###########################################################################

    # Power reflections, transmissions and losses from [1, 3]. 
    
    # R and T for highly-reflective surface
    T_bs_hr = 0.5
    L_bs_hr = (14.0 + 0.3)*1e-6           # [1] using 30 ppm 
    R_bs_hr = 1.0 - T_bs_hr - L_bs_hr    # = 0.4999914
    # R and T for anti-reflective surface
    R_bs_ar = 50e-6                      
    L_bs_ar = 1.4e-6                     # Not set in [1]
    T_bs_ar = 1.0 - R_bs_ar - L_bs_ar    # = 0.9999683
    # Mass
    M_bs = 14.251                        # Unused. [1] using 14.2 kg.
    # Thickness
    th_bs = 60.41e-3                     # Not set in [1].

    # Angle of incidence HR [deg], seen from laser side. 
    alpha_in_hr =  45
    # Angle of refraction HR
    alpha_out_hr = np.arcsin(np.sin(alpha_in_hr*np.pi/180.0)/n_silica)*180.0/np.pi
    # = 29.195033
    # Angle of incidence AR
    alpha_in_ar = alpha_out_hr

    # Thickness seen from beam
    th_bs_beam = th_bs/np.cos(alpha_out_hr*np.pi/180.0) # = 0.06887


    # Set values to kat-object
    # -------------

    # Refl, Trans, and loss.
    #  - HR
    kat.BS.R = R_bs_hr
    kat.BS.T = T_bs_hr
    kat.BS.L = L_bs_hr
    #  - AR1
    kat.BSAR1.R = 0
    kat.BSAR1.T = T_bs_ar
    kat.BSAR1.L = L_bs_ar + R_bs_ar
    #  - AR2
    kat.BSAR2.R = 0
    kat.BSAR2.T = T_bs_ar
    kat.BSAR2.L = L_bs_ar + R_bs_ar

    # Angles of incidence
    kat.BS.alpha = alpha_in_hr
    kat.BSAR1.alpha = alpha_in_ar   # Check! Signs swapped in other kat-files.
    kat.BSAR2.alpha = -alpha_in_ar
    # Propagation distance through substrate
    kat.BSsub1.L = th_bs_beam
    kat.BSsub2.L = th_bs_beam


    ###########################################################################
    # Y - arm 
    ###########################################################################

    # Lengths from [1, 2, 7]. Positions of ITMY AR and CP1 X and Y from [7].
    # Refl, Trans, loss, masses [kg], RoCs [m] from [3].
    # ITM Thickness from [8], averaged over the 4 values.
    # Compensation plates aren't included as discrete components, they are
    # included in the model as just another part of the "ITM lens"


    # CPY (CP 05) [3]
    # -------------
    # Thickness
    th_cpy = 0.10002                
    # Loss
    L_cpy1 = 0.4e-6
    L_cpy2 = 0.4e-6
    # Refl
    R_cpy1 = 27.2e-6
    R_cpy2 = 5e-6
    Rc_cpy = -1392000

    # Input test mass Y (ITM 11)
    # -------------

    # Thermal and substrate lenses in ITMY. From other kat-file. Check!
    TLY_f = np.inf           # Thermal lens ITMY, focal length
    SLY_f = 572117           # Cold substrate lens, focal length [14]
    
    # Refl, Trans, and loss.
    #  - HR
    T_itmy_hr = 0.0142                     # [3]
    L_itmy_hr = (12.0 + 0.4)*1e-6          # Scatter + absorbtion. [1] has modified
                                           # this to 65 ppm to get correct PRC gain. 
    R_itmy_hr = 1.0-T_itmy_hr-L_itmy_hr    # = 0.98518570
    #  - AR
    R_itmy_ar = 294e-6                     # [3]
    L_itmy_ar = 1e-6                       # [3]
    T_itmy_ar = 1.0-R_itmy_ar-L_itmy_ar    
    # Thickness
    th_itmy = 0.199                        # Not set in [1].
    # Radius of curvature
    RoC_itmy_hr = -1939.2
    # Mass
    M_itmy = 39.636                        # [3]

    # End test mass Y (ETM 16)
    # -------------

    # Refl, Trans, and loss.
    #  - HR
    T_etmy_hr = 3.8e-6                     # [3]
    L_etmy_hr = (12 + 0.2)*1e-6           # [3]
    
    #  - AR
    R_etmy_ar = 20e-6                       # [3]
    L_etmy_ar = 0                          # Unknown
    
    # Thickness
    th_etmy = 0.2    # Estimated from [2], and other fiducial line
                     # measurements. Consistent with other fiels.
    # Radius of curvature
    RoC_etmy_hr = 2246.9
    # Mass
    M_etmy = 39.564                        # [3]

    # Lengths [m]
    # -------------
    # BS HR to CPY1 (no CP in [1])
    L_bs_cpy = 4.8470                        
    # CPY2 to ITMY AR (no CP in [1])
    L_cpy_itmyAR = 0.0200                     
    # Optical distance from BS HR to ITMY AR. 
    L_bs_itmyAR = L_bs_cpy + th_cpy*n_silica + L_cpy_itmyAR  # = 5.0120
    # Arm length. ITMY HR to ETMY HR. 
    L_itmy_etmy = 3994.4850      # From [1, 2], [7] has 3 cm diff in ETMY position.

    # Derived
    # -------------
    # Optical distances
    L_bshr_itmyhr = L_bs_itmyAR + th_itmy*n_silica   # = 5.3009

    # Set values to kat-object
    # -------------

    # Set the lens, here we just set the cold state
    kat.ITMY_lens.f = SLY_f

    # ITMY
    kat.ITMY.R     = R_itmy_hr
    kat.ITMY.T     = None
    kat.ITMY.L     = L_itmy_hr
    kat.ITMYAR.R   = None # Putting Refl on loss for AR due to wedge as well as CP losses and refl
    kat.ITMYAR.L   = L_itmy_ar + R_itmy_ar + R_cpy1 + L_cpy1 + R_cpy2 + L_cpy2
    kat.ITMYAR.T   = 1 - kat.ITMYAR.L
    kat.ITMYsub.L  = th_itmy
    kat.ITMY.Rc    = RoC_itmy_hr
    kat.ITMY.mass  = M_itmy

    # ETMY
    kat.ETMY.R     = None
    kat.ETMY.T     = T_etmy_hr
    kat.ETMY.L     = L_etmy_hr
    kat.ETMYAR.R   = 0 # Putting Refl on loss for AR due to wedge.
    kat.ETMYAR.T   = None
    kat.ETMYAR.L   = 1 - R_etmy_ar - L_etmy_ar
    kat.ETMYsub.L  = th_etmy
    kat.ETMY.Rc    = RoC_etmy_hr
    kat.ETMY.mass  = M_etmy

    # Lengths
    kat.ly1.L = L_bs_itmyAR
    kat.LY.L  = L_itmy_etmy

    ###########################################################################
    # X - arm
    ###########################################################################


    # CPX (CP 01)
    # -------
    R_cpx1 = 10.2e-6
    L_cpx1 = 0.5e-6

    R_cpx2 = 43.7e-6
    L_cpx2 = 0.3e-6
    # Thickness, average from [3]
    th_cpx = 0.100
    Rc_cpx = 664100

    # ITMX (ITM 07)
    # --------
    # Refl, Trans, and loss
    T_itmx_hr = 0.0148                    # 1.48 ppm in [2] (missing %?)
    L_itmx_hr = (5.7 + 0.5)*1e-6         # [1] using 65 ppm, adjusted for cavity gains.
    R_itmx_hr = 1.0-T_itmx_hr-L_itmx_hr   # = 0.9851896

    R_itmx_ar = 59e-6                    # [1] using 16.4 ppm
    L_itmx_ar = 1.1e-6                    # Not set in [1]
    T_itmx_ar = 1.0-R_itmx_ar-L_itmx_ar   # = 0.9998355
    #  Thickness, average from [11]
    th_itmx = 0.19996                     # Not set in [1]
    # Radius of curvature
    RoC_itmx_hr = -1940.3
    # Mass
    M_itmx = 39.603                       # [1] using 40 kg      

    # ITM lens in cold state
    SLX_f = -310812                       # Constant ITMX substrate lens, focal length [14]

    # ETMX (ETM 13)
    # -----
    # Refl, trans, and loss
    T_etmx_hr = 3.9e-6                     
    L_etmx_hr = (21 + 0.2)*1e-6            # Scatter + absorbtion. [1] using 30 ppm for cavity gains. 

    R_etmx_ar = 8.9e-6                     
    L_etmx_ar = 0                          # Unspecified.
    
    th_etmx = 0.19985                      # Not set in [1]
    # Radius of curvature
    RoC_etmx_hr = 2244.2
    # Mass
    M_etmx = 39.630                       # [1] using 40 kg      

    # Lengths
    # -------
    # Position of BSar from [2], position of CPX1 from [2,7].
    L_bs_cpx = 4.8295
    # Position of CPX2 from [2,7]. Position of ITMX AR from [7]. Consistent with [2]. 
    L_cpx_itmx = 0.0200
    # Positions from [2].
    L_itmx_etmx = 3994.4850

    # Derived
    # -------
    # Optical distance from BS AR to ITMX AR. 
    L_bsar_itmxar = L_bs_cpx + th_cpx*n_silica + L_cpx_itmx 
    L_bshr_itmxhr = L_bsar_itmxar + (th_itmx + th_bs_beam)*n_silica

    # Set values to kat-object
    # ------------------------
    # Thermal and substrate lenses
    kat.ITMX_lens.f = SLX_f

    # Lengths
    kat.lx1.L = L_bsar_itmxar
    kat.LX.L  = L_itmx_etmx

    # ITMX
    kat.ITMX.R   = None
    kat.ITMX.T   = T_itmx_hr
    kat.ITMX.L   = L_itmx_hr
    kat.ITMXAR.R = None # Putting Refl on loss for AR due to wedge along with CP losses and refl
    kat.ITMXAR.L = L_itmx_ar + R_itmx_ar + R_cpx1 + L_cpx1 + R_cpx2 + L_cpx2
    kat.ITMXAR.T = 1 - kat.ITMXAR.L
    
    # Thickness
    kat.ITMXsub.L = th_itmx
    # RoC
    kat.ITMX.Rc = RoC_itmx_hr

    # ETMX
    kat.ETMX.R = None
    kat.ETMX.T = T_etmx_hr
    kat.ETMX.L = L_etmx_hr
    # Putting Refl on loss for AR due to wedge.
    kat.ETMXAR.R = None
    kat.ETMXAR.L = L_etmx_ar + R_etmx_ar
    kat.ETMXAR.T = 1 - kat.ETMXAR.L
    # Thickness
    kat.ETMXsub.L = th_etmx
    # RoC
    kat.ETMX.Rc = RoC_etmx_hr
    # Refl, trans and loss

    ###########################################################################
    # SRC
    ###########################################################################

    # Lengths from [2]

    # Lengths
    # -------
    # BSars to SR3
    L_bs_sr3 = 19.3659          
    # SR3 to SR2
    L_sr3_sr2 = 15.4435
    # SR2 to SRM
    L_sr2_srm = 15.7562         
    # Optical length from BS_HR to SRM_HR
    L_bshr_srm = (th_bs_beam*n_silica +
                  L_bs_sr3 + L_sr3_sr2 +
                  L_sr2_srm)               # = 50.6654 m

    # SRM (SRM 06)
    T_srm_hr = 0.3234                  # [3]
    
    #T_srm_hr = 0.3
    L_srm_hr = (8.5 + 0.7)*1e-6        # [3] Scatter + absorbtion. 
    R_srm_hr = 1.0 - T_srm_hr - L_srm_hr   

    R_srm_ar = 0                   
    L_srm_ar = (30+0.6) * 1e-6         # [3]
    T_srm_ar = 1.0-R_srm_ar-L_srm_ar
    
    th_srm = 0.074                     # [3]
    RoC_srm = -5.677                   # [3]
    M_srm = 2.899                      # [3]

    # SR2 (SR2 03)
    T_sr2_hr = 7.5e-6                  # [3]
    L_sr2_hr = (12+0.5)*1e-6           # [3]
    R_sr2_hr = 1.0-T_sr2_hr-L_sr2_hr   
    
    R_sr2_ar = 0                       # [3]  
    L_sr2_ar = 170e-6                  # [3]  
    T_sr2_ar = 1.0-R_sr2_ar-L_sr2_ar   
    
    M_sr2 = 2.888                      # [3]
    RoC_sr2 = -6.424                   # [3]


    # SR3 (SR3 02)
    T_sr3_hr = 3.5e-6
    L_sr3_hr = (21.0+0.3)*1e-6         # Scatter+absorbtion
    R_sr3_hr = 1.0-T_sr3_hr-L_sr3_hr   
    
    R_sr3_ar = 0                       # Not listed
    L_sr3_ar = 0                       # Not listed
    T_sr3_ar = 1.0-R_sr3_ar-L_sr3_ar   
    
    RoC_sr3 = 36.013                    
    M_sr3 = 12.069                     
    
    ###########################################################################
    # SRC to OMC
    ########################################################################### 
    
    SRM2OM1 = 3.45 # Zemax, Corey
    
    kat.sSRM_FI.L = 0.9046 # Distances from Zemax SRM to polarizer: Corey email
    kat.sFI_OM1.L = SRM2OM1 - kat.sSRM_FI.L
    
    # 1.39m Measured by Dan, Danny, TVo, Terra ~ 14th April 2018
    kat.sOM1_OM2.L = 1.39
    
    # Measured by Shelia ~17th April 2018
    kat.sOM2_OM3.L = 0.63
    kat.sOM3_OMC.L = 0.117 + 0.2
    
    # See T1200410-v2
    kat.OM1.Rc = 4.6  # E1100056-v2-02
    kat.OM1.setRTL(1-800e-6, 800e-6, 0)

    kat.OM2.Rc = 1.7  #   E1100056-v2-01
    kat.OM2.setRTL(1-800e-6, 800e-6, 0)
    
    kat.OM3.Rc = np.inf # E1000457-v1
    kat.OM3.setRTL(0.99, 0.01,0)
    
    # Add refl path used for tracing or outputs
    # Add AS AIR path for tracing or outputs, about guessed about 2.5m from OM1
    # as it is on the ISCT6 table that moves
    kat.parse("""
    s sOMC_REFL 2 nOMC_ICb nOMC_REFL 
    
    s sOM1_sub 0.01 1.45 nOM1c nOM1d
    m OM1AR 0 1 0 nOM1d nOM1e
    s sOM1_AS_AIR 2.5 nOM1e nAS_AIR
    """)
    kat.sOM1_sub.L = 0.01
    kat.sOM1_sub.n = n_silica
    
    ###########################################################################
    # Squeezer
    ########################################################################### 
    
    # This is a detailed path of the OPO to the OFI. This detail is only needed
    # for mode matching calculations. Once this is finalised we can simplify this
    # path with an option to this make_LHO function. This doesn't include the
    # actual OPO cavity but just the M1 cavity mirror where we start the trace from
    # as we know the cavity waist position from here between the two flat M1 and M2
    # mirrors which are 11cm apart.
    
    kat.parse("""
    %%% FTblock squeezer
    ###########################################################################
    sqz sqz 0 0 0 nSQZ
    #l sqz_clf 0 3M nSQZ 
    bs M1_OPO 0 1 0 6 nSQZ dump nM1_OPOc dump
    s subM1_OPO 6.35 $nsilica nM1_OPOc nM1_OPO_ARa
    bs M1_OPO_AR 0 1 0 6 nM1_OPO_ARa dump nM1_OPO_ARc dump
    s lM1_OPO_EDGE 44m nM1_OPO_ARc nOPO_EDGEa
    m OPO_EDGE 0 1 0 nOPO_EDGEa nOPO_EDGEb
    
    # edge of OPO block
    s lsqz_lens1 0 nOPO_EDGEb nSQZa
    
    lens sqz_lens1 1 nSQZa nSQZb
    
    s lsqz_lens1_faraday 0 nSQZb nSQZc
    
    # 20 mm TGG faraday crystal
    m msqz_faraday_a 0 1 0 nSQZc nSQZd
    s lsqz_faraday 20m 1.95 nSQZd nSQZe
    m msqz_faraday_b 0 1 0 nSQZe nSQZf
    
    s lsqz_faraday_lens2 0 nSQZf nSQZg
    
    lens sqz_lens2 1 nSQZg nSQZh
    s lsqz_lens2_zm1 0 nSQZh nZM1a 

    bs1 ZM1 0 0 0 45 nZM1a nZM1b dump dump
    s lzm1_zm2 1.0 nZM1b nZM2a
    
    bs1 ZM2 0 0 0 45 nZM2a nZM2b dump dump
    s lzm2_OFI 0 nZM2b nFI2b
    maxtem 2
    
    ###########################################################################
    %%% FTend squeezer
    """, keepComments=True)

    # VOPO - M1 to OPO block edge details - TT1700104-v2
    kat.M1_OPO.alpha    = 6 #
    kat.M1_OPO_AR.alpha = np.rad2deg(np.arcsin(np.sin(np.deg2rad(kat.M1_OPO.alpha.value))/kat.ITMXsub.n.value)) # Set indicent angle on AR surface
    kat.subM1_OPO.L     = 6.35e-3 * np.cos(np.deg2rad(kat.M1_OPO.alpha.value))
    kat.lM1_OPO_EDGE.L  = 44e-3 * np.cos(np.deg2rad(kat.M1_OPO.alpha.value)) # 44mm from Sheila - AR to OPO edge
    
    # 55mm comes from the fact that the waist is halfway between
    # M1 and M2 (separated by 110mm) and goign into M1 we are converging
    # to that waist. Waist sizes from TT1700104-v2 table 1
    # Mode slightly astigmatic coming out from OPO
    kat.sqz.nSQZ.qx = pykat.BeamParam(w0=221.2e-6, z=-55e-3)
    kat.sqz.nSQZ.qy = pykat.BeamParam(w0=206.9e-6, z=-55e-3)
    
    kat.sqz_lens1.f = 111e-3 # E1600300-v1–02
    kat.sqz_lens2.f = 334e-3 # E1600300-v1–05
    
    # Measured by Nutsinee 17th April 2018, Citing emailed picture
    # kat.lsqz_lens1.L = 31.4e-2 # from Nutsinee sketch but this is just to the edge of OPO case
    # kat.lsqz_lens1.L = 0.3531 # Lee's Alamode file, q reference inside OPO block so longer
    # Measurements done by Sheila (alog 40685) done in squeezer bay with class A ruler
    kat.data['OPO_block_to_lens1']   = 306e-3 # Measurements from Shelia (alog 40685)
    kat.data['OPO_lens1_to_faraday'] = 256e-3 - kat.lsqz_faraday.L.value/2 # Measured by Sheila 19/04/2018 - 10mm as crystal 20mm in faraday
    kat.data['OPO_faraday_to_lens2'] = 295e-3 - kat.lsqz_faraday.L.value/2 # Measured by Sheila 19/04/2018
    
    # kat.data['OPO_block_to_lens2']   = 856e-3 # Measurements from Shelia (alog 40685) 
    kat.data['OPO_block_to_lens2'] = kat.data['OPO_block_to_lens1'] + kat.data['OPO_lens1_to_faraday'] + kat.data['OPO_faraday_to_lens2'] + kat.lsqz_faraday.L.value
    
    # These are measurements taken by Nutsinee from lens2 stage to ZM1 with tape measure
    kat.data['OPO_block_to_ZM1']     = kat.data['OPO_block_to_lens2'] + 73e-3 + 132e-3
    kat.data['lens2_range']          = 42e-3/2 # Actuation range on OPO lens2 mount
        
    def shift_lenses(kat, lens1_dz, lens2_pos, lens2_dz=0):
        """
        Shifts lenses relative to those that have been measured in H1 OPO->ZM1 path
        
        lens1_dz - change in lens1 position in m, positive closer to OPO block
        lens2_pos - 0 to 1 value for position along translation stage
        lens2_dz - change in lens1 position in m, positive closer to OPO block
        """
        
        lens2_dz += (2*lens2_pos-1) * kat.data['lens2_range']
        
        kat.lsqz_lens1.L         = kat.data['OPO_block_to_lens1'] - lens1_dz
        kat.lsqz_lens1_faraday.L = kat.data['OPO_lens1_to_faraday'] + lens1_dz
        kat.lsqz_faraday_lens2.L = kat.data['OPO_faraday_to_lens2'] - lens2_dz
        kat.lsqz_lens2_zm1.L     = kat.data['OPO_block_to_ZM1'] - kat.lsqz_lens1.L.value - kat.lsqz_lens1_faraday.L.value - kat.lsqz_faraday_lens2.L.value
        
        if kat.lsqz_lens1.L < 0:         raise Exception("Uh oh")
        if kat.lsqz_lens1_faraday.L < 0: raise Exception("Uh oh")
        if kat.lsqz_faraday_lens2.L < 0: raise Exception("Uh oh")
        if kat.lsqz_lens2_zm1.L < 0:     raise Exception("Uh oh")
    
    # Store function for user to call
    kat.data['fn_OPO_shift_lenses'] = shift_lenses 

    # default 0.5 in translation stage, where distance measurements were done
    kat.data['fn_OPO_shift_lenses'](kat, 0, 0.5, 0)
    
    # From preliminary Solidworks file TVo ~18th April, could be different to as built
    kat.lzm1_zm2.L         = 2.260
    
    # From Zemax/Solidworks model distances, ZM2 -> OFI Polarizer: Cite Corey email
    kat.lzm2_OFI.L         = 120e-2
    
    
    ###########################################################################
    # Setting values to Kat-file
    ########################################################################### 
    
    # Lengths
    # ---
    kat.ls3.L = L_bs_sr3
    kat.ls2.L = L_sr3_sr2
    kat.ls1.L = L_sr2_srm

    # SRM
    # ---
    kat.SRM.R = R_srm_hr
    kat.SRM.T = T_srm_hr
    kat.SRM.L = L_srm_hr
    # Putting refl on loss due to wedge
    kat.SRMAR.R = 0
    kat.SRMAR.T = T_srm_ar
    kat.SRMAR.L = L_srm_ar + R_srm_ar
    # Tickness
    kat.SRMsub.L = th_srm
    # Mass
    kat.SRM.mass = M_srm
    # RoC
    kat.SRM.Rc = RoC_srm

    # SR2
    # ---
    kat.SR2.R = R_sr2_hr
    kat.SR2.T = T_sr2_hr
    kat.SR2.L = L_sr2_hr
    # Mass
    kat.SR2.mass = M_sr2
    # RoC
    kat.SR2.Rc = RoC_sr2

    # SR3
    # ---
    kat.SR3.R = R_sr3_hr
    kat.SR3.T = T_sr3_hr
    kat.SR3.L = L_sr3_hr
    # Mass
    kat.SR3.mass = M_sr3
    # RoC
    kat.SR3.Rc = RoC_sr3

    lprc = lpr + (L_bshr_itmyhr + L_bshr_itmxhr)/2
    
    if verbose:
        print('--------------------------------------')
        print('Lengths')
        print('--------------------------------------')
        print('lx (BShr to ITMXhr) = {:.4f} m'.format(L_bshr_itmxhr))
        print('ly (BShr to ITMYhr) = {:.4f} m'.format(L_bshr_itmyhr))
        print('Schnupp (lx - ly) = {:.4f} m'.format(L_bshr_itmxhr - L_bshr_itmyhr))
        print('lpr = {:.4f} m'.format(lpr))
        print('PRC = {:.4f} m'.format(lprc))
    
    return kat
