import numpy as np
from astropy.io import fits

#######################################################################################################################
## apjsaafad2t3_mrt.txt - main table
#######################################################################################################################
##   47- 50 F4.2   ---         Photz       ?="" Primary Photometric Redshift
##   52- 55 F4.2   ---       E_Photz       ?="" Upper Uncertainty in Photz
##   57- 60 F4.2   ---       e_Photz       ?="" Lower Uncertainty in Photz
##   62- 65 F4.2   ---         Photz2      ?="" Secondary Photometric Redshift
##   67- 70 F4.2   ---       E_Photz2      ?="" Upper Uncertainty in Photz2
##   72- 75 F4.2   ---       e_Photz2      ?="" Lower Uncertainty in Photz2
##   77- 82 F6.4   ---         Specz       ?="" Spectroscopic Redshift
##   84- 85 I2     ---       o_Specz       ?="" Number of Confirmed Cluster Members
##   87- 88 I2     ---         Richness    ?="" Richness for Primary Photometric Redshift
##   90- 90 I1     ---       e_Richness    ?="" Uncertainty in Richness for Primary Photometric Redshift
##   92- 93 I2     ---         Richness2   ?="" Richness for Secondary Photometric Redshift
##   95- 95 I1     ---       e_Richness2   ?="" Uncertainty in Richness for Secondary Photometric Redshift
##   97-100 F4.1   10+14Msun   M500        ?="" M500, enclosed mass
##  102-104 F3.1   10+14Msun E_M500        ?="" Upper Uncertainty in M500
##  106-108 F3.1   10+14Msun e_M500        ?="" Lower Uncertainty in M500

#######################################################################################################################
## apjsaafad2t4_mrt.txt
#######################################################################################################################
##   1- 14 A14    ---         Cluster     Cluster Name (1)
##  16- 17 I2     h           RAh         Hour of Right Ascension (J2000) 
##  19- 20 I2     min         RAm         Minute of Right Ascension (J2000) 
##  22- 25 F4.1   s           RAs         Second of Right Ascension (J2000) 
##  27- 27 A1     ---         DE-         Sign of the Declination (J2000)
##  28- 29 I2     deg         DEd         Degree of Declination (J2000) 
##  31- 32 I2     arcmin      DEm         Arcminute of Declination (J2000) 
##  34- 35 I2     arcsec      DEs         Arcsecond of Declination (J2000) 
##  37- 40 F4.2   ---         Peak        Peak Height (2)
##  42- 45 F4.1   ---         SNR         Signal-to-Noise
##  47- 50 F4.2   ---         Photz       ?="" Primary Photometric Redshift (3)
##  52- 55 F4.2   ---       E_Photz       ?="" Upper Uncertainty in Photz
##  57- 60 F4.2   ---       e_Photz       ?="" Lower Uncertainty in Photz
##  62- 63 I2     ---         Richness    ?="" Richness for Primary Photometric Redshift (3)
##  65- 65 I1     ---       e_Richness    ?="" Uncertainty in Richness for Primary Photometric Redshift
##  67- 69 F3.1   10+14Msun   M500        ?="" M500, enclosed mass
##  71- 73 F3.1   10+14Msun E_M500        ?="" Upper Uncertainty in M500
##  75- 77 F3.1   10+14Msun e_M500        ?="" Lower Uncertainty in M500
##  79- 79 A1     ---         Comments    Literature Names and Comments (4)
##  81- 83 A3     ---         References  References (5)


#######################################################################################################################
## apjsaafad2t5_mrt.txt
#######################################################################################################################
##   1- 14 A14    ---         Cluster     Cluster Name (1)
##  16- 17 I2     h           RAh         Hour of Right Ascension (J2000) 
##  19- 20 I2     min         RAm         Minute of Right Ascension (J2000) 
##  22- 25 F4.1   s           RAs         Second of Right Ascension (J2000) 
##  27- 27 A1     ---         DE-         Sign of the Declination (J2000)
##  28- 29 I2     deg         DEd         Degree of Declination (J2000) 
##  31- 32 I2     arcmin      DEm         Arcminute of Declination (J2000) 
##  34- 35 I2     arcsec      DEs         Arcsecond of Declination (J2000) 
##  37- 40 F4.2   ---         Photz       ?="" Primary Photometric Redshift
##  42- 45 F4.2   ---       E_Photz       ?="" Upper Uncertainty in Photz
##  47- 50 F4.2   ---       e_Photz       ?="" Lower Uncertainty in Photz
##  52- 55 F4.2   ---         Photz2      ?="" Secondary Photometric Redshift
##  57- 60 F4.2   ---       E_Photz2      ?="" Upper Uncertainty in Photz2
##  62- 65 F4.2   ---       e_Photz2      ?="" Lower Uncertainty in Photz2
##  67- 71 F5.3   ---         Specz       ?="" Spectroscopic Redshift
##  73- 74 I2     ---       o_Specz       ?="" Number of Confirmed Cluster Members
##  76- 77 I2     ---         Richness    ?="" Richness for Primary Photometric Redshift
##  79- 79 I1     ---       e_Richness    ?="" Uncertainty in Richness for Primary Photometric Redshift
##  81- 82 I2     ---         Richness2   ?="" Richness for Primary Secondary Redshift
##  84- 84 I1     ---       e_Richness2   ?="" Uncertainty in Richness for Primary Photometric Redshift
##  86- 88 F3.1   10+14Msun   M500        ?="" M500, enclosed mass
##  90- 92 F3.1   10+14Msun E_M500        ?="" Upper Uncertainty in M500
##  94- 96 F3.1   10+14Msun e_M500        ?="" Lower Uncertainty in M500
##  98-120 A23    ---         Comments    Literature Names and Comments (2)
## 122-126 A5     ---         References  References (3)

totphot1 = [np.array([]) for i in range(3)]
totphot2 = [np.array([]) for i in range(3)]
totspec1 = [np.array([]) for i in range(3)]
totrich1 = [np.array([]) for i in range(3)]
totrich2 = [np.array([]) for i in range(3)]

with open('madcows/apjsaafad2t3_mrt.txt','r') as file:
  for line in range(71): file.readline()
  for line in file.read().splitlines():
    totphot1[0] = np.append(totphot1[0],float(line[46:50].replace(' ','0')))
    totphot1[1] = np.append(totphot1[1],float(line[56:60].replace(' ','0')))
    totphot1[2] = np.append(totphot1[2],float(line[51:55].replace(' ','0')))

    totphot2[0] = np.append(totphot2[0],float(line[61:65].replace(' ','0')))
    totphot2[1] = np.append(totphot2[1],float(line[71:75].replace(' ','0')))
    totphot2[2] = np.append(totphot2[2],float(line[66:70].replace(' ','0')))

    totspec1[0] = np.append(totspec1[0],float(line[76:82].replace(' ','0')))

    totrich1[0] = np.append(totrich1[0],float(line[86:88].replace(' ','0')))
    totrich1[1] = np.append(totrich1[1],float(line[89:90].replace(' ','0')))
    totrich1[2] = np.append(totrich1[2],float(line[89:90].replace(' ','0')))

    totrich2[0] = np.append(totrich2[0],float(line[91:93].replace(' ','0')))
    totrich2[1] = np.append(totrich2[1],float(line[94:95].replace(' ','0')))
    totrich2[2] = np.append(totrich2[2],float(line[94:95].replace(' ','0')))

with open('madcows/apjsaafad2t4_mrt.txt','r') as file:
  for line in range(49): file.readline()
  for line in file.read().splitlines():
    totphot1[0] = np.append(totphot1[0],float(line[46:50].replace(' ','0')))
    totphot1[1] = np.append(totphot1[1],float(line[56:60].replace(' ','0')))
    totphot1[2] = np.append(totphot1[2],float(line[51:55].replace(' ','0')))

    totrich1[0] = np.append(totrich1[0],float(line[61:63].replace(' ','0')))
    totrich1[1] = np.append(totrich1[1],float(line[64:65].replace(' ','0')))
    totrich1[2] = np.append(totrich1[2],float(line[64:65].replace(' ','0')))

with open('madcows/apjsaafad2t5_mrt.txt','r') as file:
  for line in range(58): file.readline()
  for line in file.read().splitlines():
    totphot1[0] = np.append(totphot1[0],float(line[36:40].replace(' ','0')))
    totphot1[1] = np.append(totphot1[1],float(line[46:50].replace(' ','0')))
    totphot1[2] = np.append(totphot1[2],float(line[41:45].replace(' ','0')))

    totphot2[0] = np.append(totphot2[0],float(line[51:55].replace(' ','0')))
    totphot2[1] = np.append(totphot2[1],float(line[61:65].replace(' ','0')))
    totphot2[2] = np.append(totphot2[2],float(line[56:60].replace(' ','0')))

    totspec1[0] = np.append(totspec1[0],float(line[66:71].replace(' ','0')))

    totrich1[0] = np.append(totrich1[0],float(line[75:77].replace(' ','0')))
    totrich1[1] = np.append(totrich1[1],float(line[78:79].replace(' ','0')))
    totrich1[2] = np.append(totrich1[2],float(line[78:79].replace(' ','0')))

    totrich2[0] = np.append(totrich2[0],float(line[80:82].replace(' ','0')))
    totrich2[1] = np.append(totrich2[1],float(line[83:84].replace(' ','0')))
    totrich2[2] = np.append(totrich2[2],float(line[83:84].replace(' ','0')))

######################################################################

sptlis = ['samples/2500d_cluster_sample_Bocquet19.fits',
          'samples/sptpol100d_catalog_huang19.fits',
          'samples/sptecs_catalog_oct919.fits']
spthdu  = [fits.open(spt)[1] for s, spt in enumerate(sptlis)]
sptm500 = [np.array([spthdu[s].data['M500'],spthdu[s].data['M500_lerr'],spthdu[s].data['M500_uerr']]) for s, spt in enumerate(sptlis)]
sptzeta = [np.array([spthdu[s].data['REDSHIFT'],spthdu[s].data['REDSHIFT_UNC'],spthdu[s].data['REDSHIFT_UNC']]) for s, spt in enumerate(sptlis)]

sptm500 = np.array([np.concatenate((sptm500[0][0],sptm500[1][0],sptm500[2][0]))])
sptzeta = np.array([np.concatenate((sptzeta[0][0],sptzeta[1][0],sptzeta[2][0]))])

######################################################################

plclis = ['samples/HFI_PCCS_SZ-union_R2.08.fits']
plchdu  = [fits.open(plc)[1] for s, plc in enumerate(plclis)]
plcm500 = [np.array([plchdu[p].data['MSZ'],plchdu[p].data['MSZ_ERR_LOW'],plchdu[p].data['MSZ_ERR_UP']]) for p, plc in enumerate(plclis)]
plczeta = [np.array([plchdu[p].data['REDSHIFT'],np.zeros(plchdu[p].data['REDSHIFT'].shape[0]),np.zeros(plchdu[p].data['REDSHIFT'].shape[0])]) for p, plc in enumerate(plclis)]

plcm500 = plcm500[0]
plczeta = plczeta[0]

######################################################################

actlis = ['samples/E-D56Clusters.fits']
acthdu  = [fits.open(act)[1] for s, act in enumerate(actlis)]
actm500 = [np.array([acthdu[a].data['M500cCal'],acthdu[a].data['M500cCal_errMinus'],acthdu[a].data['M500cCal_errPlus']]) for a, act in enumerate(actlis)]
actzeta = [np.array([acthdu[a].data['z'],acthdu[a].data['zErr'],acthdu[a].data['zErr']]) for a, act in enumerate(actlis)]

actm500.append(np.array([[10.50, 3.50, 3.50, 3.40, 2.50, 4.10, 3.10, 6.70, 5.70, 5.40, 4.40, 7.00,
                           4.00, 4.00, 2.30, 3.90, 4.30, 4.60, 7.50, 1.40, 5.10,10.30, 1.70]]))
actm500.append(np.array([[3.900,3.000,5.700,2.900,3.100,5.500,4.400,2.700,5.600,2.200,3.200,5.200,2.600,3.300,3.300,2.100,5.700,
                          3.100,4.300,3.500,3.800,3.000,2.800,1.400,3.800,2.400,2.800,6.700,3.300,2.500,2.500,2.700,3.800,2.200,
                          4.700,4.200,2.700,3.300,5.500,2.500,3.700,3.100,3.100,3.000,4.600,3.800,5.100,2.500,5.300,3.500,5.500,
                          5.700,5.300,3.200,2.800,6.300,2.600,3.000,4.300,4.100,2.500,2.700,2.700,3.700,4.200,9.400,6.100,3.200]]))
actzeta.append(np.array([[0.870,0.118,0.480,0.343,0.556,0.278,0.334,0.300,0.392,0.442,0.530,0.421,
                          0.461,0.294,0.768,1.066,0.609,0.684,0.222,0.146,0.167,0.296,0.296]]))
actzeta.append(np.array([[0.360,1.360,0.533,0.211,0.750,0.805,0.650,1.110,0.545,0.690,0.760,0.786,0.277,0.720,0.379,0.700,0.230,
                          0.450,0.676,0.865,0.672,0.537,0.350,0.589,0.663,0.720,0.440,0.375,0.620,0.684,0.179,0.780,0.363,0.530,
                          0.167,0.153,0.633,0.384,0.448,1.320,1.070,0.310,0.297,0.345,0.340,0.622,0.333,0.333,0.321,0.408,0.320,
                          0.385,0.234,0.710,0.330,0.231,0.118,0.690,0.488,0.224,0.570,0.610,0.540,0.520,0.360,0.705,0.275,0.990]]))

actm500 = np.array([np.concatenate((actm500[0][0],actm500[1][0],actm500[2][0]))])
actzeta = np.array([np.concatenate((actzeta[0][0],actzeta[1][0],actzeta[2][0]))])