library(readxl)
Hausser <- readxl::read_xlsx("~/CASCADEpaper/paper/Fig5_archetype/41467_2019_13195_MOESM4_ESM.xlsx", sheet = 2)
head(Hausser)
features <- split(Hausser$`Feature Name`, Hausser$`archetype #`)
head(features)

lapply(clean_module, function(x){
  sapply(features, function(y) intersect(x, y ))
})

lapply(clean_module, function(x){
  sapply(features, function(y) length(intersect(x, y )))
})

$AR
1  2  3  4  5 
6  1  2  5 14 

$inflammation
1   2   3   4   5 
43  11  41  66 103 

$NE1
1  2  3  4  5 
12  3 10 12 48 

$NE2
1  2  3  4  5 
18  5 21 16 53 

$cycling
1  2  3  4  5 
10  6 43 56 12 

$hypoxia
1 2 3 4 5 
4 1 7 7 3 

$AR
$AR$`1`
[1] "ACACA"  "AFF3"   "AR"     "BMPR1B" "GALNT7"
[6] "UNC13B"

$AR$`2`
[1] "KLK3"

$AR$`3`
[1] "FOLH1"  "MYBPC1"

$AR$`4`
[1] "BMPR1B"  "COLEC12" "GPC6"    "PDE4D"  
[5] "PMEPA1" 

$AR$`5`
[1] "AFF3"     "AR"       "BMPR1B"   "COLEC12" 
[5] "GPC6"     "IL1RAPL1" "KLK2"     "KLK3"    
[9] "NFIA"     "PDE4D"    "PMEPA1"   "SLC4A4"  
[13] "TMEFF2"   "ZBTB16"  


$inflammation
$inflammation$`1`
[1] "AFF1"      "ALCAM"     "ARHGEF10L"
[4] "ARHGEF12"  "ATG2B"     "BAZ2B"    
[7] "C5"        "CPEB4"     "CSNK1A1"  
[10] "DNAH5"     "ELF2"      "ERBB4"    
[13] "ESRRG"     "FAAH2"     "FAM160A1" 
[16] "FNBP1L"    "FSTL4"     "GRHL2"    
[19] "LIMCH1"    "LRRIQ1"    "MAST4"    
[22] "NAALADL2"  "NBEA"      "NEBL"     
[25] "NEK10"     "PCCA"      "PDZD8"    
[28] "PPM1B"     "RALGAPA2"  "RASEF"    
[31] "RBM47"     "RPS6KA5"   "SIDT1"    
[34] "SLC11A2"   "SLC23A2"   "STRN3"    
[37] "TAOK1"     "TBC1D8"    "TC2N"     
[40] "VAV3"      "VPS13C"    "XBP1"     
[43] "ZMYND8"   

$inflammation$`2`
[1] "ASRGL1"   "BACE2"    "BAIAP2L1" "FAAH2"   
[5] "GMDS"     "IDH2"     "PITPNC1"  "SERGEF"  
[9] "SND1"     "SOD2"     "VSTM2L"  

$inflammation$`3`
[1] "ABCC4"     "ADAMTS3"   "ARHGEF10L"
[4] "ASRGL1"    "B4GALT5"   "BACE2"    
[7] "BAIAP2L1"  "CDC14B"    "CDK8"     
[10] "CP"        "EGFR"      "ETS2"     
[13] "ETV6"      "FAM135A"   "FAM160A1" 
[16] "FNDC3B"    "GMDS"      "HIVEP2"   
[19] "IDH2"      "IFNGR1"    "LTBP1"    
[22] "MAPK4"     "MGST1"     "MYO10"    
[25] "NAMPT"     "NCOA7"     "NPAS2"    
[28] "NT5C2"     "PLSCR1"    "PTPRK"    
[31] "RNF24"     "ROR1"      "SLC23A2"  
[34] "SND1"      "SOD2"      "SPINK1"   
[37] "STK3"      "TCF7L2"    "TNFRSF21" 
[40] "TTC39B"    "YAP1"     

$inflammation$`4`
[1] "ABCC4"    "ADAMTS3"  "AFF1"     "ALPK1"   
[5] "ARHGAP18" "ARHGAP26" "ARHGEF3"  "ARID5B"  
[9] "ASCC3"    "B4GALT5"  "BACE2"    "CASP10"  
[13] "CP"       "CSNK1A1"  "DIAPH1"   "DNAH5"   
[17] "DOCK4"    "EGFR"     "ETS2"     "ETV6"    
[21] "FAM135A"  "FAM160A1" "FKBP5"    "FNDC3B"  
[25] "FNIP2"    "GLIS3"    "HIF1A"    "IFNGR1"  
[29] "ITPR2"    "KITLG"    "KMO"      "LPP"     
[33] "LTBP1"    "MACF1"    "MAN1A1"   "MYO9A"   
[37] "NAMPT"    "NFKB1"    "NRP1"     "OSMR"    
[41] "PLSCR1"   "PPP2R5E"  "PRR16"    "PTPRK"   
[45] "RNF144B"  "RNF24"    "ROR1"     "RUNX1"   
[49] "SFMBT2"   "SGPP1"    "SLC4A7"   "SOD2"    
[53] "ST6GAL1"  "STAT3"    "STK3"     "SYNE2"   
[57] "TANC1"    "TAOK1"    "TNFAIP8"  "TRIO"    
[61] "TRPS1"    "TTC39B"   "UTRN"     "VPS13C"  
[65] "XPR1"     "YAP1"    

$inflammation$`5`
[1] "ABCC4"      "ADAMTS3"    "AFF1"      
[4] "ALPK1"      "ARHGEF10L"  "ARHGEF12"  
[7] "ARHGEF3"    "ARHGEF38"   "ARID5B"    
[10] "ASRGL1"     "ATP8A2"     "B4GALT5"   
[13] "BACE2"      "BAZ2B"      "BCL6"      
[16] "C5"         "CACNA1C"    "CASP10"    
[19] "CATSPERB"   "CDC14B"     "CP"        
[22] "CPEB4"      "CRIM1"      "CRTC3"     
[25] "CSGALNACT1" "CYTH3"      "DOCK4"     
[28] "EBF1"       "EGFR"       "ELF2"      
[31] "EPAS1"      "ERBB4"      "ETV6"      
[34] "FAM135A"    "FKBP5"      "FNDC3B"    
[37] "FNIP2"      "GLIS3"      "HIF1A"     
[40] "HIVEP2"     "IFNGR1"     "KCNN2"     
[43] "KMO"        "LIMCH1"     "LPP"       
[46] "LRIG1"      "LTBP1"      "MACF1"     
[49] "MAPK4"      "MAST4"      "MECOM"     
[52] "MTUS1"      "MYO9A"      "NAALADL2"  
[55] "NBEA"       "NFKB1"      "NRP1"      
[58] "NTN4"       "OSMR"       "PCCA"      
[61] "PELI2"      "PIP5K1B"    "PLXNA4"    
[64] "PRKD1"      "PRR16"      "RALGAPA2"  
[67] "RBMS3"      "RNF144B"    "RNF24"     
[70] "ROR1"       "RORA"       "RPS6KA5"   
[73] "RUFY3"      "RUNX1"      "SFMBT2"    
[76] "SGPP1"      "SHISA6"     "SIDT1"     
[79] "SIPA1L2"    "SLC1A1"     "SLC23A2"   
[82] "SLC4A7"     "SLC7A2"     "SLC9A9"    
[85] "SOD2"       "ST6GAL1"    "STAT3"     
[88] "STEAP4"     "STRN3"      "SYNE2"     
[91] "TANC1"      "TBC1D8"     "TC2N"      
[94] "TCF7L2"     "TMC5"       "TNFAIP8"   
[97] "TRIO"       "TRPS1"      "UTRN"      
[100] "VPS13C"     "YAP1"       "ZBTB20"    
[103] "ZHX2"      


$NE1
$NE1$`1`
[1] "CACNA1D" "CDH7"    "CDHR3"   "IQSEC1" 
[5] "KSR2"    "NELL1"   "PBX1"    "PRKCE"  
[9] "SAMD12"  "STAU2"   "TRAPPC9" "ZNF704" 

$NE1$`2`
[1] "IGFBP2"  "SLC18A1" "TRAPPC9"

$NE1$`3`
[1] "ADCY2"   "CADM1"   "DNER"    "ERC2"   
[5] "GABRG3"  "NOL4"    "PID1"    "PLXNA2" 
[9] "ST18"    "TMEM108"

$NE1$`4`
[1] "CDK14"  "CHST11" "DNER"   "ERC2"   "FMN1"  
[6] "JAZF1"  "NFATC2" "PBX3"   "RASSF8" "SORCS3"
[11] "TP63"   "TSHZ2" 

$NE1$`5`
[1] "ADCY2"    "ANK2"     "BTBD11"   "CACNA2D3"
[5] "CACNB2"   "CADPS"    "CDH13"    "CDK14"   
[9] "CHST11"   "CPE"      "CSMD2"    "DNER"    
[13] "ETV1"     "FGF14"    "GLCCI1"   "GPC3"    
[17] "IQSEC1"   "JAZF1"    "LSAMP"    "MAGI2"   
[21] "MGAT4C"   "MMP16"    "NFATC2"   "NHS"     
[25] "PBX1"     "PBX3"     "PDE5A"    "PID1"    
[29] "PRKCE"    "RASSF8"   "RGS7"     "RMST"    
[33] "RUNX1T1"  "SDK1"     "SETBP1"   "SLC44A5" 
[37] "SLIT2"    "SPOCK1"   "SSBP2"    "ST18"    
[41] "STIM2"    "TMEM108"  "TMEM150C" "TP63"    
[45] "TSHZ2"    "XKR6"     "ZBTB7C"   "ZNF704"  


$NE2
$NE2$`1`
[1] "ANKS1B"   "ASTN2"    "CBFA2T2"  "DACH1"   
[5] "FHIT"     "GDAP1"    "HEPACAM2" "HS3ST5"  
[9] "IMMP2L"   "MAGI1"    "MAML3"    "MLLT3"   
[13] "NRXN3"    "SIPA1L3"  "SRGAP1"   "SYBU"    
[17] "SYTL2"    "VPS13A"  

$NE2$`2`
[1] "CMTM7"  "FHIT"   "IMMP2L" "PVT1"   "SAAL1" 

$NE2$`3`
[1] "CCDC50"  "CHD7"    "CMTM7"   "FIGN"   
[5] "HDAC9"   "KALRN"   "KCNB2"   "NDST3"  
[9] "NFASC"   "NFIB"    "NPAS3"   "PLCB4"  
[13] "PRIM2"   "RANBP17" "RFX3"    "SAAL1"  
[17] "SDK2"    "SH3GL2"  "SLAIN1"  "SLCO5A1"
[21] "WASF3"  

$NE2$`4`
[1] "ACPP"   "CCDC50" "CMTM7"  "CRYBG3" "FIGN"  
[6] "HDAC9"  "LPIN2"  "MARCH1" "MCC"    "MIER1" 
[11] "PLCL2"  "PRIM2"  "PRUNE2" "RFX3"   "ROBO1" 
[16] "VPS13A"

$NE2$`5`
[1] "ACPP"     "AMPH"     "ARHGEF7"  "ARL15"   
[5] "ASTN2"    "ATP8A1"   "CACNA2D1" "CADM2"   
[9] "CBLB"     "CCDC50"   "CMTM7"    "CRYBG3"  
[13] "DACH1"    "ELMO1"    "FAT4"     "FBXL7"   
[17] "FIGN"     "HDAC9"    "HEPACAM2" "JAKMIP2" 
[21] "KALRN"    "KCNMA1"   "KIAA1211" "LHFPL3"  
[25] "LPIN2"    "MAGI1"    "MAML3"    "MARCH1"  
[29] "MCC"      "MIER1"    "MLLT3"    "MYRIP"   
[33] "NFASC"    "NFIB"     "NPAS3"    "NRXN1"   
[37] "NRXN3"    "PCDH15"   "PLCB4"    "PLCL2"   
[41] "PLCXD3"   "PRUNE2"   "ROBO1"    "SDK2"    
[45] "SLAIN1"   "SLC38A11" "SLCO5A1"  "SRGAP1"  
[49] "STXBP5L"  "SYBU"     "TNIK"     "TSHR"    
[53] "WASF3"   


$cycling
$cycling$`1`
[1] "CEP192"  "G2E3"    "HP1BP3"  "LIN52"  
[5] "MSI2"    "NEDD4L"  "PHACTR4" "PKP4"   
[9] "TTC39C"  "ZRANB3" 

$cycling$`2`
[1] "ANKRD37" "CDC25C"  "CENPP"   "RGS3"   
[5] "TACC3"   "VRK1"   

$cycling$`3`
[1] "ANLN"      "ARHGAP11B" "ASPM"     
[4] "ATAD2"     "BRIP1"     "BUB1"     
[7] "CCDC18"    "CENPI"     "CEP192"   
[10] "CIT"       "CKAP5"     "CREB3L2"  
[13] "DIAPH3"    "DNA2"      "E2F3"     
[16] "EZH2"      "GPSM2"     "KIF11"    
[19] "KIF14"     "KIF15"     "KIF18B"   
[22] "KIF4A"     "MELK"      "MID1"     
[25] "NCAPG2"    "NDC80"     "NT5DC3"   
[28] "NUCKS1"    "NUF2"      "NUSAP1"   
[31] "OSBPL3"    "PITPNM2"   "POLQ"     
[34] "PPM1E"     "PRR11"     "R3HDM1"   
[37] "RFC3"      "SMC4"      "STIL"     
[40] "TACC3"     "TPX2"      "VRK1"     
[43] "XPO1"     

$cycling$`4`
[1] "ANLN"      "ARHGAP11B" "ASPM"     
[4] "ATAD2"     "BRCA1"     "BRIP1"    
[7] "BUB1"      "CCDC18"    "CDC25C"   
[10] "CENPI"     "CENPK"     "CIT"      
[13] "CKAP5"     "DEPDC1B"   "DIAPH3"   
[16] "DLEU2"     "DNA2"      "ECT2"     
[19] "EZH2"      "GINS1"     "GPSM2"    
[22] "KIF11"     "KIF14"     "KIF15"    
[25] "KIF18B"    "KIF4A"     "KNTC1"    
[28] "LIN52"     "LMNB1"     "MELK"     
[31] "MPHOSPH9"  "NCAPG2"    "NDC80"    
[34] "NUSAP1"    "OSBPL3"    "PHACTR4"  
[37] "PITPNM2"   "POLQ"      "PRR11"    
[40] "R3HDM1"    "RACGAP1"   "RFC3"     
[43] "RRM2"      "SESN3"     "SKA2"     
[46] "SMC4"      "SPC25"     "STAG1"    
[49] "STIL"      "TACC3"     "TOP2A"    
[52] "TPX2"      "TTC39C"    "VRK1"     
[55] "XPO1"      "ZRANB3"   

$cycling$`5`
[1] "CREB3L2" "DLG2"    "HERC3"   "HP1BP3" 
[5] "LIFR"    "LIN52"   "MID1"    "OSBPL3" 
[9] "PITPNM2" "RGS3"    "SESN3"   "WSB1"   


$hypoxia
$hypoxia$`1`
[1] "ANKRD12" "KDM3A"   "PIAS2"   "PRELID2"

$hypoxia$`2`
[1] "ENO1"

$hypoxia$`3`
[1] "ABTB2"  "AMD1"   "ATP2C1" "ENO1"   "PIAS2" 
[6] "PLOD2"  "RAB1A" 

$hypoxia$`4`
[1] "AMD1"   "ATP2C1" "ENO1"   "GBE1"   "P4HA1" 
[6] "PGK1"   "PLOD2" 

$hypoxia$`5`
[1] "ANKRD12" "DTNA"    "GBE1"   

