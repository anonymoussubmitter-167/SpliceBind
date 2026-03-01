
        from pymol import cmd,stored
        
        set depth_cue, 1
        set fog_start, 0.4
        
        set_color b_col, [36,36,85]
        set_color t_col, [10,10,10]
        set bg_rgb_bottom, b_col
        set bg_rgb_top, t_col      
        set bg_gradient
        
        set  spec_power  =  200
        set  spec_refl   =  0
        
        load "data/AF-P23458-F1-model_v6.pdb", protein
        create ligands, protein and organic
        select xlig, protein and organic
        delete xlig
        
        hide everything, all
        
        color white, elem c
        color bluewhite, protein
        #show_as cartoon, protein
        show surface, protein
        #set transparency, 0.15
        
        show sticks, ligands
        set stick_color, magenta
        
        
        
        
        # SAS points
 
        load "data/AF-P23458-F1-model_v6.pdb_points.pdb.gz", points
        hide nonbonded, points
        show nb_spheres, points
        set sphere_scale, 0.2, points
        cmd.spectrum("b", "green_red", selection="points", minimum=0, maximum=0.7)
        
        
        stored.list=[]
        cmd.iterate("(resn STP)","stored.list.append(resi)")    # read info about residues STP
        lastSTP=stored.list[-1] # get the index of the last residue
        hide lines, resn STP
        
        cmd.select("rest", "resn STP and resi 0")
        
        for my_index in range(1,int(lastSTP)+1): cmd.select("pocket"+str(my_index), "resn STP and resi "+str(my_index))
        for my_index in range(1,int(lastSTP)+1): cmd.show("spheres","pocket"+str(my_index))
        for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_scale","0.4","pocket"+str(my_index))
        for my_index in range(1,int(lastSTP)+1): cmd.set("sphere_transparency","0.1","pocket"+str(my_index))
        
        
        
        set_color pcol1 = [0.361,0.576,0.902]
select surf_pocket1, protein and id [457,3524,3526,4147,3529,3816,3819,3820,3530,4101,488,4122,3368,3395,901,3366,3398,3399,3402,3512,3514,3536,3538,456,719,721,454,449,459,461,463,3339,3341,723,724,705,3340,3432,3434,3435,4521,3562,3610,3814,3810,3619,558,3856,3852,3854,589,591,592,587,590,3618,559,3408,3410,616,606,480,614,618,585,582,462,486,487,494,3865,900,4522,3363,3523,3549,3551,3553,3554] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.353,0.702]
select surf_pocket2, protein and id [5841,4836,7687,6108,6109,6110,6111,6116,6120,7685,5840,4827,4828,4829,4830,4854,4872,6363,6148,6150,6153,6154,6418,6500,6501,5836,5503,5839,7095,7091,5825,7072,7070,7686,7055,7670,6456,7088,7079,7668,7680,7681,6417,6419,7075,6431,6436,6439,5848,5849,6003,6005,5806,5811,5812,6122,6125,6121,5999,5071,5072,5073,5997,6000,6008,6009,5807,6144,5992,6012] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.388,0.361,0.902]
select surf_pocket3, protein and id [5461,5476,5478,5487,5865,5486,5837,7092,4872,5055,4871,4806,4809,5840,4827,4815,4819,5054,5056,5455,5456,5468,7301,7123,7299,5510,5507,4807,4812,7122,7089,7090,4810,5864,5992,5436,5999,5072,5073,5997,6008,6009,5849,6003,6005,5551] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.396,0.278,0.702]
select surf_pocket4, protein and id [2407,5652,5673,5646,5648,5650,5674,6873,6899,6907,6905,5680,2408,2409,5682,5670,5676,5668,5649,5675,2398,2402,2403,3164,3211,3214,3215,3148,3216,3080,3081,3160,3093,3096,3146,3163,3095,3169,3194,3198,3192,3193,3195,3071,5702,6838,2395] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.631,0.361,0.902]
select surf_pocket5, protein and id [2313,2315,697,884,905,876,699,2281,2307,2300,2301,2302,673,676,634,635,637,1386,2294,1384,1351,2296,871,1349,868,1317,1401,2272,1385,1389,1388,2299,907,909,2298,2316,674] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.584,0.278,0.702]
select surf_pocket6, protein and id [7168,7169,7166,7824,7823,7825,7747,7748,8271,8273,8275,7769,7357,7753,7765,7767,7164,7165,7153,7224,7226,7227,7171,7172,7174,8189,8160,8168,7794,8169,8265,8266,8268,8171,8174,8188,8187,7773,7614,7615,7791] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.875,0.361,0.902]
select surf_pocket7, protein and id [3570,6696,6700,6762,3568,3453,6790,6792,6794,6779,6781,6785,6821,6795,6778,3464,3467,3544,3577,6849,6855,6856,6857,6754,6760,4449,6793,6771,3440,3565,3569,4453,4455,4454,4456,3571,3540,3543] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.627]
select surf_pocket8, protein and id [7380,7385,7378,7679,7234,7236,7218,7235,4821,7083,7723,7724,7084,7678,7684,4863,7127,7130,4822,4825,7111,7110,7113,7364,7250,5109,4864,5096,4842,4844,5118,5119,4829,4837,4838,4835,6090,7683,6088,6089] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.361,0.682]
select surf_pocket9, protein and id [1326,1332,1333,1300,1302,1607,1308,1309,1331,318,1272,334,337,339,327,333,346,1366,1772,1785,1787,1788,1790,1771,1738,1740,1618,1620,1774,1732,1606,1610,1612,1616,1637,1639,1640,1335] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.278,0.439]
select surf_pocket10, protein and id [4448,6765,4298,4299,3576,6766,4242,4268,4432,4433,4434,3780,3582,3583,3585,3589,4240,4262,4264,4267] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.361,0.439]
select surf_pocket11, protein and id [8455,8470,8506,8508,8534,8585,8586,8404,8407,8566,8129,8915,8918,8543,8917,8540,8437,8457,8406,8421,8423,8297,8298] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.702,0.310,0.278]
select surf_pocket12, protein and id [569,1207,1208,4510,4512,3830,3631,3633,4513] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.522,0.361]
select surf_pocket13, protein and id [6254,6256,6257,6250,6251,6248,5749,4471,4472,6791,6793,6283,6284,6771,6747,6258,6261,4462,4463,4464] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.502,0.278]
select surf_pocket14, protein and id [8166,8495,8496,7815,7816,8707,8745,8794,8797,8798,8781,8783,8161,8162,7813,7814,8747] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.765,0.361]
select surf_pocket15, protein and id [1850,1852,1516,1506,1434,1507,1511,1481,1414,1435,1821,1823,1825,1828,1407,1408,1539,1405] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.702,0.690,0.278]
select surf_pocket16, protein and id [9168,9167,9012,8822] 
set surface_color,  pcol16, surf_pocket16 
   
        
        deselect
        
        orient
        