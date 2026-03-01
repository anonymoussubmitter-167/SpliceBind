
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
        
        load "data/AF-P29597-F1-model_v6.pdb", protein
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
 
        load "data/AF-P29597-F1-model_v6.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [5850,5851,5490,5820,5821,5822,5823,6380,6329,5530,6367,5473,5474,5054,4763,5053,5827,5438,5464,5465,5437,5442,5449,5812,5798,5799,5975,5835,5988,6096,6098,6132,6088,6086,5992,6099,6100,5070,5071,4761,4716,4722,4764,4749,4731,4723,4725,4728,4710,4706,5435,4703,4704,4707,4702,4709,6131,6121,6122,6128] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.404,0.702]
select surf_pocket2, protein and id [1213,1244,1245,1236,1237,794,798,571,585,774,776,608,811,2181,2182,2184,2188,2196,2198,629,803,627,572,573,605,609,1184,836,1241,1243,1239,1268,1269,2162,2180,1246,1284,1287,2178,2192,2194] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.361,0.463,0.902]
select surf_pocket3, protein and id [7718,8144,8145,8146,7723,5253,5255,5716,5718,5254,5720,5751,5753,5754,8153,7713,8196,8197,3163,5257,7597,7598,7599,3025,3026,3159,5741,7540,7578,8180,3027,3150,3167,3165,3154,3160,2969,3023] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.278,0.314,0.702]
select surf_pocket4, protein and id [7075,7076,7077,7231,7262,7263,7265,7072,7074,5029,5034,5422,5426,5430,5439,5322,5323,5418,7094,7276,7209,7211,5961,5423,5424,7293,5025,5030] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.373,0.361,0.902]
select surf_pocket5, protein and id [7757,8130,7760,7762,7789,7753,7754,7742,7117,7120,7178,7304,7735,7121,7728,7785,7116,7122,7118,7123,7176,7135,7179,7319,7709,7320,7321,7136,8115,8213,8101,7128,7130,8129,7569,7708,8210,7453] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.333,0.278,0.702]
select surf_pocket6, protein and id [3528,3494,3495,3499,495,3579,3588,3587,3758,3760,3761,3762,3764,524,526,527,528,529,523,515,503,496,500,4433,3478,3480,392,519,521,3802,3400,3501,3503,3506,3526,3363,3364,3367,3399,3407,4432,4429,4431,3578,3564] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.486,0.361,0.902]
select surf_pocket7, protein and id [2286,2300,2307,2308,2309,2310,2297,6988,6990,6989,6991,7005,7220,6976,7282,3164,7280,3171,5269,5271,5272,3183,3184,7283,7285,7598,7018,7284,7595,3181,7274,7275,7270,7272,7216,7218,7223,7235,7237,7238,7277] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.424,0.278,0.702]
select surf_pocket8, protein and id [5782,5789,6217,6216,5771,5218,5775,4544,4547,6014,6015,5160,5191,5192,5193,6018,5188,5187,5189,5989,6003,6004,6006,5200,5987] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.600,0.361,0.902]
select surf_pocket9, protein and id [1458,1935,1937,1936,1930,1932,1933,1928,1931,1459,1460,1480,1490,1491,1492,1482,1484,2013,2014,1972,2075,2076] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.514,0.278,0.702]
select surf_pocket10, protein and id [8393,8391,8406,8440,8442,8504,8506,8509,8476,8479,8865,8363,8379,8841,8843,8866,8345,8528,8529,8470] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.718,0.361,0.902]
select surf_pocket11, protein and id [4352,4179,4205,4206,4207,4208,4204,3557,3552,4354,4349,4353,4367,4368,4231,4232,4233,3541,4180,3556,3550,4172,3744,6724] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.600,0.278,0.702]
select surf_pocket12, protein and id [8362,8365,8367,7382,8363,7384,7376,7377,7378,7383,7147,7148,8240,7411,7159,7157,7390,7391,8372,8376,8377,8253,8232,8254,8366,7445] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.831,0.361,0.902]
select surf_pocket13, protein and id [2965,2977,3150,3182,3167,3179,3181,2980,2270,2273,2274,2584,2982,2279,7576,7521,7596,3184,7598,7595,6990,6992,6991,7583,6999] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.690,0.278,0.702]
select surf_pocket14, protein and id [970,972,1031,1094,990,991,993,1002,998,1016,1028,989,987,4408,9276,9277,9309,4494,4495,4497,4410,4411,4482,4457] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.361,0.855]
select surf_pocket15, protein and id [3536,3533,3535,3537,3538,3406,4373,6805,3429,6740,6719,6721,6730,6737,6738,6714,6749,6744,3415,4390,4391,6777,6778,6658,6660] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.702,0.278,0.620]
select surf_pocket16, protein and id [2244,3219,3221,3205,3206,3207,2247,2266,2332,3106,3108,2337,2341,2345,2348,3101,3102] 
set surface_color,  pcol16, surf_pocket16 
set_color pcol17 = [0.902,0.361,0.741]
select surf_pocket17, protein and id [3775,3593,3753,3780,3781,3841,3997,3998,3999,3774,3787,3843,3638] 
set surface_color,  pcol17, surf_pocket17 
set_color pcol18 = [0.702,0.278,0.533]
select surf_pocket18, protein and id [2292,5669,2300,2297,5667,5639,3173,3174,5282,5655,5663,3171,5271,2296,2291,2293,2294,2295,3175,5687,6794,5657] 
set surface_color,  pcol18, surf_pocket18 
set_color pcol19 = [0.902,0.361,0.627]
select surf_pocket19, protein and id [9198,9200,8003,8545,8579,8542,8544,8554,8555,8557,7999,8549,8551,8547,9196,9202] 
set surface_color,  pcol19, surf_pocket19 
set_color pcol20 = [0.702,0.278,0.443]
select surf_pocket20, protein and id [638,645,3275,668,669,3254,3255,3257,616,599,747,670,728,729,581,748] 
set surface_color,  pcol20, surf_pocket20 
set_color pcol21 = [0.902,0.361,0.514]
select surf_pocket21, protein and id [7983,7531,7534,7954,7979,7951,7953,7955,7535,7538,7543,4462,7956,7958,4470,947,948] 
set surface_color,  pcol21, surf_pocket21 
set_color pcol22 = [0.702,0.278,0.353]
select surf_pocket22, protein and id [553,647,649,651,652,827,828,543,386,3478,3480,390,392,524,527,516,415,413,3304,3305,3367] 
set surface_color,  pcol22, surf_pocket22 
set_color pcol23 = [0.902,0.361,0.400]
select surf_pocket23, protein and id [605,611,623,3280,3282,615,2196,2199,2204,2216,2217,2224,2194,3260,3263,3262] 
set surface_color,  pcol23, surf_pocket23 
set_color pcol24 = [0.702,0.294,0.278]
select surf_pocket24, protein and id [1384,1423,1425,1764,1765,1802,1801,1392,1867,1866,1828,1829] 
set surface_color,  pcol24, surf_pocket24 
set_color pcol25 = [0.902,0.435,0.361]
select surf_pocket25, protein and id [7742,7746,5231,5233,4577,7747,7749,7753,7754,7730,7735,7121,7728,4595,4579,4580,4589,4590,7102,7105,7110,7118] 
set surface_color,  pcol25, surf_pocket25 
set_color pcol26 = [0.702,0.384,0.278]
select surf_pocket26, protein and id [7778,7781,8729,8731,7776,8431,8656,8103,8107,8745,8744,8432,8694,8695] 
set surface_color,  pcol26, surf_pocket26 
set_color pcol27 = [0.902,0.549,0.361]
select surf_pocket27, protein and id [1031,1094,1090,1040,1037,1032,952,958,964,966,967,969,875,1067,1091,1096,1089,9277] 
set surface_color,  pcol27, surf_pocket27 
set_color pcol28 = [0.702,0.475,0.278]
select surf_pocket28, protein and id [875,880,884,1098,1100,1116,1151,1149,1152,1154,2033,889,886,2067,2068] 
set surface_color,  pcol28, surf_pocket28 
set_color pcol29 = [0.902,0.667,0.361]
select surf_pocket29, protein and id [7778,7781,7130,8103,8104,8107,8431,8383,8395,8397,8432,8394] 
set surface_color,  pcol29, surf_pocket29 
set_color pcol30 = [0.702,0.561,0.278]
select surf_pocket30, protein and id [6799,6827,6824,6826,6829,6830,6833,6834,3057,3065,3066,3082,3080,2322,6800,6836,3452,3453,3071] 
set surface_color,  pcol30, surf_pocket30 
set_color pcol31 = [0.902,0.780,0.361]
select surf_pocket31, protein and id [4224,4093,4117,4119,4121,4241,4095,4258,4248,4255,4257,4253,4260,4264,4114,4115,4100,4112] 
set surface_color,  pcol31, surf_pocket31 
set_color pcol32 = [0.702,0.651,0.278]
select surf_pocket32, protein and id [6608,6610,6844,6862,5595,5631,5632,5633,6843,6594,6605,5592,5598] 
set surface_color,  pcol32, surf_pocket32 
set_color pcol33 = [0.902,0.894,0.361]
select surf_pocket33, protein and id [822,1111,838,877,1185,805,836,1119] 
set surface_color,  pcol33, surf_pocket33 
set_color pcol34 = [0.659,0.702,0.278]
select surf_pocket34, protein and id [8599,8527,8439,8563,8525] 
set surface_color,  pcol34, surf_pocket34 
   
        
        deselect
        
        orient
        