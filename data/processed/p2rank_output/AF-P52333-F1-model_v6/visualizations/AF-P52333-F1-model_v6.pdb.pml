
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
        
        load "data/AF-P52333-F1-model_v6.pdb", protein
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
 
        load "data/AF-P52333-F1-model_v6.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [5072,5220,5033,5070,5071,5203,4706,5087,5057,4699,4700,4703,4730,5223,4113,5059,5060,5063,5058,5227,4111,4123,5324,5327,5086,5047,5034,4130,4131,4318,4319,4105,4148,4086,4108,4114,4121,4098,4302,4653,4678,4672,4301,4685,4085,4089,4088] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.400,0.702]
select surf_pocket2, protein and id [6264,6239,6261,6262,6258,6259,6231,2512,2531,2532,2324,2464,2466,2467,2339,2481,2548,2453,2354,2355,2356,6238,6240,6244,4898,4900,4869,5153,4924,4922,4890,4891,4892,4926,6060,6081,6248,4520,4523,4533,5169,5171,5154,5111,6575,6292,6294,6260,4532,6603,6601,6599,6605] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.361,0.451,0.902]
select surf_pocket3, protein and id [7520,2289,2290,2291,3823,3825,2299,2496,2497,2498,2499,2501,6865,7521,6863,6864,2302,6900,2271,4980,4984,3810,3811,3820,7511,3821,7489,4990,4992,7496,7498,4993,4508,4955,4957,4499,2504,7516,7490,7491,7492,7538,2507,6899,6923,7039,7044,2503,4495,7064,7049,7037,7045] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.278,0.302,0.702]
select surf_pocket4, protein and id [1935,1964,1952,1953,1954,1938,1958,1904,1906,1908,775,778,1972,1959,1957,1970,599,799,801,803,804,807,1104,1105,1106,1129,1131,1127,575,1163,1162,1070,1095,770,762,765,601,572,573,1174,1159,1168,1169,1171,1173] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.396,0.361,0.902]
select surf_pocket5, protein and id [7454,7087,7446,7085,7460,7553,7555,7557,7559,6506,6504,6507,6646,7035,6631,6648,6647,6456,6457,6468,6469,6494,6483,6487,6491,7080,7082,6451,7079,6453,6450,7054,7476,7550,6891,7475,7552,7057,7052,7034,7061] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.357,0.278,0.702]
select surf_pocket6, protein and id [503,505,506,507,2718,2720,2705,3747,794,795,2676,626,2708,2709,2712,2831,2829,627,2704,2706,3749,2941,2931,2932,2883,2884,2851,2854,2747,3746,501,624,376,520,496,502,470,2940,498,2843,2850,3124] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.522,0.361,0.902]
select surf_pocket7, protein and id [2523,6919,6286,6290,6288,6289,6291,6292,6294,6920,6923,2504,4516,4519,4520,4523,6324,6326,6327,6546,6548,6550,6568,6553,6341,6611,6617,6626,6353,6355,6354,6608,2510,2516,2511,2512,6262,2521] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.455,0.278,0.702]
select surf_pocket8, protein and id [3877,3880,3881,3882,5011,5012,4463,5010,5447,5017,3885,4997,5005,5440,5443,5235,5246,5428,4435,4441,4442,4444,5231,5233,5222,5217,3901,3892,3893,3897,4439,4438,4440,5253,3908,3910,3912,4486,4487,4457,4458,4460,5444] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.647,0.361,0.902]
select surf_pocket9, protein and id [8613,8640,8580,3788,8579,8617,957,958,8619,3818,3816,3819,3783,3785,3787,3777,3801,943,950,955,978,927,930,941,956,903,904,3771,3727,923,925,926] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.549,0.278,0.702]
select surf_pocket10, protein and id [5972,5974,5894,5895,5896,6043,6044,5964,3681,5955,2890,2889,2768,2776,2778,6036,2751,2752,2757,2753,2758,2759,2891,2866,2761,2765,6009,3700,2741,2742,5978,5983] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.773,0.361,0.902]
select surf_pocket11, protein and id [6448,6437,6454,6524,7056,6442,6430,3980,6451,7069,3978,7153,4555,4561,4566,4556,4476,3961,3974,3976,3977,4569,4576,4572,3971,4477,4492,4469] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.647,0.278,0.702]
select surf_pocket12, protein and id [8613,8612,8609,8636,8640,3789,7276,8582,8583,8579,8581,8617,3818,3819,3783,3787,3793,3801,7251,7308,3812,3832,925,926,3771] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.898,0.361,0.902]
select surf_pocket13, protein and id [7719,7722,6488,6713,6714,6711,6715,6717,7706,7705,7710,7715,6483,6486,6484,6665,7574,7576,7582,6786,6741,6747,6748,7415,7689,7595,7596,6744,6746,6742] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.278,0.655]
select surf_pocket14, protein and id [2248,2258,2260,2524,2529,2169,2167,2171,2025,2172,2534,2541,2026,2027,2539,2166,2016,2168,6326,6328,6335,2523] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.361,0.776]
select surf_pocket15, protein and id [520,528,530,532,397,399,500,501,624,376,622,496,498,794,626,2831,2829,627,503,505,506] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.702,0.278,0.557]
select surf_pocket16, protein and id [6151,6144,5799,5578,5579,5581,4852,4809,4822,4819,6153,5815,6145,6143,5821] 
set surface_color,  pcol16, surf_pocket16 
set_color pcol17 = [0.902,0.361,0.651]
select surf_pocket17, protein and id [5045,5346,5350,5351,5352,5047,5355,5053,5566,4719,4759,5034,5322,5323,5320,5313,5356,5608,5059,4723,5056,4762] 
set surface_color,  pcol17, surf_pocket17 
set_color pcol18 = [0.702,0.278,0.459]
select surf_pocket18, protein and id [2255,2492,2248,2249,2524,2169,6852,6848,6851,2243,6843,6893,6905,6901,6904,6838,6917,6907,6895,6897,6898,6899,6923,2523,6919] 
set surface_color,  pcol18, surf_pocket18 
set_color pcol19 = [0.902,0.361,0.525]
select surf_pocket19, protein and id [5958,2912,2914,3675,2901,3551,3555,3556,3550,3552,3660,3091,3523,3519,3494,3520,3490,3493,3525,3526] 
set surface_color,  pcol19, surf_pocket19 
set_color pcol20 = [0.702,0.278,0.361]
select surf_pocket20, protein and id [983,971,2955,2956,3738,478,985,3139,3149,479,489,491,3138] 
set surface_color,  pcol20, surf_pocket20 
set_color pcol21 = [0.902,0.361,0.400]
select surf_pocket21, protein and id [4975,5468,5469,5460,4973,4972,5467,5451,3853,3858,3844,3847,3687,3689,3698,3699,3701,3690,4997] 
set surface_color,  pcol21, surf_pocket21 
set_color pcol22 = [0.702,0.298,0.278]
select surf_pocket22, protein and id [221,223,1121,1123,557,559,270,1150,1152,214,220,206,208,207,216,742,210,226,1089,1092,1130,1164,1167,760,765,741] 
set surface_color,  pcol22, surf_pocket22 
set_color pcol23 = [0.902,0.443,0.361]
select surf_pocket23, protein and id [7106,7776,8031,7102,7988,7448,7775,8066,8079,8029,8030,8032] 
set surface_color,  pcol23, surf_pocket23 
set_color pcol24 = [0.702,0.392,0.278]
select surf_pocket24, protein and id [6865,6858,6863,6864,6844,6852,6853,2301,2302,2255,2492,2248,2249,6846,2262,2263,2264,2268,2141,2247,3805,2271,6899] 
set surface_color,  pcol24, surf_pocket24 
set_color pcol25 = [0.902,0.569,0.361]
select surf_pocket25, protein and id [4574,4655,4657,4659,4661,6594,6595,4667,4663,4277,4273,4279,4282,4283,6416,6429,6406,6409,6541,6596] 
set surface_color,  pcol25, surf_pocket25 
set_color pcol26 = [0.702,0.490,0.278]
select surf_pocket26, protein and id [609,611,641,642,2621,620,2595,2597,2599,588] 
set surface_color,  pcol26, surf_pocket26 
set_color pcol27 = [0.902,0.694,0.361]
select surf_pocket27, protein and id [4456,4486,4487,3909,7190,7163,7506,3940,3959] 
set surface_color,  pcol27, surf_pocket27 
set_color pcol28 = [0.702,0.588,0.278]
select surf_pocket28, protein and id [3163,3176,3353,3357,2846,2845,2875,2850,2847,3128] 
set surface_color,  pcol28, surf_pocket28 
set_color pcol29 = [0.902,0.820,0.361]
select surf_pocket29, protein and id [2928,3012,3732,3734,2948,3724,3716,3717,3719,3714,3715,2748,2754,2755,2757,3749,3688] 
set surface_color,  pcol29, surf_pocket29 
set_color pcol30 = [0.702,0.686,0.278]
select surf_pocket30, protein and id [1740,1741,1742,1410,1411,1406,1408,1409,1405,1395,1739,1767,1704,1738,1769] 
set surface_color,  pcol30, surf_pocket30 
set_color pcol31 = [0.855,0.902,0.361]
select surf_pocket31, protein and id [8106,8307,8312,7924,8302,8445,8443,8446,8444] 
set surface_color,  pcol31, surf_pocket31 
   
        
        deselect
        
        orient
        