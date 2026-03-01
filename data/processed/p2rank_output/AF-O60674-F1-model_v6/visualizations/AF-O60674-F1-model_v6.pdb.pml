
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
        
        load "data/AF-O60674-F1-model_v6.pdb", protein
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
 
        load "data/AF-O60674-F1-model_v6.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [5464,5465,5466,6842,4463,4469,6839,7046,7047,7048,4459,4513,4456,4450,4451,4454,4482,4472,4473,4475,4477,4478,4479,4495,4486,5663,5700,4694,5474,5489,5636,5638,5490,5097,5100,5073,5096,5662,5658,5438,5655,4693,4512,4677,5048,4453,4676,5067,5791,5793,5798,5452,5795,5799,6016,5462,5461,5463,5103,5105,5765,5766,5767,5770,5678,5437,5439,5473,6060,6061,6062,6825,6840,6841,6822,5753,5754,5755,5756,5757,5759,6834,6835,6838,7438,7435,7436,6057,6079,6059,6139,6082] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.365,0.702]
select surf_pocket2, protein and id [491,3158,488,490,522,613,521,748,750,581,3468,3470,3472,610,497,3482,547,3150,3151,3155,3156,3183,3680,3682,3715,3681,3500,3503,3028,3144,3146,3061,3064,3166,3063,4111,3248,3249,4112,620,3037,615,617,618,619,582,2996,3159,3187,3190,3191,3433,3240,3431,3071,3074] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.361,0.361,0.902]
select surf_pocket3, protein and id [6949,6951,6953,6948,6950,6952,8086,8087,8088,8063,8064,8065,7922,7961,8235,7921,7950,8219,7118,7119,7263,7188,7191,7219,7224,7227,7231,7218,8078,8079,8080,7223,7220,6945,6954,6946,6935,6943,6957,6922,6924,6926,6931,6932,6919,6955,7122,7124,7136,6956,6970,7116,7117,7190,8103,8104] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.365,0.278,0.702]
select surf_pocket4, protein and id [6967,6970,7117,6913,6969,6916,6917,7531,7536,7540,6912,6914,7564,7567,7595,7599,7563,6918,6919,7102,7514,7519,8061,8063,8064,8065,7958,7955,7961,7977,7947,7570,7118,7119,8059,7513,8058,7368,7976,7369,8073,8055,8056,8075] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.576,0.361,0.902]
select surf_pocket5, protein and id [6834,6832,6833,4465,4466,4467,7489,7490,4741,7438,7439,7472,7471,7473,7432,7434,7435,4740,4730,6978,6979,7127,7128,7126,6876,6878,6880,6882,6883,6993,7123,7133,4713,4715,4716,4717,4503,4504,4505,4471,4476,4480,4481,4738] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.533,0.278,0.702]
select surf_pocket6, protein and id [6431,6436,6462,3082,4063,4064,6496,6498,4058,4046,3080,3059,6454,3093,3174,3199,3096,6490,6497,3117,3118,3069,3193,3198,3200,6339,6341,6342,6315,6286,6425,6427,6429,6406,6417,6307,3119,3196,3201,3222] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.792,0.361,0.902]
select surf_pocket7, protein and id [6240,6241,6007,6669,6166,6704,6070,6071,6077,6078,6064,6008,6058,6083,6746,6717,6687,6688,6745,6094,6096,6087,6735,6084] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.702]
select surf_pocket8, protein and id [907,1245,661,902,1289,696,699,1343,722,724,697,1331,2159,2187,2164,2172,2173,915,2195,2197,2193,2191] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.361,0.792]
select surf_pocket9, protein and id [1272,1273,1275,1267,1268,1271,1316,1230,2054,1511,1514,1520,1522,1543,1544,1545,322,1234,1236,335,337,338,339,362,371,368,369,378,1237,1278,1200,2053] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.278,0.533]
select surf_pocket10, protein and id [5269,6537,5289,5290,6534,6535,2672,2673,2795,2842,2656,2796,2777,2778,2859,2671,2654,5326,5267,5268] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.361,0.576]
select surf_pocket11, protein and id [665,667,734,737,853,2916,870,861,783,769,774,778,788,856,792,855,2913,768] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.702,0.278,0.365]
select surf_pocket12, protein and id [4100,3279,591,4103,1129,588,590,601,603,1137,1136,3443,3261,3263,3445,3289,1116,1073] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.361,0.361]
select surf_pocket13, protein and id [2045,1967,1969,1970,1971,1502,1992,2005,1451,1479,1470,1449,2068] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.365,0.278]
select surf_pocket14, protein and id [3445,3271,3446,3641,3427,3439,3513,3514,3450,3452,3453,3454,3512] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.576,0.361]
select surf_pocket15, protein and id [4524,4520,4522,4444,6877,4428,4503,4505,6978,6876,6878,6879,6880,6881,6882,6883] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.702,0.533,0.278]
select surf_pocket16, protein and id [710,2942,697,718,2944,709,2203,2195,2193,2191,2198,2214,2216,2218,2922,2924] 
set surface_color,  pcol16, surf_pocket16 
set_color pcol17 = [0.902,0.792,0.361]
select surf_pocket17, protein and id [3204,3221,3894,3833,3859,6411,3868] 
set surface_color,  pcol17, surf_pocket17 
set_color pcol18 = [0.702,0.702,0.278]
select surf_pocket18, protein and id [4981,4347,4349,4350,4332,4333,4314,4315,4813,4815,4812,4814] 
set surface_color,  pcol18, surf_pocket18 
   
        
        deselect
        
        orient
        