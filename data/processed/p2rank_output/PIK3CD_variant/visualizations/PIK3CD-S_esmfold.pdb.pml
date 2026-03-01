
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
        
        load "data/PIK3CD-S_esmfold.pdb", protein
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
 
        load "data/PIK3CD-S_esmfold.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [7069,7074,7079,7080,7062,7065,7068,7245,7266,7667,7244,7514,7512,7666,7508,7545,7547,7541,7543,7072,8144,7458,7460,7461,7466,7476,7724,7687,7658,7664,7685,7686,7695,7700,7660,7690,7657,7546,7692,7263,7265,6511,6510,6516,7249,7247,6512,7242,7238,7241,7228,6374,6375,6436,6491,6489,6486,6488,6492,7260,7261] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.380,0.702]
select surf_pocket2, protein and id [4888,5201,5203,1386,4890,4892,5171,5172,5205,5919,6539,6628,5926,5928,5929,4886,5248,4891,1387,2106,2107,5247,6345,6347,6379,6381,6382,6522,6383,6384,6523,6356,6358,4894,4898,6352,6546,6346,6530,6538,6542,4925,4929,6325,6324,6351,6354] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.361,0.400,0.902]
select surf_pocket3, protein and id [6618,6619,7275,7276,7278,7282,7274,6235,6118,6119,6681,6318,6532,6534,7189,7270,7273,6624,6636,6642,7255,6640,6666,7190,6049,6050,6062,6064,6255,6254,6275,6096,6634] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.318,0.278,0.702]
select surf_pocket4, protein and id [4737,3440,3441,5038,5017,5019,5021,3417,5060,5045,5043,5048,3096,4720,4723,3088,5005,5006,5007,4990,4758,4759,4687,4712,4730,4735,4739,4676,4679,4682,5057,4760,4761,4788,5051,4738,4764,4748,5053] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.498,0.361,0.902]
select surf_pocket5, protein and id [101,5812,304,312,98,100,736,759,87,102,5623,119,120,118,5811,5625,5559,5588,5589,5590,5784,5838,5567,5569,5834,5836,5825,5560,5591] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.455,0.278,0.702]
select surf_pocket6, protein and id [1402,1403,2105,2106,2107,2102,2103,2135,2137,2139,2157,1401,1404,1438,2131,2159,2161,2162,2133,6408,6445,6446,6494,6483,6484,6499,2163,1473,1435,6524,6525,6508,7236,6383,6384] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.675,0.361,0.902]
select surf_pocket7, protein and id [3674,3675,3676,3679,3466,3467,3464,3495,3498,3376,3465,3512,3494,3366,3666,3680,3683,3685,1073,1074,1075,1089,1091,5020,5022,3406,3408,3412,3413,3419,5031,5032,5012,3668] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.592,0.278,0.702]
select surf_pocket8, protein and id [2013,2015,1953,1954,1992,5935,1990,1993,1991,1994,1351,1359,1350,1361,1363,2328,2288,2291,2292,2320,5934,5163,5158,2234,5900,5458,5459,5460,5901] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.851,0.361,0.902]
select surf_pocket9, protein and id [958,960,999,1017,1018,959,647,933,645,937,935,653,991,998,1027,599,635,641,630,1005,1010,1013,927,928,968,648] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.278,0.671]
select surf_pocket10, protein and id [6201,6202,1966,6183,6184,6200,1981,1984,1980,2137,2139,2145,2115,2117,2119,2177,2180,2186,2203] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.361,0.773]
select surf_pocket11, protein and id [7346,7348,7404,7349,7421,7441,7330,7998,8000,8002,8171,7090,7092,8029,8031,8132,8170,7439,7445,7448,7447,8022,8023,8024,8027,8028] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.702,0.278,0.533]
select surf_pocket12, protein and id [5314,3542,3546,5318,5320,3540,3655,3658,3660,5306,4981,4982,5326,3648,4984,7897,3635,3652,3355,3551] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.361,0.596]
select surf_pocket13, protein and id [1387,1392,4863,2106,1410,1411,4861,4860,4888,4890,5172,6381,6382,6383,6384,4898,4900,6415] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.278,0.396]
select surf_pocket14, protein and id [2031,2036,2037,2038,2039,1295,1296,1297,1751,1923,2027,1719,1723,1906,1726,1729,2045,1325,1289,1291,1934,1937,1293,1298,1333] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.361,0.420]
select surf_pocket15, protein and id [7981,7828,7864,7980,7826,7833,7892,4686,4688,4690,4717,4707,8018,8019,7954,4663,4665,7827,7830,7838] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.702,0.302,0.278]
select surf_pocket16, protein and id [4663,4665,4419,4391,7838,7841,4149,4393,4686,4688,4690,4717,4707,4718,7864,8019] 
set surface_color,  pcol16, surf_pocket16 
set_color pcol17 = [0.902,0.478,0.361]
select surf_pocket17, protein and id [6596,6599,5549,5550,5553,6569,5532,5530,5238,6567,6571,973,5262,5264,5267,5555,5556,5551,6598,5584] 
set surface_color,  pcol17, surf_pocket17 
set_color pcol18 = [0.702,0.439,0.278]
select surf_pocket18, protein and id [4913,6391,4901,4900,6414,6454,6455,7812,6423,7790,7791,4870,4872,4627,4629,4634,7811,4911,4630,4841,4638] 
set surface_color,  pcol18, surf_pocket18 
set_color pcol19 = [0.902,0.655,0.361]
select surf_pocket19, protein and id [6610,6611,5576,5524,5545,5891,6147,5868,6226,5965,6227,5913] 
set surface_color,  pcol19, surf_pocket19 
set_color pcol20 = [0.702,0.576,0.278]
select surf_pocket20, protein and id [6255,6254,6271,6273,6274,6275,6096,6062,6064,7280,7281,7282,6319] 
set surface_color,  pcol20, surf_pocket20 
set_color pcol21 = [0.902,0.831,0.361]
select surf_pocket21, protein and id [7555,7556,7583,7631,7632,7634,7636,8209,7493,8208,8211,8212,7521,7522] 
set surface_color,  pcol21, surf_pocket21 
set_color pcol22 = [0.686,0.702,0.278]
select surf_pocket22, protein and id [4572,4747,3884,3885,4564,4566,4563,4568,4743,4745,4744] 
set surface_color,  pcol22, surf_pocket22 
   
        
        deselect
        
        orient
        