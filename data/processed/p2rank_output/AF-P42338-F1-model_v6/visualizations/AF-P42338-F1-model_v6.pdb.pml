
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
        
        load "data/AF-P42338-F1-model_v6.pdb", protein
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
 
        load "data/AF-P42338-F1-model_v6.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [5462,5463,5464,5465,5466,5467,6159,6162,6169,5432,5435,6171,6807,6812,6803,6806,5460,5155,5154,5156,6585,5189,5509,5510,6615,6617,6588,6607,6610,6799,6800,1453,6172,6784,6644,6647,5164,6643] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.380,0.702]
select surf_pocket2, protein and id [7521,7524,7525,7526,6934,7416,7442,7443,6372,6892,6918,6373,6956,6284,6894,6900,6876,6882,6888,6496,6796,6798,7529,7533,6318,6303,6304,6516,6515,6350,6877] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.361,0.400,0.902]
select surf_pocket3, protein and id [8121,5248,5245,8129,3850,8126,8128,8130,3849,3852,3853,3838,3795,5580,3741,3735,3737,3837,3855,3858,3863,3559,3731,3763,3742,3746,3761,3557,5587,5246,5568,5576] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.318,0.278,0.702]
select surf_pocket4, protein and id [8071,8099,4403,4919,4920,8076,8123,4414,8192,8191,4404,4408,4977,4921,4672,4945,4946,4964,4950,2924,3852,3853,4953,4942,4944,4965,8068] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.498,0.361,0.902]
select surf_pocket5, protein and id [5822,5823,5824,6078,6080,6057,156,5847,5848,807,367,370,806,808,135,830] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.455,0.278,0.702]
select surf_pocket6, protein and id [1149,3879,3880,3886,3869,1132,1133,1139,1140,3871,1144,3588,3583,3668,3579,3582,3612,3617,5287,3623,3625,5276,3868,3586] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.675,0.361,0.902]
select surf_pocket7, protein and id [3459,3478,3513,3514,2856,2858,2876,2987,3460,2984,2986,2854,2992,3000,3005,3007,3480] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.592,0.278,0.702]
select surf_pocket8, protein and id [6708,6671,6709,6710,6705,1472,1473,2194,2164,2166,2204,2205,2206,2191,6785,6786,6768,6707,6744,1500,1525,6743,2223,2225,2234,2227,7487,6754] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.851,0.361,0.902]
select surf_pocket9, protein and id [5712,5746,5737,5738,5424,2410,1385,2371,2408,2370,2372,2413,1418,1428,1429,1431,2374,2053,2054,6142,6143,6144,6201] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.278,0.671]
select surf_pocket10, protein and id [8234,8236,8238,7653,8264,8231,8258,8259,8260,8267,8265,7596,7597,7693,7695,8404,8405,7341,7343,7345,8366,7696] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.361,0.773]
select surf_pocket11, protein and id [3214,4996,4997,4999,5000,5001,5037,5039,3206,3630,3632,5271,5282,5286,5302,5312,5315,5278,5284,3621,5307,5309,5024] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.702,0.278,0.533]
select surf_pocket12, protein and id [6135,6208,6400,6156,6157,5831,6380,6487,6113,6868,5798,5805,5807,5800] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.361,0.596]
select surf_pocket13, protein and id [2094,2095,2097,1784,1786,1788,1970,1972,1811,1988,2086,1777,1787,2104,1364] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.278,0.396]
select surf_pocket14, protein and id [1020,1022,5760,5758,812,778,1023,997,721,783,780,779,756,1001,813,5814,722,821,5756,5786,5789,5759,5791,5794,5792,5793] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.361,0.420]
select surf_pocket15, protein and id [688,689,700,672,1054,1072,1041,1042,1015,1016,986,1017,1023,991,994,995,1018,1021,993,719,706] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.702,0.302,0.278]
select surf_pocket16, protein and id [1407,1419,2412,1379,1383,5701,5417,5421,5424,1316,2440,2442,5702,1406,1408,1412,5398,1317] 
set surface_color,  pcol16, surf_pocket16 
set_color pcol17 = [0.902,0.478,0.361]
select surf_pocket17, protein and id [6541,6552,6580,6532,6534,6535,6516,6515,6536,6350,6316,7546,7545,7531,7532,7533,6318,6579] 
set surface_color,  pcol17, surf_pocket17 
set_color pcol18 = [0.702,0.439,0.278]
select surf_pocket18, protein and id [6826,6830,5500,5803,5806,5802,5804,5815,5783,5784,6857,5839,5841,5836,6833,5526,5529,5814,1021] 
set surface_color,  pcol18, surf_pocket18 
set_color pcol19 = [0.902,0.655,0.361]
select surf_pocket19, protein and id [6470,6267,6469,6906,6372,6892,6265,6370,6388,6478,6479,6387] 
set surface_color,  pcol19, surf_pocket19 
set_color pcol20 = [0.702,0.576,0.278]
select surf_pocket20, protein and id [3758,3795,3761,3553,3557,3754,3753,3837,3559,3834,3819,3820] 
set surface_color,  pcol20, surf_pocket20 
set_color pcol21 = [0.902,0.831,0.361]
select surf_pocket21, protein and id [8041,4898,4900,4614,4616,4649,4369,4405,8038,4897,4899,4639,4615,4921,4896,4893,4895,4922,8076] 
set surface_color,  pcol21, surf_pocket21 
set_color pcol22 = [0.686,0.702,0.278]
select surf_pocket22, protein and id [7848,7852,7804,7832,7805,7875,7878,7880,7882,7770,8445,8446,8442,7771,7741,7739] 
set surface_color,  pcol22, surf_pocket22 
   
        
        deselect
        
        orient
        