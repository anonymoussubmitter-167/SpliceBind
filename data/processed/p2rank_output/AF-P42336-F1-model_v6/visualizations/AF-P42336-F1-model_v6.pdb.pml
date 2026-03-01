
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
        
        load "data/AF-P42336-F1-model_v6.pdb", protein
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
 
        load "data/AF-P42336-F1-model_v6.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [5102,5104,5106,5107,5108,5076,6686,6687,6690,1409,2217,2219,1408,6662,5135,5137,6625,6627,6660,6658,6628,5105,6831,6838,6850,6853,6854,6847,6931,5383,5417,5419,5420,6163,6170,5098,5100,5097,5099,1393,5384,5418,6862,5413,5415,5416,5466,6631,6685,6723,6655,6650,6656] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.361,0.702]
select surf_pocket2, protein and id [6945,7564,7579,7582,6969,6970,7500,7501,7583,6986,7011,7012,6845,6943,6368,6933,6937,6939,6365,6369,6539,6927,6921,6622,7584,7585,7587,7589,7590,7591,7476,6303,6315,6346,6538,6922,6557,6559] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.373,0.361,0.902]
select surf_pocket3, protein and id [3638,3634,3636,3640,8216,8224,8210,8225,8226,8227,8263,8265,8266,8267,3716,3712,3714,3718,8264,3702,3701,3661,3660,8215,5537,5539,3443,3614,5523,5524,5540,3608,3623,8217,8222,5189,5192,5546,3722,3727,3728,3719] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.380,0.278,0.702]
select surf_pocket4, protein and id [5456,6874,5779,5817,5816,52,5760,5756,5789,5791,30,33,5479,5502,993,998,5501,5503,5507,5508,1000,5483,42,47,38,5848,23,26,27,21,19,6879,6882,6885,6886,6899,6901,6565,6566,6567,5847,5488,5514,5510] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.600,0.361,0.902]
select surf_pocket5, protein and id [3744,3740,3743,1088,1086,5244,1102,1103,5237,5242,1099,5243,1087,5561,5562,3731,3736,3587,3737,3738,5220,3733,3735,3548,3574,3575,5231,5233,3493,3495,1104,3465,3468,3458,3462,3463,3469,3471] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.557,0.278,0.702]
select surf_pocket6, protein and id [6685,6713,6714,6715,6832,6833,6803,6805,6815,6721,6722,6723,2274,2277,2278,6753,2285,2286,6806,2284,6749,1409,2214,2217,2219,1407,1408,2215,2254,2256,6686,6690] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.831,0.361,0.902]
select surf_pocket7, protein and id [5687,5689,6144,6145,6179,6203,6205,5715,5717,2434,5675,5718,5375,5376,5680,2451,1359,1358,1368,1371,2094,2433,2415,2416,2428,2435,2420,2098,6204] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.667]
select surf_pocket8, protein and id [2220,2221,2254,2256,2262,2294,2297,2306,2090,6196,6478,6477,6472,6485,6487,1392,2106,2107,2223,2227,2229,6467,6195,6173,2285,2301] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.361,0.741]
select surf_pocket9, protein and id [6930,6952,6928,6929,6931,6936,6194,6195,6515,6488,6508,6512,6514,6503,6196,6477,6485,6487,6489,6151,6162,6165,6529,6531,6519] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.278,0.486]
select surf_pocket10, protein and id [708,709,710,742,743,686,5768,5769,5770,951,981,983,979,982,5729,5735,5741,957,5730,5732,5733,959] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.361,0.514]
select surf_pocket11, protein and id [1019,1021,1056,1023,939,1026,626,913,915,935,936,941,942,937,1016,632] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.702,0.278,0.310]
select surf_pocket12, protein and id [4067,4682,4059,4686,4465,4692,4694,4160,4667,4145,4066,4172] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.435,0.361]
select surf_pocket13, protein and id [5461,5489,6601,6602,5169,5497,6599,5460,5464,5136] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.427,0.278]
select surf_pocket14, protein and id [4039,4041,4042,4005,4007,4712,4740,4006,4015,4742,4504] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.667,0.361]
select surf_pocket15, protein and id [5799,64,5801,5824,6081,117,6073,5827,5797,5821,5826] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.702,0.608,0.278]
select surf_pocket16, protein and id [6913,6914,5783,6531,6398,6399,6376] 
set surface_color,  pcol16, surf_pocket16 
set_color pcol17 = [0.902,0.894,0.361]
select surf_pocket17, protein and id [6979,7211,7212,6980,7044,7171,7175,7042,7043] 
set surface_color,  pcol17, surf_pocket17 
   
        
        deselect
        
        orient
        