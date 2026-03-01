
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
        
        load "data/AF-P48736-F1-model_v6.pdb", protein
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
 
        load "data/AF-P48736-F1-model_v6.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [5318,5331,5332,5298,6879,6881,6882,6883,1652,1654,1655,5319,6858,6859,5322,5632,5633,5634,5636,5635,5637,6913,6915,6943,6947,2369,2373,2370,2372,2374,2399,2402,6948,1734,1727,6826,6822,5679,5681,6886,6411,7021,7022,6413,6911,2415,2437,2438,2439,2447,2398,2407,2408,2410,2412,2414,2400,2434,7040,7044,6395,6397,6398,6407,6409,6401,5296,5317,5323,5330,1682,1683,5299,1668,6906,6909,7005,6910,6820,6823,6846,6880,6916,6912,6905,6907,5600,5601,6412,2258,6420,2384,2388] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.365,0.702]
select surf_pocket2, protein and id [7035,7769,7772,7773,7774,7692,7693,7161,7114,6732,7120,7126,6591,7131,7137,7133,7139,7775,7677,7174,7665,7666,7115,6568,6731,6752,6535,6751,7781,6533,7777,7779,6524,6817,7033] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.361,0.361,0.902]
select surf_pocket3, protein and id [5138,5139,5140,5126,3186,8430,4545,4548,4793,4552,4563,4553,4555,4547,5097,5099,5102,4781,4780,4782,4816,4818,4580,4579,4790,8353,8356,8420,4558,8455,8422,8423,8427,8428,8429,8418,8364,8308,8331,8333,8299,3185,4032,4034,4036,4038,4033,8411,8413] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.365,0.278,0.702]
select surf_pocket4, protein and id [5171,5172,5173,5174,3815,3477,3841,3842,2704,2674,2677,2669,2691,2662,5470,5503,5506,5467,5469,1474,1476,5471,5449,5472,5473,5445,5468,3814,1450,1441,1439,5474,5478,2679] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.576,0.361,0.902]
select surf_pocket5, protein and id [6417,6418,2266,6379,6441,2265,5948,5949,5588,5591,6419,1634,5922,2574,5593,1616,2535,2543,2562,2565,5912,5914] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.533,0.278,0.702]
select surf_pocket6, protein and id [4088,4090,4091,4087,4089,4092,4094,4096,1389,1388,4102,3785,1390,5462,5463,5794,1370,1379,1375,1377,1381,1415,1410,1413,1397,1398,5789] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.792,0.361,0.902]
select surf_pocket7, protein and id [6752,7780,6535,6751,7781,7779,7853,7855,7628,7632,8780,8817,7778,7676,7677,7665,7653,7652,7794,7796,7803,7806,6779,6818,6550,6774,6544,6772,6773,6770] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.702]
select surf_pocket8, protein and id [7187,7220,7222,7280,7282,8788,8790,8782,8784,7668,7185,7188,8824,7669,8747,7657,8781,8813,8816,8817] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.361,0.792]
select surf_pocket9, protein and id [5305,5306,5307,6917,1683,5066,5342,5344,5274,5326,5328,6915,8260,8259,6914,6916,6918,6919,6954,6923,5063,8284,6893,5051,5053,5060,5048] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.278,0.533]
select surf_pocket10, protein and id [6403,6405,6702,6431,6711,6706,7123,6694,7145,7147,6697,7122,7124,7121,6710,6634,6432,2245,2247,6632,6691,6692,6693,6665,6679,6678,2242,6667,6675,6663,7021,6413,2415,7037,6398,6409,6401,2258,2253] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.361,0.576]
select surf_pocket11, protein and id [7824,7842,7844,7846,8472,8474,7914,7932,7934,7587,8476,7610,7611,7918,7941,7825,7826,7899,7592,8647,7940,8648,7938,8499,8469,8498,8504,8505,8612,8610] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.702,0.278,0.365]
select surf_pocket12, protein and id [5743,5741,5739,5744,4056,4060,5762,5763,5764,3752,3755,5747,4046,3969,4001,3970,3957,3973,3975,5733,4066,8370,8374,8375,8376,6801] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.361,0.361]
select surf_pocket13, protein and id [5166,3472,3475,5169,3469,3463,5174,3840,3842,4301,4304,4286,4298,4290,4307,4899,4875,4281,4933,5165,5181,4931] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.365,0.278]
select surf_pocket14, protein and id [4068,3526,4014,4069,3527,3162,3202,3161,3163,3509,3751,4064,4065,4066,3747,3749,4077,4044,4035,4037,3177,3179,3183,3182] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.576,0.361]
select surf_pocket15, protein and id [1624,2211,2299,2301,2206,2212,2292,1971,1975,1976,1998,2198,2528,2530,1589,1599,2303,2304] 
set surface_color,  pcol15, surf_pocket15 
set_color pcol16 = [0.702,0.533,0.278]
select surf_pocket16, protein and id [2865,2899,2866,2778,4250,2908,4229,4232,4266,4268] 
set surface_color,  pcol16, surf_pocket16 
set_color pcol17 = [0.902,0.792,0.361]
select surf_pocket17, protein and id [3113,3115,3117,3264,3097,3099,3703,3706,3690,3243,3247,3252,3267,3266,3554,3555] 
set surface_color,  pcol17, surf_pocket17 
set_color pcol18 = [0.702,0.702,0.278]
select surf_pocket18, protein and id [2890,2894,4878,4904,4882,2886,2888,2731,2884,4308,2745,4305,4870,4872,4873,4877,4299,4905,4941,4312,4329,4330,4311] 
set surface_color,  pcol18, surf_pocket18 
   
        
        deselect
        
        orient
        