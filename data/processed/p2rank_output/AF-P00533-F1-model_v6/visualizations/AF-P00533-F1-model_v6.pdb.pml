
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
        
        load "data/AF-P00533-F1-model_v6.pdb", protein
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
 
        load "data/AF-P00533-F1-model_v6.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [5879,6589,6586,6588,6590,6487,5538,8303,8305,8304,5710,5711,5534,5542,5522,5523,5540,5941,5943,6501,6502,6582,8294,8296,6106,5707,5694,6054,6055,5558,5559,8291,8300] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.278,0.349,0.702]
select surf_pocket2, protein and id [2218,2221,835,836,2224,649,463,2213,2232,1959,1962,2231,2217,459,9069,460,660,9083,9084,850,2304,2377,2380,2361,9068,2375,188,448,187,207,447,184] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.404,0.361,0.902]
select surf_pocket3, protein and id [5511,5496,8253,6085,8283,7773,7736,7737,7738,6071,6072,6077,6082,5578,5580,5494] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.416,0.278,0.702]
select surf_pocket4, protein and id [1127,1135,1137,1138,1925,835,2224,1949,1933,1936,1937,1938,1939,1940,1951,850,849,851,9083,9109,9110,9112,9114,9085,9086] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.663,0.361,0.902]
select surf_pocket5, protein and id [6069,6072,5680,5682,7774,7862,7781,7782,7746,7747,7748,7750,7744,7745,6087,6084,6086,6517,7773,7734,7736,7855,7858,7860,7846] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.616,0.278,0.702]
select surf_pocket6, protein and id [2392,2504,2327,9104,2523,2505,2336,9058,9066,2366,2369,2385,2387,2406,2407,2408,2414,2393,2343,2346,2325,2368] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.361,0.878]
select surf_pocket7, protein and id [5538,6475,6476,6479,8303,8305,8315,8313,8326,6589,6586,6587,6588,6486,6487,6450,5540,5541,6744,6746,6749,6750,6752,8317,8318,8319] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.278,0.584]
select surf_pocket8, protein and id [7719,7731,7733,7734,7735,6096,6098,7678,7674,7682,7684,7685,7717,7710,6085,8287,6094,6171,6144,6178,6180,6181,6531,7693,7690,7683] 
set surface_color,  pcol8, surf_pocket8 
set_color pcol9 = [0.902,0.361,0.620]
select surf_pocket9, protein and id [247,370,488,337,338,336,487,492,233,288,360,361,362,364,365,367,333,507,244] 
set surface_color,  pcol9, surf_pocket9 
set_color pcol10 = [0.702,0.278,0.380]
select surf_pocket10, protein and id [3204,3206,3451,3987,4016,3980,3981,3198,4004,4015] 
set surface_color,  pcol10, surf_pocket10 
set_color pcol11 = [0.902,0.361,0.361]
select surf_pocket11, protein and id [2134,4594,4532,4533,2073,2075,4515,4548,4577,4575,4546,4563,4573,2133,2111,4593,4426,4541,4534,4530] 
set surface_color,  pcol11, surf_pocket11 
set_color pcol12 = [0.702,0.380,0.278]
select surf_pocket12, protein and id [9060,9062,9063,9065,9076,858,859,9053,857,9069,458,460,662,847,2380,240,474,219,464,465] 
set surface_color,  pcol12, surf_pocket12 
set_color pcol13 = [0.902,0.620,0.361]
select surf_pocket13, protein and id [4071,4073,4266,4088,4090,4091,4027,3476,4021,4002,4093,4031,4033,4036,4041] 
set surface_color,  pcol13, surf_pocket13 
set_color pcol14 = [0.702,0.584,0.278]
select surf_pocket14, protein and id [3563,3350,3338,3393,3349,3372,3136,3408,3562,3374] 
set surface_color,  pcol14, surf_pocket14 
set_color pcol15 = [0.902,0.878,0.361]
select surf_pocket15, protein and id [1125,1321,1323,1325,1111,1112,1716,1645,1646,1647,1644,1128,1927,1706,1707,1708,1715,1820] 
set surface_color,  pcol15, surf_pocket15 
   
        
        deselect
        
        orient
        