
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
        
        load "data/AF-P00519-F1-model_v6.pdb", protein
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
 
        load "data/AF-P00519-F1-model_v6.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [2940,2513,2536,2537,2539,2542,2114,1949,1954,2003,2005,2006,2939,1952,2126,1948,1950,1955,2283,2357,2492,2486,2487,2471,2130,3021,1962,1963,3020,3022,1990,2129,2920,2926,2912,3012,3015,3016,3018,3013,2257,2506,2507,2508,2512,2535,1953,1982,1975] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.302,0.278,0.702]
select surf_pocket2, protein and id [1014,1018,1024,1023,1022,1029,1052,2711,2713,2715,1000,999,1003,1085,1086,1087,1724,2985,1707,1718,1721,1723,1699,2328,2764,1197,3910,3911,2762,2342,2343,2737,2741,2744,2733,2734,2708,2731,2732,2736,2332,2334,1196,1198,3940,3942,1027,3944,2786] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.631,0.361,0.902]
select surf_pocket3, protein and id [1730,884,1732,549,551,568,570,2095,566,2074,2075,2955,2957,1721,1739,1743,1747,1748,1749,3002,1753,2383,2096,1777,1755,1752,1771,1773,2496,2514,2952,2986,2489,2491,2497,2501,2101,543,553,557,538,1734,2362,2364] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.678,0.278,0.702]
select surf_pocket4, protein and id [2670,2672,4039,3681,4038,4073,2702,2696,2697,4009,4014,2727,2690,2691,2692,3428,2694,3451,3453,3683,3685,3687,3694,3673,3680,3934,3718,3719,3395,3430,3431,3433,3434,2662,2664,2669,2629,4011] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.361,0.682]
select surf_pocket5, protein and id [2288,2317,2318,724,725,732,1838,1842,1837,1839,1841,1822,1812,1813,1820,2316,745,722,729,1809,2299,2300,2301,714,1783,1796,1823,746,1760,1824,1782] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.278,0.341]
select surf_pocket6, protein and id [1174,1185,1189,1245,1247,1033,1035,1036,1339,1184,1187,1225,1358,1361,1363,1021,1025,1177,1179,1030,3980,3977,3976,1027,3941,1204,1210] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.522,0.361]
select surf_pocket7, protein and id [935,963,1526,425,409,435,436,918,952,953,1514,1509,1511,1513] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.596,0.278]
select surf_pocket8, protein and id [3179,3157,3168,2881,3062,3151,3170,3116,3125,3127,3135,3044,2875,3045,2220,2874] 
set surface_color,  pcol8, surf_pocket8 
   
        
        deselect
        
        orient
        