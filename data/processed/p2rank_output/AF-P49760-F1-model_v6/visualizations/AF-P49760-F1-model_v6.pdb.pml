
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
        
        load "data/AF-P49760-F1-model_v6.pdb", protein
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
 
        load "data/AF-P49760-F1-model_v6.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [1949,1950,2544,2785,2783,2786,2118,2099,1489,1492,1654,2110,1494,1549,2116,2087,2089,2090,2091,2092,2135,2132,2542,2543,2516,2123,1488,1490,2093,2094,1672,1528,1670,2787,2788,2790,2792,2793,1797,1547,1550,1667,1495,1496,1498,2528] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.416,0.278,0.702]
select surf_pocket2, protein and id [4094,4098,2749,4169,4170,2310,2312,2313,2314,2315,2316,4164,4166,4122,4123,4168,4172,4177,4181,4215,4083,4201,2279,2273,2274,2276,2302,2336,2338,2346,2347,2349,2340,4095,4128,4150,4135,4132,2372] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.902,0.361,0.878]
select surf_pocket3, protein and id [2945,2946,3219,3244,3246,672,674,683,686,681,673,3263,3264,715,716,2519,2925,2959,3165,3218,684,2152,2153,2154,2155] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.702,0.278,0.380]
select surf_pocket4, protein and id [1764,1735,1737,1521,1672,1798,1523,1524,1526,1688,2792,2793,2808,2815,2807,2794] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.620,0.361]
select surf_pocket5, protein and id [613,3209,3211,3692,3759,3217,3247,3660,3691,3701,3331,3242,3224,3239] 
set surface_color,  pcol5, surf_pocket5 
   
        
        deselect
        
        orient
        