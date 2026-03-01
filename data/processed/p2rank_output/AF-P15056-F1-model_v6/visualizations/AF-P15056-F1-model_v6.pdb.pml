
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
        
        load "data/AF-P15056-F1-model_v6.pdb", protein
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
 
        load "data/AF-P15056-F1-model_v6.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [3922,3919,3920,3923,4038,4039,4044,3660,3672,2810,3578,3579,4568,2809,3675,4572,4573,4596,4564,4567,4583,3844,3843,4483,4484,4486,2817,4487,4488,4485,4063,4024,4025] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.302,0.278,0.702]
select surf_pocket2, protein and id [4623,4625,4733,4736,4601,4426,4430,3560,3561,4428,4429,4447,4416,4731,4738,3558,3676,3677,4574,3541,3565,3542,2802,2792,2794,3579] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.631,0.361,0.902]
select surf_pocket3, protein and id [2106,2083,2104,2105,2041,2042,1963,1965,1967,2082,2060,1771,1775,1772,1777,2112,4370,4385] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.678,0.278,0.702]
select surf_pocket4, protein and id [4706,4716,4705,4729,4731,4887,4888,4890,4710,4780,4781,4782,4738,4601,4426,3561,4623,4625,4733,4428,4429,4416,4620,4708,4602] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.361,0.682]
select surf_pocket5, protein and id [1758,2006,2012,2018,2020,2113,2137,2143,2140,2161,2163,1879,2005,2159,2155,1877,2111,2127] 
set surface_color,  pcol5, surf_pocket5 
set_color pcol6 = [0.702,0.278,0.341]
select surf_pocket6, protein and id [2221,2223,2219,2253,2254,2255,2256,2239,2240,2283,2296,1695,1696,2243,1552] 
set surface_color,  pcol6, surf_pocket6 
set_color pcol7 = [0.902,0.522,0.361]
select surf_pocket7, protein and id [4111,4113,4740,4744,4774,4739,5056,4750,4759,4760,4752,5150,5116,4758,5125,5118,5120] 
set surface_color,  pcol7, surf_pocket7 
set_color pcol8 = [0.702,0.596,0.278]
select surf_pocket8, protein and id [1189,1196,1199,1207,1209,3931,3894,3899,1201,4047,4503,4506,4514,4515,4046,4048,4554,4537,4553,4538] 
set surface_color,  pcol8, surf_pocket8 
   
        
        deselect
        
        orient
        