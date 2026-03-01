
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
        
        load "data/AF-Q9HAZ1-F1-model_v6.pdb", protein
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
 
        load "data/AF-Q9HAZ1-F1-model_v6.pdb_points.pdb.gz", points
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
select surf_pocket1, protein and id [1470,1471,1587,1588,1589,8,1421,1422,1423,1425,1427,1449,1451,1997,1999,2001,2002,1998,2000,2414,2415,2426,2439,2704,2706,2708,2700,2701,2702,1861,2697,2716,2699,1862,1710,2455,2023,2025,2026,2427,4,1417,1572,2020,2007,6,7] 
set surface_color,  pcol1, surf_pocket1 
set_color pcol2 = [0.416,0.278,0.702]
select surf_pocket2, protein and id [640,642,643,607,608,650,679,683,2062,2064,2134,2130,2132,664,3169,3170,676,678,3148,3149,3150,3151,3152,3071,3120,3125,3112,3110,3564,3167,3147,3154,2860,2859] 
set surface_color,  pcol2, surf_pocket2 
set_color pcol3 = [0.902,0.361,0.878]
select surf_pocket3, protein and id [541,542,543,571,509,510,3117,3101,3118,3119,3104,3679,3598,535,531,534,3600,3587,3661,3617,3618,3619,3659] 
set surface_color,  pcol3, surf_pocket3 
set_color pcol4 = [0.702,0.278,0.380]
select surf_pocket4, protein and id [1444,1446,1604,2721,2722,2723,2729,1711,1645,1648,1679,1681,1684,1682,1647,1650,1588,1589,1447,2707] 
set surface_color,  pcol4, surf_pocket4 
set_color pcol5 = [0.902,0.620,0.361]
select surf_pocket5, protein and id [3168,3169,678,703,704,680,2430,677,702,717,2859,2838,2843,2846,2847] 
set surface_color,  pcol5, surf_pocket5 
   
        
        deselect
        
        orient
        