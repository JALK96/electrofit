#!/usr/bin/env python
# coding: utf-8
import mdtraj as md
import os

home = os.getcwd()
sim_dir = "/home/lf1071fu/project_b3/ip6_cation/gaff/ip6_IPL/param_2/sample_confs"

raw_traj = md.load(f"{ sim_dir }/fitted_traj.xtc", top=f"{ sim_dir }/md.gro")
ipl = raw_traj.top.select("resname MOL")
traj = raw_traj.atom_slice(ipl)

# each frame is 100 ps
configs = [c for c, t in zip(traj, traj.time) if t % 1000 == 0]

for i, c in enumerate(configs):
    name = "config_" + str(i) 
    if not os.path.isdir(name):
        os.mkdir(name)
    file = name + "/conform.pdb"
    c.save_pdb(file)

os.chdir(home)