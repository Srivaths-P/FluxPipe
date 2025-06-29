#!/usr/bin/env python3
import sys
from openems import OpenEMS, Box, Cylinder, Via, Port, Metal, Dielectric, geometries
import numpy as np
em = OpenEMS('via_1_3_oshpark', EndCriteria = 1e-4, fmin = 0, fmax = 40e9,
             boundaries = ['PEC', 'PEC', 'PEC', 'PEC', 'PEC', 'PEC'],
             fsteps=801)
copper = Metal(em, 'copper')
pcopper = Metal(em, 'pcopper')
sub1 = Dielectric(em, 'substrate', eps_r=3.2)
sub2 = Dielectric(em, 'substrate', eps_r=4)

sub1t = 0.166e-3
sub2t = 47*25.4e-6
ifoil = 0.0125e-3
ofoil = 0.035e-3
port_length = 0.1e-3
box_length = 5e-3
box_width = 2e-3
sl_width = 0.24e-3
ms_width = 0.35e-3
airspace = 1e-3

bt = sub1t + ofoil
bb = -1*(ofoil+sub2t+sub1t)

em.resolution = 50e-6

planar = geometries.planar_full_box(
    x=[-0.5*box_length, 0.5*box_length],
    y=[-0.5*box_width, 0.5*box_width])

clearance_r = 0.86e-3 * 0.5

em.mesh.AddLine('z', sub1t+airspace)
em.mesh.AddLine('z', -1.0*(2*sub1t+sub2t+airspace))

planar.add(sub1, [0, sub1t]) # sub1 top
planar.add_center_hole(pcopper, [0, ifoil], clearance_r, priority=2) # inner 1 foil
planar.add(sub2, [0, -sub2t]) # core
planar.add(sub1, [-sub2t, -(sub2t+sub1t)]) # sub1 bot
planar.add_center_hole(pcopper, [bb, -(sub2t+sub1t)], clearance_r, priority=2) # bottom foil

# line (sl)
start = np.array([0, 0.5*sl_width, -sub2t])
stop  = np.array([0.5*box_length-port_length, -0.5*sl_width, -sub2t-ifoil])
Box(copper, 9, start, stop)

# line (ms)
Box(copper, 9, [-0.5*box_length+port_length, 0.5*ms_width, sub1t], [0, -0.5*ms_width, bt])

# port (ms)
start = [-0.5*box_length, ms_width/2.0, sub1t]
stop  = [-0.5*box_length + port_length, ms_width/-2.0, bt]
Port(em, start, stop, direction='x', z=50)

# port (sl)
start = [0.5*box_length, sl_width/2.0, -sub2t]
stop  = [0.5*box_length - port_length, sl_width/-2.0, -sub2t-ifoil]
Port(em, start, stop, direction='x', z=50)

Via(copper, priority=9, x=0, y=0,
            z=[[bb, bt], [bt, bt-ofoil], [0, ifoil], [-sub2t, -sub2t-ifoil], [bb+ofoil, bb]],
            drillradius = 0.5*0.25e-3,
            wall_thickness = 25e-6,
            padradius = 0.46e-3*0.5,
            padname = '1')

for x in range(-3,4):
    x *= 0.5e-3
    for y in [-0.75e-3, 0.75e-3]:
        Cylinder(copper, 9, [x, y, bb], [x, y, bt], 0.3e-3*0.5)
        Cylinder(copper, 9, [x, y, bt], [x, y, bt-ofoil], 0.46e-3*0.5)

command = 'view solve'
if len(sys.argv) > 1:
    command = sys.argv[1]
em.run_openems(command)
