#!/usr/bin/env python3
import sys
mm = 0.001
import openems
import openems.geometries
import numpy as np

em = openems.OpenEMS('vert_connector_ms_oshpark', EndCriteria = 1e-5, fmin = 0e6, fmax = 50e9,
                     boundaries = ['PEC', 'PEC', 'PEC', 'PEC', 'PML_12', 'PEC'])
em.fsteps = 1601
copper = openems.Metal(em, 'copper')
pcopper = openems.Metal(em, 'pcopper')
sub1 = openems.Dielectric(em, 'substrate', eps_r=3.2)
sub2 = openems.Dielectric(em, 'substrate', eps_r=4.0)

sub1t = 0.19*mm
sub2t = 1.0*mm

ifoil = 0.0125*mm
ofoil = 0.035*mm
port_length = 0.1*mm
box_length = 2*mm
box_width = 2*mm
ms_width = 0.42*mm
airspace = 1*mm
via_pad = 0.5*mm
via_clearance = 0.275*mm

bt = sub1t + ofoil
bb = -1*(ofoil+sub2t+sub1t)

em.resolution = 25e-6

em.mesh.AddLine('z', sub1t+airspace)
zmin = bb - 1*mm
em.mesh.AddLine('z', zmin)

planar = openems.geometries.planar_full_box(x=[-0.5*box_length, 0.5*box_length],
                                            y=[-0.5*box_width, 0.5*box_width])

clearance_r = via_pad*0.5 + via_clearance
planar.add(sub1, [0, sub1t], priority=1) # sub1 top
planar.add_center_hole(pcopper, [0, ifoil], clearance_r, priority=2) # inner 1 foil
planar.add(sub2, [0, -sub2t], priority=1) # sub2
planar.add(sub1, [-sub2t, -(sub2t+sub1t)], priority=1) # sub1 bottom
planar.add_center_hole(pcopper, [-sub2t, -(sub2t+ifoil)], clearance_r, priority=2) # inner2 foil
planar.add_center_hole(pcopper, [bb, bb+ofoil], 0.75*mm, priority=1) # bottom foil

# ms line
start = np.array([-0.5*box_length+port_length, 0.5*ms_width, sub1t])
stop  = np.array([0, -0.5*ms_width, bt])
copper.AddBox(start, stop, priority=9)

# ms port
start = [-0.5*box_length, ms_width/2.0, sub1t]
stop  = [-0.5*box_length + port_length, ms_width/-2.0, bt]
openems.Port(em, start, stop, direction='x', z=50)

via_z = [[bt,bb],[bt, bt-ofoil], [0, ifoil], [-sub2t, -sub2t-ifoil], [bb+ofoil, bb]]

# ground vias
for n in range(-3,4):
    r = 1 * mm
    c = np.exp(1j*2*np.pi*n*22.0/180.0) * r
    openems.Via(copper, priority=9, x=np.real(c), y=np.imag(c), z=via_z,
                drillradius = 0.25*mm*0.5,
                wall_thickness = 25e-6,
                padradius = via_pad * 0.5,
                padname='2')

# signal via
openems.Via(copper, priority=9, x=0, y=0,
            z=[[bt,bb],[bt, bt-ofoil], [0, ifoil], [-sub2t, -sub2t-ifoil], [bb+ofoil, bb]],
            drillradius = 0.25*mm*0.5,
            wall_thickness = 25e-6,
            padradius = via_pad * 0.5,
            padname='1')

# coax shield
planar.add_center_hole(copper, [bb, zmin], r=1.5*mm/2.0, priority=1)

pin_diameter = 0.695*mm
coax_port_length = 0.2*mm

# pin
start = np.array([0, 0, bb])
stop  = np.array([0, 0, zmin])
copper.AddCylinder(start, stop, 0.5*pin_diameter, priority=9)

# coax goes into Z- PML

command = 'view solve'
if len(sys.argv) > 1:
    command = sys.argv[1]
em.write_kicad(em.name)
em.run_openems(command)
