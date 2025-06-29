import numpy as np
import openems

class Polygon(openems.Object):
    def __init__(self,
                 material,
                 points,
                 elevation,
                 priority=1,
                 normal_direction = 'z',
                 pcb_layer = 'F.Cu',
                 mirror = '',
                 **kwargs
    ):
        """
        Add a polygon, generally assumed to be in xy plane, but can be changed to yz or xz
        name: any unique name
        material: the name of a previously defined material
        priority: in the case of overlapping materials, the one with higher priority will be used
        points: pairs of points (xy, yz or xz) [[x1, y1], [x2, y2], ...]
        elevation: start and stop points in the normal direction
        normal_direction: optional, default = z, direction normal to the polygon - 'x', 'y' or 'z'
        """
        self.points = np.array(points)
        if mirror != '':
            assert(normal_direction == 'z')
        if 'x' in mirror:
            self.points[:,0] *= -1.0
        if 'y' in mirror:
            self.points[:,1] *= -1.0
        self.priority = priority
        self.material = material
        self.elevation = elevation
        self.pcb_layer = pcb_layer
        self.normal_direction = normal_direction
        self.em = material.em
        name = self.em.get_name(None)
        self.em.objects[name] = self
        self.kwargs = kwargs

    def generate_kicad(self, g):
        if self.material.__class__.__name__ == 'Dielectric':
            return
        if self.pcb_layer == None:
            return
        if 'is_custom_pad' in self.kwargs.keys():
            name = self.kwargs["pad_name"]
            x = self.kwargs["x"]*1000.0
            y = self.kwargs["y"]*1000.0
            g.add_custom_pad(name, x, y, [1000.0*self.points], layer=self.pcb_layer)
        else:
            g.add_polygon(points = 1000.0 * self.points, layer = self.pcb_layer, width=0)

    def generate_octave(self):
        height = self.elevation[1] - self.elevation[0]
        self.material.material.AddLinPoly(np.swapaxes(self.points, 0, 1),
                                          self.normal_direction,
                                          self.elevation[0],
                                          height,
                                          priority=self.priority)
