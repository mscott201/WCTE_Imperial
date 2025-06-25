from Geometry.Device import Device


class PMT(Device):
    """ A photomultiplier tube """

    # class properties:
    design_mean = {}  # dictionary of design properties of different kinds of PMTs
    design_scale = {}  # scales of variations of properties used to create objects
    design_var = {}  # distribution of variations (Normal by default, "uniform")

    # default design properties of PMTs:
    # all properties are defined by primitives, so shallow dictionary copy works
    # if new properties are needed, be sure to add to design_desc, def_design_mean
    # and def_design_scale

    design_desc = {'fov': 'radian field of view (cone half angle)',
                   'size': 'mm diameter',
                   'qe': 'quantum efficiency',
                   'gain': 'mV per photoelectron',
                   'delay': 'ns mean delay of signal wrt pulse',
                   'tts': 'ns transit time spread: standard deviation of transit delay'
                   }

    def_design_mean = {'fov': 1.0,  # radian field of view (cone half angle)
                       'size': 75,  # mm diameter
                       'qe': 0.5,  # quantum efficiency
                       'gain': 2.,  # mV per photoelectron
                       'delay': 4.0,  # ns mean delay of signal wrt pulse
                       'tts': 1.0  # ns transit time spread: standard deviation of transit delay
                       }
    def_design_scale = {'fov': 0.001,
                        'size': 0.1,
                        'qe': 0.01,
                        'gain': 0.1,
                        'delay': 0.5,
                        'tts': 0.1
                        }
    def_design_var = {'gain': 'gamma',
                      'delay': 'gamma',
                      'tts': 'gamma'}

    # Standard 3inch PMT:
    p3_design_mean = def_design_mean.copy()
    p3_design_mean['gain'] = 2.1  # 2.1 mV per photoelectron
    p3_design_scale = def_design_scale.copy()
    p3_design_var = def_design_var.copy()

    design_mean['P3'] = p3_design_mean
    design_scale['P3'] = p3_design_scale
    design_var['P3'] = p3_design_var

    # Modified 3inch PMT:
    p31_design_mean = def_design_mean.copy()
    p31_design_mean['gain'] = 2.5  # 2.5 mV per photoelectron
    p31_design_scale = def_design_scale.copy()
    p31_design_var = def_design_var.copy()

    design_mean['P31'] = p31_design_mean
    design_scale['P31'] = p31_design_scale
    design_var['P31'] = p31_design_var

    def get_xy_points(self, place_info, device_for_coordinate_system=None):
        """Return set of points that shows extent on x-y plane (z=0)"""
        return self.get_circle_points(20, place_info, device_for_coordinate_system)

    def __init__(self, name, container=None, kind='P3', place_design={}, place_true={}):
        super().__init__(PMT, name, container, kind, place_design, place_true)
