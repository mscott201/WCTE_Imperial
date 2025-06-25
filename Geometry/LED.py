from Geometry.Device import Device


class LED(Device):
    """ A LED-fibre light source """

    # class properties:
    design_mean = {}  # dictionary of mean properties of different kinds of LEDs
    design_scale = {}  # scale of variations of properties used to create objects
    design_var = {}  # distribution for variations

    # default properties of LEDs
    # all properties are defined by primitives, so shallow dictionary copy works
    # if new properties are needed, be sure to add to def_design_mean, design_desc,
    # and def_design_scale

    design_desc = {'cone_angle': 'radian cone angle',
                   'cone_ring_width': 'if zero, then full cone, otherwise specifies ring width (radians)',
                   'intensity': 'billion photons per ns per volt',
                   'intensity_rel_sig': 'relative standard deviation',
                   'delay': 'ns delay for leading edge',
                   'delay_jitter': 'ns jitter',
                   'rise_time': 'ns rise time',
                   'pulse_width': 'ns mean pulse width (flat top)',
                   'pulse_width_jitter': 'standard deviation of above',
                   'wavelength': 'nm peak wavelength'
                   }

    def_design_mean = {'cone_angle': 1.0,  # radian cone angle
                       'cone_ring_width': 0.0,  # if zero, then full cone, otherwise specifies ring width (radians)
                       'intensity': 1.E8,  # photons per ns per volt
                       'intensity_rel_sig': 0.01,  # relative standard deviation
                       'delay': 4.,  # ns delay for leading edge
                       'delay_jitter': 0.05,  # ns jitter
                       'rise_time': 0.5,  # ns rise time
                       'pulse_width': 1.0,  # ns mean pulse width (flat top)
                       'pulse_width_jitter': 0.05,  # standard deviation of above
                       'wavelength': 470.  # nm peak wavelength
                       }
    def_design_scale = {'cone_angle': 0.01,
                        'cone_ring_width': 0.0,
                        'intensity': 2.E7,
                        'intensity_rel_sig': 0.,
                        'delay': 1.0,
                        'delay_jitter': 0.01,
                        'rise_time': 0.1,
                        'pulse_width': 0.1,
                        'pulse_width_jitter': 0.01,
                        'wavelength': 2.
                        }
    def_design_var = {'intensity': 'gamma',
                      'delay': 'gamma',
                      'rise_time': 'gamma',
                      'pulse_width_jitter': 'gamma'}

    # Diffuse LED:
    ld_design_mean = def_design_mean.copy()
    ld_design_scale = def_design_scale.copy()
    ld_design_var = def_design_var.copy()

    design_mean['LD'] = ld_design_mean
    design_scale['LD'] = ld_design_scale
    design_var['LD'] = ld_design_var

    # Collimated LEDs:
    lc_design_mean = def_design_mean.copy()
    lc_design_mean['cone_angle'] = 0.1745  # 10 degrees
    lc_design_mean['intensity'] = 1.E7  # 0.01 billion photons per ns per volt
    lc_design_scale = def_design_scale.copy()
    lc_design_scale['intensity'] = 2.E6
    lc_design_var = def_design_var.copy()

    design_mean['LC'] = lc_design_mean
    design_scale['LC'] = lc_design_scale
    design_var['LC'] = lc_design_var

    lc15_design_mean = lc_design_mean.copy()
    lc15_design_mean['cone_angle'] = 0.2618  # 15 degrees
    design_mean['LC15'] = lc15_design_mean
    design_scale['LC15'] = lc_design_scale
    design_var['LC15'] = lc_design_var

    lc30_design_mean = lc_design_mean.copy()
    lc30_design_mean['cone_angle'] = 0.5236  # 30 degrees
    design_mean['LC30'] = lc30_design_mean
    design_scale['LC30'] = lc_design_scale
    design_var['LC30'] = lc_design_var

    def __init__(self, name, container=None, kind='LD', place_design={}, place_true={}):
        super().__init__(LED, name, container, kind, place_design, place_true)
