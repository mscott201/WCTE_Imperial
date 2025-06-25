from Geometry.Device import Device


class TARGET(Device):
    """ A survey target: not a device, but a target location used as a survey reference. No properties.
    Note: these are not to be used for fiducials that are surveyed to measure placements of devices. Fiducials
    are defined in the CLASS that uses them, e.g. CAMERA, SM, etc. """

    # no class properties:
    design_mean = {'T': {}}  # dictionary of mean properties of different kinds of LEDs
    design_scale = {'T': {}}  # scale of variations of properties used to create objects
    design_var = {'T': {}}  # distribution for variations
    design_desc = {}

    def __init__(self, name, container=None, kind='T', place_design={}, place_true={}):
        super().__init__(TARGET, name, container, kind, place_design, place_true)
