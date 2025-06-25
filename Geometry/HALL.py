from Geometry.Device import Device
from Geometry.MPMT import MPMT
from Geometry.SM import SM
from Geometry.WCD import WCD

import numpy as np
from scipy.spatial.transform import Rotation


class HALL(Device):
    """A hall contains WCDs, supermodules, or mPMTs. It only has geometry information, no other properties.
    The hall is the top level of the geometry hierarchy, defining the coordinate system for the experiment or survey.
    """

    mpmts_design = {}
    sms_design = {}

    loc_sig = [10.0, 10.0, 10.0]  # mm positioning accuracy
    rot_angles_sig = [0.01, 0.01, 0.01]  # rad rotational angle positioning accuracy

    wcte_wcds = []

    wcte_wcds.append({'name': 'wcte', 'kind': 'WCTE',
                     'loc': [0., 0., 0.],
                     'loc_sig': loc_sig,
                     'rot_axes': 'ZYX',
                     'rot_angles': [0., 0., 0.],
                     'rot_angles_sig': rot_angles_sig})

    wcds_design = {'WCTE':wcte_wcds}

    def __init__(self, name, container=None, kind='WCTE', place_design={}, place_true={}, device_type=WCD):
        super().__init__(HALL, name, container, kind, place_design, place_true)

        # create and place the set of devices
        if device_type == MPMT:
            self.mpmts = self.place_devices(device_type, self.mpmts_design, kind)
            self.sms = None
            self.wcds = None
        elif device_type == SM:
            self.sms = self.place_devices(device_type, self.sms_design, kind)
            # collect all the mpmts in the WCD
            self.mpmts = []
            for sm in self.sms:
                sm.get_mpmts(self.mpmts)
            # give each mpmt a unique name
            for i, mpmt in enumerate(self.mpmts):
                mpmt.name = str(i)
            self.wcds = None
        elif device_type == WCD:
            self.wcds = self.place_devices(device_type, self.wcds_design, kind)
            self.sms = None
