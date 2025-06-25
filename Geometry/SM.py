from Geometry.Device import Device
from Geometry.MPMT import MPMT
from Geometry.CAMERA import CAMERA
from Geometry.TARGET import TARGET

import numpy as np
from scipy.spatial.transform import Rotation


class SM(Device):
    """A super module consists of a set of either mPMTs or super modules.
    It can also contain cameras and targets. It only has geometry information, no other properties.
    """

    # A dictionary of device kinds and placements in the super module:
    devices_design = {}
    # A dictionary of camera kinds and placements in the super module:
    cameras_design = {}
    # A dictionary of target kinds and placements in the super module:
    targets_design = {}

    ssm_mpmts = []
    # 3 x 2 rectangular pattern of MPMTs (for testing):
    for i in range(-1, 2):
        for j in range(-1, 1):
            ssm_mpmts.append({'name': str(len(ssm_mpmts)), 'kind': 'MR',
                              'loc': [600. * i, 600. * (j + 0.5), -100., ],
                              'loc_sig': [1.0, 1.0, 1.0],
                              'rot_axes': 'XZ',
                              'rot_angles': [0., 0.],
                              'rot_angles_sig': [0.01, 0.01]})

    devices_design['SSM'] = ssm_mpmts

    ssm_cameras = []
    # 2 x 2 rectangular pattern of CAMERAs (for testing):
    for i in np.arange(-0.5, 0.51, 1.0):
        for j in np.arange(-0.5, 0.51, 1.0):
            ssm_cameras.append({'name': str(len(ssm_cameras)), 'kind': 'C',
                                'loc': [600. * i, 600. * (j + 0.5), -100., ],
                                'loc_sig': [1.0, 1.0, 1.0],
                                'rot_axes': 'XZ',
                                'rot_angles': [0., 0.],
                                'rot_angles_sig': [0.01, 0.01]})

    cameras_design['SSM'] = ssm_cameras

    # Super modules that make up the WCTE WCD:
    # The z-axis for each super module is pointing upwards during each super module construction and their initial
    # surveys. Once assembled the WCTE z-axis is aligned with the beam direction

    tb_pitch = 580.  # mm separation of mPMT centres on top and bottom (x and y the same)
    wcte_diameter = 3422.166  # mm separation of mPMT baseplates on opposite barrel locations
    barrel_vertical_pitch = 580.  # mm separation of mPMT centres for rows

    loc_sig = [1.0, 1.0, 1.0]  # mm positioning accuracy
    rot_angle_sig = 0.01  # radians

    # Bottom super module:
    ######################
    # Origin and z-axis coincides with that of centre mPMT origin.
    # Ordering of MPTs increases with phi. Second mPMT is displaced along SM +x axis.

    bottom_mpmts = []

    # Start with mPMT at centre of bottom
    offsets = [[0., 0., 0.]]

    offs = [-tb_pitch, 0, tb_pitch]
    for i in range(8):
        offset = [offs[[2, 2, 1, 0, 0, 0, 1, 2][i]], offs[[1, 2, 2, 2, 1, 0, 0, 0][i]], 0.]
        offsets.append(offset)

    offs = [-2. * tb_pitch, -tb_pitch, 0, tb_pitch, 2. * tb_pitch]
    for i in range(12):
        offset = [offs[[4, 4, 3, 2, 1, 0, 0, 0, 1, 2, 3, 4][i]], offs[[2, 3, 4, 4, 4, 3, 2, 1, 0, 0, 0, 1][i]], 0.]
        offsets.append(offset)

    # rotate mPMTs to put feed-throughs in correct orientations... pi multiplier starting at #0
    y_rots = [0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0]
    y_rots = [0.5, 0.5, 0.5, 1.5, 1.5, 1.5, 1.5, 1.5, 0.5, 0.5, 0.5, 0.5, 0.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 0.5, 0.5]

    mpmt_kinds = ['ME']*12 + ['FD','ME','FD','ME','FD','ME','FD','ME','ME']

    for i, (offset, y_rot, mpmt_kind) in enumerate(zip(offsets, y_rots, mpmt_kinds)):
        location = offset.copy()
        if i in [12, 14, 16, 18]:
            location[2] = 7.775  # mm offset for FD mPMTs
        bottom_mpmts.append({
            'kind': mpmt_kind,
            'loc': location,
            'loc_sig': loc_sig,
            'rot_axes': 'ZYX',
            'rot_angles': [(y_rot + 0.5) * np.pi, 0., 0.],
            'rot_angles_sig': [rot_angle_sig] * 3
        })

    devices_design['bottom'] = bottom_mpmts

    bottom_cameras = []
    # located at 4 corners of bottom super module
    # they are movable: for now, a typical location is given (same for all)
    camera_radius = 990 * np.sqrt(2)  # typically 990 mm transverse from centre
    camera_z_bot = 200.  # typically 200 mm
    camera_angle = 0.92729  # radians (sin(0.92729) = 0.8)

    for i_cam in range(4):
        # start with a camera located on the bottom endcap x axis
        loc = [camera_radius, 0., camera_z_bot]
        # now rotate it around the bottom endcap z axis
        phi_angle = np.pi / 4. + 2. * np.pi * i_cam / 4
        rot_phi = Rotation.from_euler('Z', phi_angle)
        rot_loc = rot_phi.apply(loc)
        # rotations of the normal defined by 3 rotations
        rot_angles = [phi_angle + np.pi / 2., 0., -camera_angle]

        bottom_cameras.append({'name': str(i_cam), 'kind': 'C',
                               'loc': rot_loc,
                               'loc_sig': [1.0, 1.0, 1.0],
                               'rot_axes': 'ZYX',
                               'rot_angles': rot_angles,
                               'rot_angles_sig': [rot_angle_sig] * 3})

    cameras_design['bottom'] = bottom_cameras

    bottom_targets = []
    # located at a point extending from the support beams of bottom super module
    # they are not to be used to define the coordinate system of the super module
    # instead, they serve as references for determining the placement of the super module after full assembly
    z_inner_top = -20.  # mm of the inner top surface of the support beams
    height_inner_beam = 110.  # mm height of the inner volume of the support beams
    dist_mpmt_to_beam_centre = 267.5 + 45. / 2.  # mm distance from an mPMT to the centre of the support beam
    width_outer_beam = 45.  # mm width of the outer surface of the support beams
    length_beams = [1920., 3030., 3570.]  # lengths of the small, medium, and long beams
    target_extension = 48  # mm extension of the target from the end of the support beams (est by minimizing rms)

    target_z = z_inner_top - height_inner_beam / 2.
    target_xs = [length_beams[2] / 2. + target_extension,
                 length_beams[1] / 2. + target_extension,
                 length_beams[0] / 2. + target_extension,
                 dist_mpmt_to_beam_centre,
                 -dist_mpmt_to_beam_centre,
                 -length_beams[0] / 2. - target_extension,
                 -length_beams[1] / 2. - target_extension,
                 -length_beams[2] / 2. - target_extension,
                 -length_beams[2] / 2. - target_extension,
                 -length_beams[1] / 2. - target_extension,
                 -length_beams[0] / 2. - target_extension,
                 -dist_mpmt_to_beam_centre,
                 dist_mpmt_to_beam_centre,
                 length_beams[0] / 2. + target_extension,
                 length_beams[1] / 2. + target_extension,
                 length_beams[2] / 2. + target_extension]
    target_ys = [dist_mpmt_to_beam_centre,
                 1. * tb_pitch + dist_mpmt_to_beam_centre,
                 2. * tb_pitch + dist_mpmt_to_beam_centre,
                 length_beams[2] / 2. + target_extension,
                 length_beams[2] / 2. + target_extension,
                 2. * tb_pitch + dist_mpmt_to_beam_centre,
                 1. * tb_pitch + dist_mpmt_to_beam_centre,
                 dist_mpmt_to_beam_centre,
                 -dist_mpmt_to_beam_centre,
                 -1. * tb_pitch - dist_mpmt_to_beam_centre,
                 -2. * tb_pitch - dist_mpmt_to_beam_centre,
                 -length_beams[2] / 2. - target_extension,
                 -length_beams[2] / 2. - target_extension,
                 -2. * tb_pitch - dist_mpmt_to_beam_centre,
                 -1. * tb_pitch - dist_mpmt_to_beam_centre,
                 -dist_mpmt_to_beam_centre]

    for i_target in range(16):
        loc = [target_xs[i_target], target_ys[i_target], target_z]

        bottom_targets.append({'name': str(i_target), 'kind': 'T',
                               'loc': loc,
                               'loc_sig': [1.0, 1.0, 1.0],
                               'rot_axes': 'Z',
                               'rot_angles': 0.,
                               'rot_angles_sig': 0.})

    targets_design['bottom'] = bottom_targets

    # Top super module:
    ###################
    # Constructed with the mPMTs z-axes pointing downwards, unlike the bottom super module.
    # Second mPMT is displaced along SM +x axis.

    top_mpmts = []

    # as built:
    y_rots = [0, 0.5, 0.5, 1.5, 1.5, 1.5, 0., 1.5, 0.5, 0.5, 0.5, 0.5, 1.5, 1.5, 0., 1.5, 1.5, 1.5, 1.5, 0.5, 0.5]

    for offset, y_rot in zip(offsets, y_rots):
        location = offset
        top_mpmts.append({
            'kind': 'ME',
            'loc': location,
            'loc_sig': loc_sig,
            'rot_axes': 'ZYX',
            'rot_angles': [(y_rot - 0.5) * np.pi, np.pi, 0.],
            'rot_angles_sig': [rot_angle_sig] * 3
        })

    devices_design['top'] = top_mpmts

    top_cameras = []
    # located at 4 corners of top super module
    camera_z_top = 270.  # typically 270 mm

    for i_cam in range(4):
        # start with a camera located on the bottom endcap x axis
        loc = [camera_radius, 0., -camera_z_top]
        # now rotate it around the bottom endcap z axis
        phi_angle = np.pi / 4. + 2. * np.pi * i_cam / 4
        rot_phi = Rotation.from_euler('Z', phi_angle)
        rot_loc = rot_phi.apply(loc)
        # rotations of the normal defined by 2 extrinsic rotations
        rot_angles = [phi_angle + np.pi / 2., 0., -np.pi / 2. - camera_angle]

        top_cameras.append({'name': str(i_cam + 4), 'kind': 'C',
                            'loc': rot_loc,
                            'loc_sig': [1.0, 1.0, 1.0],
                            'rot_axes': 'ZYX',
                            'rot_angles': rot_angles,
                            'rot_angles_sig': [rot_angle_sig] * 3})

    cameras_design['top'] = top_cameras

    top_targets = []
    # located at a point extending from the support beams of bottom super module
    # they are not to be used to define the coordinate system of the super module
    # instead, they serve as references for determining the placement of the super module after full assembly
    z_inner_top = -10.  # mm of the inner top surface of the support beams (different from bottom)

    target_z = z_inner_top - height_inner_beam / 2.

    for i_target in range(16):
        loc = [target_xs[i_target], target_ys[i_target], target_z]

        top_targets.append({'name': str(i_target), 'kind': 'T',
                            'loc': loc,
                            'loc_sig': [1.0, 1.0, 1.0],
                            'rot_axes': 'Z',
                            'rot_angles': 0.,
                            'rot_angles_sig': 0.})

    targets_design['top'] = top_targets

    # Barrel super module:
    ####################
    # Origin on cylinder axis, z-axis along that axis, z=0 at centre of second from bottom row of
    # mPMTs. Ordering of mPMTs increases with phi. Second mPMT is displaced along SM +x axis.

    barrel_mpmts = []

    n_col = 16
    for i_row in range(-1, 3):
        loc = [wcte_diameter / 2., 0., i_row * barrel_vertical_pitch]
        for j_col in range(n_col):
            phi_angle = 2. * np.pi * j_col / n_col
            rot_phi = Rotation.from_euler('Z', phi_angle)
            rot_loc = rot_phi.apply(loc)
            # rotations of the normal defined by 3 extrinsic rotations
            rot_angles = [np.pi, np.pi / 2., -phi_angle]
            barrel_mpmts.append({'kind': 'ME',
                                 'loc': rot_loc,
                                 'loc_sig': loc_sig,
                                 'rot_axes': 'ZYX',
                                 'rot_angles': rot_angles,
                                 'rot_angles_sig': [rot_angle_sig] * 3})

    devices_design['barrel'] = barrel_mpmts

    def __init__(self, name, container=None, kind='SSM', place_design={}, place_true={}, device_type=MPMT):
        super().__init__(SM, name, container, kind, place_design, place_true)

        # create and place the set of devices
        if device_type == MPMT:
            self.mpmts = self.place_devices(device_type, self.devices_design, kind)
            self.sms = None
        elif device_type == SM:
            self.sms = self.place_devices(device_type, self.devices_design, kind)
            # collect all the mpmts in the WCD
            self.mpmts = []
            for sm in self.sms:
                sm.get_mpmts(self.mpmts)
            # give each mpmt a unique name
            for i, mpmt in enumerate(self.mpmts):
                mpmt.name = str(i)

        # create and place the set of cameras
        self.cameras = None
        if self.cameras_design.get(kind, None) is not None:
            self.cameras = self.place_devices(CAMERA, self.cameras_design, kind)

        # create and place the set of targets
        self.targets = None
        if self.targets_design.get(kind, None) is not None:
            self.targets = self.place_devices(TARGET, self.targets_design, kind)

    def get_mpmts(self, mpmt_list):
        # recursive discovery of mpmts in super modules
        if self.mpmts is None:
            for sm in self.sms:
                sm.get_mpmts(mpmt_list)
        else:
            mpmt_list.extend(self.mpmts)

    def get_cameras(self, camera_list):
        # recursive discovery of cameras in super modules
        if self.sms is not None:
            for sm in self.sms:
                sm.get_cameras(camera_list)
        else:
            if self.cameras is not None:
                camera_list.extend(self.cameras)

    def get_targets(self, target_list):
        # recursive discovery of targets in super modules
        if self.sms is not None:
            for sm in self.sms:
                sm.get_targets(target_list)
        else:
            if self.targets is not None:
                target_list.extend(self.targets)
