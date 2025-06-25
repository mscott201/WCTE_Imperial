from Geometry.Device import Device
import numpy as np

class CAMERA(Device):
    """ A camera system"""

    # class properties:
    design_mean = {}  # dictionary of design properties of different kinds of Cameras
    design_scale = {}  # scales of variations of properties used to create objects
    design_var = {}  # distribution of variations (Normal by default, "uniform")

    # default design properties of Cameras:
    # all properties are defined by primitives, so shallow dictionary copy works
    # if new properties are needed, be sure to add to design_desc, def_design_mean
    # and def_design_scale

    design_desc = {'fov': 'radian field of view (cone half angle)',
                   'size': 'mm diameter of dome of camera housing at front surface',
                   'housing_size': 'mm diameter of camera housing at front surface',
                   'housing_length': 'mm length of camera housing from front surface to back surface',
                   }

    def_design_mean = {'fov': 1.0,  # radian field of view (cone half angle)
                       'size': 170,  # mm diameter
                       'housing_size': 254,  # mm diameter
                       'housing_length': 162.275  # mm length
                       }
    def_design_scale = {'fov': 0.001,
                        'size': 0.1,
                        'housing_size': 0.1,
                        'housing_length': 0.1
                        }
    def_design_var = {}

    # Standard Camera:
    c_design_mean = def_design_mean.copy()
    c_design_scale = def_design_scale.copy()
    c_design_var = def_design_var.copy()

    design_mean['C'] = c_design_mean
    design_scale['C'] = c_design_scale
    design_var['C'] = c_design_var

    # Survey holes in camera housing:
    # 16 holes, some are available for surveying the camera housing
    # The holes are located at diameter 236.22 mm
    # The holes are located at the front surface of the camera housing
    n_holes = 16
    survey_radius = 236.22 / 2
    survey_holes_xy_points = []
    for i in range(n_holes):
        angle = 2 * i * np.pi / n_holes
        survey_holes_xy_points.append([survey_radius * np.cos(angle), survey_radius * np.sin(angle)])

    # Fiducial points (centres of corner cube reflectors) are offset in zm by some distance: set at 0 mm to show holes
    fiducial_z_offset = 60.  # mm
    fiducials = []
    for i in range(n_holes):
        fiducials.append([survey_holes_xy_points[i][0], survey_holes_xy_points[i][1], fiducial_z_offset])

    def get_xy_points(self, place_info, feature='dome', device_for_coordinate_system=None):
        """Return set of points that shows features
        To show outer housing front surface, set feature='front_housing'
        To show outer housing back surface, set feature='back_housing'
        """

        if feature == 'dome':
            return self.get_circle_points(20, place_info, device_for_coordinate_system)
        else:
            xy_points = []
            diameter = self.get_properties('design')['housing_size']
            n_points = 20
            if feature == 'front_housing':
                for i in range(n_points):
                    angle = 2 * i * np.pi / n_points
                    xy_points.append([diameter / 2 * np.cos(angle), diameter / 2 * np.sin(angle), 0])
            elif feature == 'back_housing':
                z_offset = -1.*self.get_properties('design')['housing_length']
                for i in range(n_points):
                    angle = 2 * i * np.pi / n_points
                    xy_points.append([diameter / 2 * np.cos(angle), diameter / 2 * np.sin(angle), z_offset])

            return self.get_transformed_points(xy_points, place_info, device_for_coordinate_system)

    def get_fiducials(self, place_info, device_for_coordinate_system=None):
        """Return the set of fiducial points for surveying (locations of the corner cube reflectors)"""
        return self.get_transformed_points(self.fiducials, place_info, device_for_coordinate_system)

    def __init__(self, name, container=None, kind='C', place_design={}, place_true={}):
        super().__init__(CAMERA, name, container, kind, place_design, place_true)
