import numpy as np

class pm:

    def __init__(self, name):
        # Parametersettings#
        if name == 'work':

            self.C_on = 17.0  # signalconcentration of activ cel#
            self.K = 18.0  # threshold c#
            self.feedback = 1  # positiv(1) or negative(0) feedback#

            self.radius = 1  # set intial radius of cell
            self.maxrad = 2 * self.radius
            self.maxvolume = 2 * 4 * np.pi * self.radius ** 3  # max r which cell can have before division#

            self.k = 1.0  # factor how close neighboring cells can be 1.0 just touching
            self.g_rate = 0.2  # growthrate of radius

            self.F_g = 0.002  # gravity constant
            self.threeD = 1  # 3d or 2d

            self.thres_skr = 0.3  # threshold for skrinking

            self.cell_age = 0
            self.maxcell_age = 6

            self.gamma_ = 7.0
            self.diff_const = 7.0
            bruch = float(self.diff_const / self.gamma_)
            self.lambda_ = float(np.sqrt(bruch))

        if name == 'non':

            self.C_on = 17.0  # signalconcentration of activ cel#
            self.K = 18.0  # threshold c#
            self.feedback = 1  # positiv(1) or negative(0) feedback#

            self.radius = 1  # set intial radius of cell
            self.maxrad = 2 * self.radius
            self.maxvolume = 2 * 4 * np.pi * self.radius ** 3  # max r which cell can have before division#

            self.k = 1.0  # factor how close neighboring cells can be 1.0 just touching
            self.g_rate = 0.2  # growthrate of radius

            self.F_g = 0.002  # gravity constant
            self.threeD = 1  # 3d or 2d

            self.thres_skr = 0.3  # threshold for skrinking

            self.cell_age = 0
            self.maxcell_age = 6
            self.gamma_ = 7.0
            self.diff_const = 0.1
            bruch = float(self.diff_const / self.gamma_)
            self.lambda_ = float(np.sqrt(bruch))