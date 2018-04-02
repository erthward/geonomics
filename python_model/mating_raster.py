class MatingRaster:
    def __init__(self, params):
        self.dims = params['dims']
        mating_radius = params['mating_radius']

        x = dims[0] / (2 * mating_radius)
        y = dims[0] / (2 * mating_radius)

        self.offset1 = [[set() for i in range(x)] for j in range(y)]
        self.offset2 = [[set() for i in range(x)] for j in range(y)]

    def add(self, ind):

        return None

    def remove(self, ind):

        return None

    def move(self, ind):

        return None

