class MatingRaster:
    def __init__(self, params, population):
        self.dims = params['dims']
        mating_radius = params['mating_radius']

        x = dims[0] / (2 * mating_radius)
        y = dims[0] / (2 * mating_radius)

        self.offset1 = [[set() for i in range(x)] for j in range(y)]
        self.offset2 = [[set() for i in range(x)] for j in range(y)]

    def add(self, ind):

        x = round(ind.x, mating_radius - int(mating_radius))
        y = round(ind.y, mating_radius - int(mating_radius))
        offset1[int(x)][int(y)] = ind

        x2 = round(ind.x, mating_radius - int(0.5 * mating_radius))
        y2 = round(ind.y, mating_radius - int(0.5 * mating_radius))
        offset2[int(x)][int(y)] = ind

    def remove(self, ind):

        return None

    def move(self, ind):

        return None




