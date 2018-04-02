class MatingRaster:
    def __init__(self, params):
        self.dims = params['dims']
        mating_radius = params['mating_radius']

        #
        x = self.dims[0] / (2 * mating_radius)
        y = self.dims[0] / (2 * mating_radius)

        self.offset1 = [[set() for i in range(int(x))] for j in range(int(y))]
        self.offset2 = [[set() for i in range(int(x))] for j in range(int(y))]


    def add(self, ind):
        self.__add(ind.x, ind.y)

    def remove(self, ind):
        self.__remove(ind.x, ind.y)
        return None

    def move(self, old_pos, new_pos):
        self.__remove(old_pos[0], old_pos[1])
        self.__add(new_pos[0], new_pos[1])
        return None

    #########################
    # Private Methods Below #
    #########################

    def __add(self, x_pos, y_pos):
        # private method
        # TODO: NotImplemented
        return None

    def __remove(self, x_pos, y_pos):
        # private method
        # TODO: NotImplemented
        return None
