class MatingGrid:

    # creates a mating grid with two offsets to find closest pairs
    def __init__(self, params):
        self.dims = params['dims']
        self.mating_radius = params['mating_radius']

        # may need to change grid size so that we don't cover too much individuals
        self.grid_size = 2 * self.mating_radius

        self.x = self.dims[0] / self.grid_size
        self.y = self.dims[0] / self.grid_size

        self.offset1 = [[set() for _ in range(int(self.x) + 1)] for j in range(int(self.y) + 1)]
        self.offset2 = [[set() for _ in range(int(self.x) + 1)] for j in range(int(self.y) + 1)]
        self.offset3 = [[set() for _ in range(int(self.x) + 1)] for j in range(int(self.y) + 1)]
        self.offset4 = [[set() for _ in range(int(self.x) + 1)] for j in range(int(self.y) + 1)]

    # adds an individual to the grid
    def add(self, ind):
        self.__add(ind.x, ind.y, ind)

    #removes an individual from the grid
    def remove(self, ind):
        try:
            self.__remove(ind.x, ind.y, ind)
        except KeyError:
            pass

        return None

    #moves an individual to a new position
    def move(self, old_pos, new_pos, ind):
        self.__remove(old_pos[0], old_pos[1], ind)
        self.__add(new_pos[0], new_pos[1], ind)

    #retrieves the individual from the grid
    def get(self, ind):
        """
        :param ind:
        :return set: a set containing all the individuals possible as a mate
        """
        x_pos = ind.x
        y_pos = ind.y
        set1, set2, set3, set4 = self.__get_set(x_pos, y_pos)
        return set1 | set2 | set3 | set4 - {ind}

    #########################
    # Private Methods Below #
    #########################

    def __get_set(self, x_pos, y_pos):
        x_index = int(x_pos // (self.grid_size))
        y_index = int(y_pos // (self.grid_size))
        set1 = self.offset1[int(x_index)][int(y_index)]

        x2_index = int((x_pos + self.grid_size * 0.5) // (self.grid_size))
        y2_index = int((y_pos + self.grid_size * 0.5) // (self.grid_size))
        set2 = self.offset2[int(x2_index)][int(y2_index)]

        x3_index = int((x_pos + self.grid_size * 0.5) // (self.grid_size))
        y3_index = int((y_pos) // (self.grid_size))
        set3 = self.offset2[int(x3_index)][int(y3_index)]

        x4_index = int((x_pos) // (self.grid_size))
        y4_index = int((y_pos + self.grid_size * 0.5) // (self.grid_size))
        set4 = self.offset2[int(x4_index)][int(y4_index)]

        return set1, set2, set3, set4,

    def __add(self, x_pos, y_pos, ind):
        # private method
        set1, set2, set3, set4 = self.__get_set(x_pos, y_pos)
        set1.add(ind)
        set2.add(ind)
        set3.add(ind)
        set4.add(ind)
        return None

    def __remove(self, x_pos, y_pos, ind):
        # private method
        set1, set2, set3, set4 = self.__get_set(x_pos, y_pos)

        set1.remove(ind)
        set2.remove(ind)
        set3.remove(ind)
        set4.remove(ind)
        return None
