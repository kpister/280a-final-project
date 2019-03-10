


class Island:
    def __init__(self, x, y):
        self.grid = [[0 for _ in range(x)] for _ in range(y)]
        self.width = x
        self.height = y
        self.min_x = x
        self.max_x = 0
        self.min_y = y
        self.max_y = 0

    def set(self, x, y, val):
        self.grid[x][y] = val
        self.max_x = max(x, self.max_x)
        self.max_y = max(y, self.max_y)
        self.min_x = min(x, self.min_x)
        self.min_y = min(y, self.min_y)

    def get_best(self):
        best_i, best_j, val = -1,-1, 0
        for i, row in enumerate(self.grid):
            for j, cell in enumerate(row):
                if cell <= val:
                    best_i, best_j, val = i, j, cell

        self.max_x = best_i
        self.max_y = best_j
        #print(f'Min: {self.min_x}, {self.min_y}')
        #print(f'Max: {self.max_x}, {self.max_y}')
        #print()
        return best_i, best_j, val

    def conflicts(self, o_i):
        if (self.min_x < o_i.min_x < self.max_x):
            return True
        if (self.min_y < o_i.min_y < self.max_y): 
            return True
        if (o_i.min_x < self.min_x < o_i.max_x):
            return True
        if (o_i.min_y < self.min_y < o_i.max_y):
            return True
        return False

    def list_conflicts(self, ilist):
        for el in ilist:
            if self.conflicts(el):
                return True

        return False
