


class Island:
    def __init__(self, x, y):
        self.grid = [[0 for _ in range(x)] for _ in range(y)]
        self.width = x
        self.height = y
        self.min_x = x
        self.max_x = 0
        self.min_y = y
        self.max_y = 0
        self.best = 0

    def set(self, x, y, val):
        self.grid[x][y] = val
        self.min_x = min(x, self.min_x)
        self.min_y = min(y, self.min_y)

        if val <= self.best:
            self.best = val
            self.max_x = x
            self.max_y = y

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

    def explore(self, x, y, mem):
        if mem.get(x, y):
            self.set(x, y, mem.clear(x, y))
            return True
        return False

def CreateIslands(mem):
    islands = []
    for i, row in enumerate(mem.table):
        for j, cell in enumerate(row):
            if cell != 0:
                island = Island(mem.width, mem.height)
                island.set(i, j, mem.clear(i, j))
                queue = [(i, j)]
                while queue:
                    x, y = queue.pop()
                    if island.explore(x-1, y, mem):
                        queue.insert(0, (x-1, y))
                    if island.explore(x+1, y, mem):
                        queue.insert(0, (x+1, y))
                    if island.explore(x, y-1, mem):
                        queue.insert(0, (x, y-1))
                    if island.explore(x, y+1, mem):
                        queue.insert(0, (x, y+1))
                    if island.explore(x+1, y+1, mem):
                        queue.insert(0, (x+1, y+1))

                if island.best < -1:
                    islands.append(island)
    return islands
