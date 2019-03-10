

class Table:
    def __init__(self, width, height, default):
        self.table = [[default for _ in range(width)] for _ in range(height)]
        self.width = width
        self.height = height
        self.default = default

    def get(self, i, j):
        if i < 0 or j < 0:
            return self.default
        if i >= self.height or j >= self.width:
            return self.default

        return self.table[i][j]

    def set(self, i, j, val):
        if i < 0 or j < 0 or i >= self.height or j >= self.width:
            return 'err'
        self.table[i][j] = val

    def clear(self, i, j):
        if i < 0 or j < 0 or i >= self.height or j >= self.width:
            return 'err'
        val = self.table[i][j]
        self.table[i][j] = self.default
        return val
        
