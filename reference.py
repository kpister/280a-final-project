"""Reference
Usage: 
    reference.py [options]

Options:
    -h, --help          print this message and exit
    -n, --cnum INT      set the chromosome number to test [default: 10]
    -f, --cmap FILE     set the cmap file to read [default: hg19_BspQI.cmap]

Authors: kaiser, steven
"""
from docopt import docopt 
import random 
import numpy as np
from island import Island

class Reference:
    def __init__(self, filename, cnumber):
        self.cnumber = cnumber
        # read file, build positions and distances
        self.reference = self.parse_cmap(filename, cnumber)
        self.distances = self.distances[:50]

    def parse_cmap(self, filename, cnumber):
        # cnumber = the chromosome we are targeting
        # read the postions of the cuts from restriction enzyme
        self.positions = [0]
        self.distances = []

        with open(filename) as cmap:
            # the cmap file has the following tab separated headers
            #h CMapId	ContigLength	NumSites	SiteID	LabelChannel	Position	StdDev	Coverage	Occurrence
            for row in cmap:
                # '#' defines a comment in cmap
                if row[0] == '#':
                    continue

                values = row.split('\t')
                chrom = int(values[0])
                if chrom == cnumber and int(values[-1]) >= 1:
                    self.clen = int(float(values[1]))
                    self.num_sites = int(values[2])

                    # store distances in array
                    position = int(float(values[5]))
                    self.positions.append(position)

                    # distance between last cut and current cut
                    self.distances.append(position - self.positions[-2])

                # save time and exit early
                elif chrom > cnumber:
                    break


    # create an optical read based on reference
    def generate_read(self):
        # randomized on length, and structural variations present
        min_length = 5
        max_length = 10
        min_delete = 20
        max_delete = 30
        expected_merges = 1 # np.random.poisson() # TODO Find a good value
        expected_splits = 1 # np.random.poisson() # TODO Find a good value
        expected_inversions = 1
        long_delete_prob = 1
        long_delete_length = 0
        long_delete_start = 0

        length = random.randint(min_length, max_length)

        if random.random() < long_delete_prob:
            long_delete_length = random.randint(min_delete, max_delete)
            long_delete_start = random.randint(length//4, 3*length//4)

        start_pos = random.randint(0, len(self.distances)-length - long_delete_length)

        read = Read(start_pos)
        i = start_pos
        while i < start_pos+length+long_delete_length:
            if i > start_pos + long_delete_start and i < start_pos + long_delete_start + long_delete_length:
                i += 1
                continue

            if random.random() > 1.0 / expected_merges and i + 1 < start_pos + length:
                read.cuts.append(self.distances[i] + self.distances[i+1])
                expected_merges -= 1
                i += 2
            elif random.random() > 1.0 / expected_splits:
                split = self.distances[i] - random.randint(0, self.distances[i])
                read.cuts.append(split)
                read.cuts.append(self.distances[i] - split)
                expected_splits -= 1
                i += 1
            elif random.random() > 1.0 / expected_inversions and i + 1 < start_pos + length:
                read.cuts.append(self.distances[i+1])
                read.cuts.append(self.distances[i])
                expected_inversions -= 1
                i += 2
            else:
                read.cuts.append(self.distances[i])
                i += 1

        # return Read object
        return read

    # C(s_i, r_j) = min(C(s_{i-1}, r_{j-1}) + d_{ij}, C(s_{i-1}, r_j) + d_gap, C(s_i, r_{j-1}) + d_gap)
    # With memoization (using a table or something of the sort): O(len(s) * len(r))
    # also space of O(len(s)*len(r)). With d+c, we can reduce to O(len(s) + len(r))
    def locate_read(self, read):
        # find where in reference the read fits
        # local alignment means no penalty for shifting
        table = [[0 for _ in read.cuts] for _ in self.distances]
        path = [["" for _ in read.cuts] for _ in self.distances]
        ref = self.distances

        def dist(x, y):
            if x == y:
                return -1
            return 1

        for i, ref_i in enumerate(ref):
            for j, read_j in enumerate(read.cuts):
                opt1 = opt2 = opt3 = opt4 = opt5 = opt6 = 0
                left = right = diag = diag2 = i2j1 = i1j2 = 0
                if j > 0:
                    left = table[i][j-1]
                if i > 0:
                    right = table[i-1][j]
                if i > 0 and j > 0:
                    diag = table[i-1][j-1]
                if i > 1 and j > 1:
                    diag2 = table[i-2][j-2]
                if i > 1 and j > 0:
                    i2j1 = table[i-2][j-1]
                if i > 0 and j > 1:
                    i1j2 = table[i-1][j-2]

                opt1 = dist(ref_i, read_j) + diag
                opt2 = int(right / 1.5)
                opt3 = int(left / 1.5)

                # missing site
                opt4 = dist(ref_i + ref[i-1], read_j) + i2j1
                
                # extra site
                opt5 = dist(ref_i, read_j + read.cuts[j-1]) + i1j2

                # inversion
                opt6 = dist(ref_i, read.cuts[j-1]) + dist(ref[i-1], read_j) + diag2

                best_val = min(opt1, opt2, opt3, opt4, opt5, opt6, 0)

                if 0 == best_val:
                    table[i][j] = 0
                    path[i][j] = "s"
                elif opt1 == best_val:
                    table[i][j] = opt1
                    path[i][j] = "d"
                elif opt2 == best_val:
                    table[i][j] = opt2
                    path[i][j] = "|"
                elif opt3 == best_val:
                    table[i][j] = opt3
                    path[i][j] = "-"
                elif opt4 == best_val:
                    table[i][j] = opt4
                    path[i][j] = "missing_site"
                elif opt5 == best_val:
                    table[i][j] = opt5
                    path[i][j] = "extra_site"
                elif opt6 == best_val:
                    table[i-1][j-1] = 1 
                    table[i][j] = opt6
                    path[i][j] = "inversion"
                    path[i-1][j-1] = 'x'
                else:
                    table[i][j] = 0
                    path[i][j] = "s"


        islands = []
        for i, row in enumerate(table):
            for j, cell in enumerate(row):
                island = Island(len(row), len(table))
                if cell != 0:
                    table[i][j] = 0
                    island.set(i, j, cell)
                    queue = [(i, j)]
                    while queue:
                        x, y = queue.pop()
                        if x > 0 and table[x-1][y] != 0:
                            island.set(x-1, y, table[x-1][y])
                            table[x-1][y] = 0
                            queue.insert(0,(x-1,y))
                        if x+1 < len(table) and table[x+1][y] != 0:
                            island.set(x+1, y, table[x+1][y])
                            table[x+1][y] = 0
                            queue.insert(0,(x+1,y))
                        if y > 0 and table[x][y-1] != 0:
                            island.set(x, y-1, table[x][y-1])
                            table[x][y-1] = 0
                            queue.insert(0,(x,y-1))
                        if y+1 < len(row) and table[x][y+1] != 0:
                            island.set(x, y+1, table[x][y+1])
                            table[x][y+1] = 0
                            queue.insert(0,(x,y+1))
                        if y+1 < len(row)  and x+1 < len(table) and table[x+1][y+1] != 0:
                            island.set(x+1, y+1, table[x+1][y+1])
                            table[x+1][y+1] = 0
                            queue.insert(0,(x+1,y+1))
                    islands.append(island)

        for l in path:
            for e in l:
                print(e[0], end=' ')
            print()

        sorted_islands = sorted(islands, key=lambda i: i.get_best()[2])

        best_islands = []
        for island in sorted_islands:
            if not island.list_conflicts(best_islands):
                best_islands.append(island)

        chains = [self.backtrack_island(island, path, read, ref) for island in best_islands]

        print(chains)
        final = ["" for _ in range(len(ref))]
        for index in range(len(final)):
            final[index] = ''.join([chain[index] for chain in chains])

        counter = 0
        for index, distance in enumerate(final):
            if counter >= len(read.cuts):
                break

            if distance == str(read.cuts[counter]):
                counter += 1
            elif counter > 0 and distance == '':
                final[index] = 'd'
            
        print(final)
        return final


    def backtrack_island(self, island, path, read, ref):
        best_i, best_j, val = island.get_best()

        if val == 0:
            return None

        for l in island.grid:
            for e in l:
                print(f'{e:2}', end=' ')
            print()

        print()

        out = []
        while best_i >=0 and best_j >= 0:
            if path[best_i][best_j] == 'd':
                out.append(str(read.cuts[best_j]))
                best_i -= 1
                best_j -= 1
            elif path[best_i][best_j] == '-':
                out.append('*')
                best_j -= 1
            elif path[best_i][best_j] == '|':
                out.append('-')
                best_i -= 1
            elif path[best_i][best_j] == 'missing_site':
                out.append(f'{ref[best_i]}s')
                out.append(f'{ref[best_i-1]}s')
                best_i -= 2
                best_j -= 1
            elif path[best_i][best_j] == 'extra_site':
                out.append(f'{read.cuts[best_j-1]}+{read.cuts[best_j]}')
                best_i -= 1
                best_j -= 2
            elif path[best_i][best_j] == 'inversion':
                out.append(f'{read.cuts[best_j]}r')
                out.append(f'{read.cuts[best_j-1]}r')
                best_i -= 2
                best_j -= 2
            else:
                break

        out = out + [""] * (best_i + 1)
        out = out[::-1] + [""] * (len(ref) - len(out))
        return out

        # return postion and list of structural variations present

class Read:
    def __init__(self, start_pos):
        # list of cuts. 
        # This is all we have initially to map with
        self.cuts = []

        # list of structural variations (for comparison to guess)
        self.svs = []
        # actual location on reference (for comparison to guess)
        self.start_pos = start_pos

    def add_cuts(self, cuts):
        self.cuts += cuts


if __name__ == '__main__':
    args = docopt(__doc__)
    print(f"Creating a reference for chromosome {args['--cnum']}")
    ref = Reference(args['--cmap'], int(args['--cnum']))
    print(f'Found {ref.num_sites} sites and parsed {len(ref.positions)} sites')

    # generate some random reads:
    reads = [ref.generate_read() for _ in range(10)]
    print(f'Randomly generated {len(reads)} reads!')
    for r in reads:
        print(f'This read has {len(r.cuts)} cuts')

    #reads[0].cuts = [0, 10, 20, 30, 70, 80, 90]
    #ref.distances = [0, 10, 20, 30, 40, 50, 70, 80, 90]
    print(reads[0].cuts)
    print(ref.distances)
    final = ref.locate_read(reads[0])

    started = False
    for i in range(len(final)):
        if final[i] == 'd' and started:
            continue
        started = False
        if final[i] == 'd' and not started:
            print(f'...\t{final[i]}')
            started= True
            continue
        print(f'{ref.distances[i]}\t{final[i]}')
