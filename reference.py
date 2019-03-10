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
from island import Island, CreateIslands
from table import Table
from read import Read, SV

class Reference:
    def __init__(self, filename, cnumber):
        self.cnumber = cnumber
        # read file, build positions and distances
        self.reference = self.parse_cmap(filename, cnumber)
        #self.distances = self.distances[:80]

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
        min_length, max_length = 20, 100
        min_delete, max_delete = 100, 1000
        expected_merges = 2 
        expected_splits = 2 
        expected_inversions = 2
        long_delete_prob = 1
        long_delete_length = 0
        long_delete_start = 0

        length = random.randint(min_length, max_length)

        if random.random() < long_delete_prob:
            long_delete_length = random.randint(min_delete, max_delete)
            long_delete_start = random.randint(length//4, 3*length//4)

        start_pos = random.randint(0, len(self.distances) - length - long_delete_length)

        read = Read(start_pos)

        if long_delete_start != 0:
            read.add_sv(SV('long_delete', long_delete_start, long_delete_length))

        i = 0
        ref = self.distances[start_pos:start_pos + length + long_delete_length]
        while i < len(ref):
            if i in range(long_delete_start, long_delete_start + long_delete_length):
                i = long_delete_start + long_delete_length
                continue

            seq = []
            if random.random() > 1.0 / expected_merges and i + 1 < length:
                merge = ref[i] + ref[i+1]
                seq = [merge]
                read.add_sv(SV('missing_site', i, merge ,{'pos': ref[i]}))
                expected_merges -= 1

            elif random.random() > 1.0 / expected_splits:
                split = ref[i] - random.randint(0, ref[i])
                seq = [split, ref[i] - split]
                read.add_sv(SV('extra_site', i, ref[i], {'pos': split}))
                expected_splits -= 1

            elif random.random() > 1.0 / expected_inversions and i + 1 < start_pos + length:
                seq = [ref[i+1], ref[i]]
                read.add_sv(SV('inversion', i, 2))
                expected_inversions -= 1

            else:
                seq = [ref[i]]

            read.add_cuts(seq)
            i += len(seq)

        # return Read object
        return read

    def locate_read(self, read):
        # find where in reference the read fits
        # local alignment means no penalty for shifting
        width, height = len(read.cuts), len(self.distances)
        mem = Table(width, height, 0)
        path = Table(width, height, '')

        ref = self.distances

        def dist(x, y):
            if x == y:
                return -1
            return 1

        for i, ref_i in enumerate(ref):
            for j, read_j in enumerate(read.cuts):
                opt4 = opt5 = opt6 = 0

                opt1 = dist(ref_i, read_j) + mem.get(i-1, j-1)
                opt2 = int(mem.get(i-1, j) / 1.5)
                opt3 = int(mem.get(i, j-1) / 1.5)

                # missing site
                if i > 0:
                    opt4 = dist(ref_i + ref[i-1], read_j) + mem.get(i-2, j-1)
                
                # extra site
                if j > 0:
                    opt5 = dist(ref_i, read_j + read.cuts[j-1]) + mem.get(i-1, j-1)

                # inversion
                if i > 0 and j > 0:
                    opt6 = dist(ref_i, read.cuts[j-1]) + dist(ref[i-1], read_j) + mem.get(i-2, i-2)

                best_val = min(opt1, opt2, opt3, opt4, opt5, opt6, 0)

                if 0 == best_val:
                    p = 's'
                elif opt1 == best_val:
                    p = 'd'
                elif opt2 == best_val:
                    p = '|'
                elif opt3 == best_val:
                    p = '-'
                elif opt4 == best_val:
                    if mem.get(i-1, j) == 0:
                        mem.set(i-1, j, 1)
                        path.set(i-1, j, 'x')
                    p = 'missing_site'
                elif opt5 == best_val:
                    if mem.get(i, j-1) == 0:
                        mem.set(i, j-1, 1)
                        path.set(i, j-1, 'x')
                    p = 'extra_site'
                elif opt6 == best_val:
                    if mem.get(i-1, j-1) == 0:
                        mem.set(i-1, j-1, 1)
                        path.set(i-1, j-1, 'x')
                    p = 'inversion'

                mem.set(i, j, best_val)
                path.set(i, j, p)

        islands = sorted(CreateIslands(mem), key=lambda i: i.best)

        best_islands = []
        for island in islands:
            if not island.list_conflicts(best_islands):
                best_islands.append(island)

        chains = [self.backtrack(i.max_x, i.max_y, path, read.cuts) for i in best_islands]
        final = [''.join([c[i] for c in chains]) for i in range(height)]

        counter = 0
        s = False
        for index, distance in enumerate(final):
            if counter >= len(read.cuts):
                break

            if distance != '':
                if 's' in distance:
                    s = not s
                    if s:
                        counter += 1
                elif '+' in distance:
                    counter += 2
                else:
                    counter += 1
            elif counter > 0 and distance == '':
                final[index] = 'd'
            
        read.add_svs(final)
        print(read.start_pos)
        print(read.start_guess)
        print(read.svs_guesses)
        return final
                           
    def backtrack(self, max_x, max_y, path, cuts):
        ref = self.distances

        translate = {'d': lambda: [str(cuts[j])],
                     '-': lambda: ['*'],
                     '|': lambda: ['-'],
                     'm': lambda: [f'{ref[i]}s', f'{ref[i+1]}s'],
                     'e': lambda: [f'{cuts[j-1]}+{cuts[j]}'],
                     'i': lambda: [f'{cuts[j]}r', f'{cuts[j-1]}r']}
        decrement = {'d' : (1, 1),
                     '-' : (0, 1),
                     '|' : (1, 0),
                     'm' : (2, 1),
                     'e' : (1, 2),
                     'i' : (2, 2)}

        out = []
        i, j = max_x, max_y
        while i >= 0 and j >= 0:
            c = path.get(i, j)[0]
            if c in ['s', 'x']:
                break

            out = translate[c]() + out
            di, dj = decrement[c]
            i, j = i - di, j - dj

        out = [""] * (i + 1) + out
        return out + [""] * (len(self.distances) - len(out))

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

    final = ref.locate_read(reads[0])
    print(reads[0].svs)

    started = False
    for i in range(len(final)):
        if (final[i] == 'd' or final[i] == '') and started:
            continue
        started = False
        if (final[i] == 'd' or final[i] == '') and not started:
            print(f'...\t{final[i]}')
            started= True
            continue
        print(f'{ref.distances[i]}\t{final[i]}')
