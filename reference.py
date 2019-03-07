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

class Reference:
    def __init__(self, filename, cnumber):
        self.cnumber = cnumber
        # read file, build positions and distances
        self.reference = self.parse_cmap(filename, cnumber)

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
        length = random.randint(10, 100)
        start_pos = random.randint(0, len(self.distances))

        read = Read(start_pos)
        read.add_cuts(self.distances[start_pos:start_pos+length])

        # return Read object
        return read

    def locate_read(self, read):
        # find where in reference the read fits

        # return postion and list of structural variations present
        pass

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
