
class Read:
    def __init__(self, start_pos):
        # list of distances
        self.cuts = []

        # list of structural variations (for comparison to guess)
        self.svs = []
        self.svs_guesses = []
        # actual location on reference (for comparison to guess)
        self.start_pos = start_pos
        self.start_guess = 0

    def add_cuts(self, cuts):
        self.cuts += cuts

    def add_sv(self, sv, guess=False):

        if guess:
            sv.contents['global'] = self.start_guess + sv.contents['local']
            self.svs_guesses.append(sv)
        else:
            sv.contents['global'] = self.start_pos + sv.contents['local']
            self.svs.append(sv)

    def add_svs(self, final):
        i = 0
        start = -1
        while i < len(final):
            if '' == final[i]:
                i += 1
                continue

            if start == -1:
                self.start_guess = i

            if 's' in final[i]:
                first = int(final[i][:-1]) 
                length = first + int(final[i+1][:-1])
                self.add_sv(SV('missing_site', i-start, length, {'pos': first}), True)
                i += 2

            elif '+' in final[i]:
                first, second = final[i].split('+')
                self.add_sv(SV('extra_site', i-start, int(first) + int(second), {'pos': int(first)}), True)
                i += 1

            elif 'r' in final[i]:
                saved_i = i
                length = 0
                while i < len(final) and 'r' in final[i]:
                    length += 1
                    i += 1

                self.add_sv(SV('inversion', saved_i-start, length), True)

            elif 'd' == final[i]:
                saved_i = i
                length = 0
                while i < len(final) and 'd' == final[i]:
                    length += 1
                    i += 1

                self.add_sv(SV('long_delete', saved_i-start, length), True)
            else:
                i += 1
  
class SV:
    def __init__(self, name, local_pos, length, dic={}):
        self.name = name
        self.contents = dic
        self.contents['local'] = local_pos
        self.contents['length'] = length
