from itertools import groupby
from quantum_codon_opt.src.constants import *


class SeqScorer(object):
    def __init__(self, n_seq):
        self.n_seq = n_seq
        self.gc_constant = gc_constant
        self.repeat_constant = repeat_constant
        self.rarity_constant = rarity_constant
        self.execute()

    def __repr__(self):
        return 'Classical scoring function.'

    def execute(self):
        self._gc_score()
        self._repeat_score()
        self._rarity_score()
        self.score = self.gc_constant * self.gc_score + self.repeat_constant * self.rep_score + self.rarity_constant * self.rarity_score

    def _gc_score(self):
        self.gc_score = ((self.n_seq.count('C') + self.n_seq.count('G')) /
                         float(len(self.n_seq)) - 0.5)**2

    def _repeat_score(self):
        score = 0
        for _i, c in enumerate(list(self.n_seq)[:-3][::3]):
            ind = 3 * _i
            s = self.n_seq[ind:ind + 6]
            repeats = sorted([(letter, len(list(group)))
                              for letter, group in groupby(s)],
                             key=lambda i: i[1],
                             reverse=True)
            score += max(repeats, key=lambda i: i[1])[1]**2 - 1
        self.rep_score = score

    def _rarity_score(self):
        self.rarity_score = sum([
            codon_scores[self.n_seq[i:i + 3]]
            for i in range(len(self.n_seq))[::3]
        ])
