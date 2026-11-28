import itertools
import unittest

import numpy as np

from early_pole_endpoint_warm_streaming_pin_budget import solve_budget


class PinBudget(unittest.TestCase):
    def test_sorted_complement_is_exact_cardinality_minimum(self):
        loss=np.array([.8,.2,.5,.1,.3])
        rank,tail,minimum,sufficient=solve_budget(loss,.605)
        self.assertEqual(minimum,2)
        self.assertEqual(sufficient,2)
        self.assertEqual(rank.tolist(),[0,2,4,1,3])
        for i in range(len(loss)+1):
            best=min(sum(x for j,x in enumerate(loss) if j not in p)
                     for p in itertools.combinations(range(len(loss)),i))
            self.assertAlmostEqual(tail[i],best,places=14)

    def test_bounds_keep_numeric_near_zero_undecided(self):
        _,_,lo,hi=solve_budget([.75,.25],.25)
        self.assertEqual((lo,hi),(1,2))
        self.assertEqual(solve_budget([0.,0.],.125)[2:],(0,0))


if __name__=='__main__':
    unittest.main()
