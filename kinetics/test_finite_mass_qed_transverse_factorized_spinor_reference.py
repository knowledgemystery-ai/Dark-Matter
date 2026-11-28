import unittest

import numpy as np

from early_pole_hard_yukawa_photoproduction_reference import norm_spinors
from finite_mass_qed_timelike_pole_batched_absorption_reference import reduce_vertex
from finite_mass_qed_transverse_factorized_spinor_reference import transverse_pair


class FactorizedTransverseSpinSum(unittest.TestCase):
    def test_replays_both_transverse_polarizations_at_each_node(self):
        rng=np.random.default_rng(417)
        M,m=.05,.045
        for w,k in ((.002304225,.0015),(.00150000001,.0015),(.006,.004)):
            a,b=[rng.normal(size=(12,8,3))*.5 for _ in range(2)]
            a[0]=[1e-10,1e-13,12.]
            b[0]=[-1e-11,1e-13,11.]
            p=np.concatenate((np.hypot(np.linalg.norm(a,axis=-1),M)[...,None],a),axis=-1)
            f=np.concatenate((np.hypot(np.linalg.norm(b,axis=-1),m)[...,None],b),axis=-1)
            K=np.array([w,0.,0.,k])
            u=sum(norm_spinors(p,f,M,m,reduce_vertex(p,f,K,M,m,e)) for e in np.eye(4)[1:3])
            v=transverse_pair(p,f,K,M,m)
            np.testing.assert_allclose(v[1:],u[1:],rtol=3e-12,atol=2e-20)
            # Deliberately extreme nearly collinear rows subtract large amplitude terms.
            np.testing.assert_allclose(v[0],u[0],rtol=3e-10,atol=2e-20)

    def test_rejects_misaligned_or_nontimelike_pole(self):
        y=np.zeros((1,1,4))
        for k in ([.01,.001,0.,.004],[.01,0.,0.,.01]):
            with self.assertRaises(ValueError):
                transverse_pair(y,y,np.array(k),.05,.045)


if __name__=='__main__':
    unittest.main()
