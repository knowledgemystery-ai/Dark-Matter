import unittest

import numpy as np

from finite_mass_qed_transverse_factorized_spinor_reference import transverse_pair
from finite_mass_qed_transverse_real_pauli_reference import transverse_real_pauli


class RealPauliSpinSum(unittest.TestCase):
    def test_replays_complex_amplitudes_without_summing_amplitudes_over_spin(self):
        rng=np.random.default_rng(1254)
        M,m=.05,.045
        for w,k in ((.002304225,.0015),(.00150000001,.0015),(.006,.004)):
            a,b=[rng.normal(size=(27,8,3))*.65 for _ in range(2)]
            a[0]=[1e-10,1e-13,12.]
            b[0]=[-1e-11,1e-13,11.]
            p=np.concatenate((np.hypot(np.linalg.norm(a,axis=-1),M)[...,None],a),axis=-1)
            f=np.concatenate((np.hypot(np.linalg.norm(b,axis=-1),m)[...,None],b),axis=-1)
            K=np.array([w,0.,0.,k])
            u=transverse_pair(p,f,K,M,m)
            v=transverse_real_pauli(p,f,K,M,m)
            np.testing.assert_allclose(v[1:],u[1:],rtol=3e-12,atol=2e-20)
            np.testing.assert_allclose(v[0],u[0],rtol=1e-9,atol=2e-20)

    def test_rejects_singular_branch_scope(self):
        y=np.zeros((1,1,4))
        with self.assertRaises(ValueError):
            transverse_real_pauli(y,y,np.array([.2,0.,0.,.01]),.05,.045)


if __name__=='__main__':
    unittest.main()
