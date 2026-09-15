# Written by Arunav Kumar, MIT Plasma Science and Fusion center, 10th May, 2026
from pathlib import Path
import tempfile
import unittest
import warnings

import f90nml
import numpy as np
from eqdsk2vmec import convert_eqdsk_to_vmec, read_geqdsk, eqdisk2vmec_inputfile
from eqdsk2vmec.comparison import reconstruct_boundary, plot_comparison
from eqdsk2vmec.vmec_namelist import VMECNamelist, write_vmec_input

ROOT = Path(__file__).resolve().parents[1]
FIXTURES = sorted((ROOT / 'Device').glob('*/*.geq'))


def synthetic():
    u = np.linspace(0, 1, 33)
    theta = np.linspace(0, 2*np.pi, 1024, endpoint=False)
    theta += .35*np.sin(theta)  # Deliberately nonuniform boundary sampling.
    return dict(sf=np.ones(33), sp=10*(1-u), qpsi=1+u, psiaxis=0., psilim=1.,
                btor=2., totcur=1e6, xaxis=2., zaxis=0.,
                xbndry=2+.5*np.cos(theta), zbndry=.5*np.sin(theta))


class TestNumerics(unittest.TestCase):
    def convert(self, data, **kwargs):
        return eqdisk2vmec_inputfile('synthetic', data, cocos=5, **kwargs)

    def test_analytic_flux_integral_both_directions(self):
        d = synthetic()
        u = np.linspace(0, 1, 33)
        for axis, edge in [(0., 1.), (3., 2.)]:
            d.update(psiaxis=axis, psilim=edge)
            a = self.convert(d)
            np.testing.assert_allclose(a.tflux, (u+.5*u*u)/1.5, atol=1e-14)
            self.assertAlmostEqual(a.phiedge, 3*np.pi)

    def test_polynomial_endpoints(self):
        d = synthetic()
        d['qpsi'][:] = 2
        a = self.convert(d)
        for coefficients, values in [(a.am, a.am_aux_f), (a.ai, a.ai_aux_f)]:
            self.assertAlmostEqual(coefficients[0], values[0], places=10)
            self.assertAlmostEqual(sum(coefficients), values[-1], places=10)
        np.testing.assert_allclose(np.polynomial.polynomial.polyval([0, .5, 1], a.am), [10, 5, 0], atol=1e-10)
        np.testing.assert_allclose(np.polynomial.polynomial.polyval([0, .5, 1], a.ai), .5, atol=1e-10)

    def test_nonuniform_circle_and_reversed_order(self):
        d = synthetic()
        for reverse in (False, True):
            if reverse:
                d['xbndry'] = d['xbndry'][::-1]
                d['zbndry'] = d['zbndry'][::-1]
            a = self.convert(d)
            self.assertAlmostEqual(a.refou[0, 0], 2, places=5)
            self.assertAlmostEqual(a.refou[1, 0], .5, places=5)
            self.assertAlmostEqual(a.zefou[1, 0], .5, places=5)
            self.assertLess(np.max(np.abs(a.refou[2:])), 1e-5)

    def test_offset_boundary_preserved_only_when_asymmetric(self):
        d = synthetic()
        d['zbndry'] += .2
        d['zaxis'] = .2
        a = self.convert(d, lasym=True)
        self.assertAlmostEqual(a.zefou2[0, 0], .2, places=5)
        with self.assertWarnsRegex(UserWarning, 'LASYM=0'):
            symmetric = self.convert(d)
        self.assertEqual(np.max(np.abs(symmetric.zefou2)), 0)

    def test_cocos_equivalent_representations(self):
        reference = eqdisk2vmec_inputfile('test', synthetic(), cocos=5)
        for cocos in (*range(1, 9), *range(11, 19)):
            with self.subTest(cocos=cocos):
                d = synthetic()
                base = cocos % 10
                if base in (1, 2, 7, 8):
                    d['qpsi'] *= -1
                if base % 2 == 0:
                    d['btor'] *= -1
                    d['totcur'] *= -1
                if base in (3, 4, 7, 8):
                    d['psilim'] *= -1
                if cocos >= 10:
                    d['psilim'] *= 2*np.pi
                a = eqdisk2vmec_inputfile('test', d, cocos=cocos)
                self.assertAlmostEqual(a.phiedge, reference.phiedge)
                self.assertAlmostEqual(a.curtor, reference.curtor)
                np.testing.assert_allclose(a.iotaf, reference.iotaf)
                np.testing.assert_allclose(a.tflux, reference.tflux)

    def test_invalid_profiles_rejected(self):
        for key, value in [('psilim', 0.), ('btor', 0.), ('qpsi', np.zeros(33)),
                           ('sp', np.full(33, np.nan)), ('qpsi', np.linspace(-1, 1, 33))]:
            with self.subTest(key=key, value=str(value)[:20]):
                d = synthetic()
                d[key] = value
                with self.assertRaises(ValueError):
                    self.convert(d)

    def test_unspecified_cocos_warns(self):
        with self.assertWarnsRegex(UserWarning, 'COCOS is unspecified'):
            eqdisk2vmec_inputfile('test', synthetic())


class TestSerialization(unittest.TestCase):
    def test_axis_and_pure_asymmetric_mode(self):
        v = VMECNamelist()
        v.lasym = True
        v.raxis_cc = [2.]
        v.zaxis_cc = [.2]
        v.zaxis_cs = [0.]
        v.rbc = np.zeros((1, v.mpol)); v.rbc[0, 0] = 2
        v.zbc = np.zeros((1, v.mpol)); v.zbc[0, 0] = .2
        with tempfile.TemporaryDirectory() as t:
            p = Path(t)/'input.test'
            write_vmec_input(p, v)
            n = f90nml.read(p)['indata']
            self.assertEqual(n['raxis_cc'], 2.)
            self.assertEqual(n['zaxis_cc'], .2)
            r, z = reconstruct_boundary(n, np.linspace(0, 2*np.pi, 50))
            np.testing.assert_allclose(r, 2)
            np.testing.assert_allclose(z, .2)

    def test_invalid_current_mode_and_array_shape(self):
        with tempfile.TemporaryDirectory() as t:
            p = Path(t)/'input.test'
            v = VMECNamelist(); v.ncurr = 1
            with self.assertRaisesRegex(ValueError, 'current profile'):
                write_vmec_input(p, v)
            self.assertFalse(p.exists())
            v.ncurr = 0; v.rbc = np.zeros((1, 1))
            with self.assertRaisesRegex(ValueError, 'shape'):
                write_vmec_input(p, v)

    def test_real_fixtures_written_and_compared(self):
        for fixture in FIXTURES:
            for lasym in (False, True):
                with self.subTest(file=fixture.name, lasym=lasym), tempfile.TemporaryDirectory() as t:
                    p = Path(t)/'input.test'
                    with warnings.catch_warnings():
                        warnings.simplefilter('ignore', UserWarning)
                        convert_eqdsk_to_vmec(fixture, output_path=p, cocos=5, lasym=lasym)
                        d = read_geqdsk(fixture)
                        a = eqdisk2vmec_inputfile(fixture, d, cocos=5, lasym=lasym)
                    n = f90nml.read(p)['indata']
                    self.assertEqual(n['ntor'], 0)
                    self.assertEqual(n['lasym'], lasym)
                    self.assertEqual(n['ncurr'], 0)
                    self.assertNotIn('ac_aux_f', n)
                    self.assertLessEqual(len(n['am_aux_s']), 101)
                    self.assertAlmostEqual(n['raxis_cc'], d['xaxis'])
                    np.testing.assert_allclose(n['am_aux_f'], d['sp'][a.profile_indices], rtol=1e-11)
                    np.testing.assert_allclose(1/np.array(n['ai_aux_f']), d['qpsi'][a.profile_indices], rtol=1e-11)
                    for name in ('rbc', 'zbs', 'rbs', 'zbc'):
                        if name in n:
                            self.assertLess(n.start_index[name][1]+len(n[name]), n['mpol']+1)
                    plot = Path(t)/'comparison.png'
                    stats = plot_comparison(d, p, a, plot_path=plot)
                    self.assertTrue(plot.exists())
                    self.assertTrue(plot.with_suffix('.json').exists())
                    self.assertLess(stats['q_max_knot_error'], 1e-9)
                    # Numerical fit acceptance, not equilibrium-solver validation.
                    self.assertLess(stats['boundary_rms_distance_m'], .04)


class TestReader(unittest.TestCase):
    def test_truncation_and_bad_numeric_field(self):
        source = FIXTURES[0].read_text()
        for text in [source[:90], source[:81]+'BADVALUE'+source[89:]]:
            with tempfile.TemporaryDirectory() as t:
                p = Path(t)/'bad.geq'; p.write_text(text)
                with self.assertRaises(ValueError):
                    read_geqdsk(p)

    def test_fortran_d_exponents(self):
        original = read_geqdsk(FIXTURES[0])
        with tempfile.TemporaryDirectory() as t:
            p = Path(t)/'d.geq'
            lines = FIXTURES[0].read_text().splitlines(keepends=True)
            p.write_text(lines[0]+''.join(line.replace('E', 'D') for line in lines[1:]))
            parsed = read_geqdsk(p)
            np.testing.assert_array_equal(parsed['psixz'], original['psixz'])
            np.testing.assert_array_equal(parsed['qpsi'], original['qpsi'])


if __name__ == '__main__':
    unittest.main()
