"""Synthetic contract tests; no network or claims about real cohorts."""
import csv
import json
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

COURSE = Path(__file__).resolve().parents[1]


@unittest.skipUnless(shutil.which('Rscript'), 'Rscript is required')
class AnalysisContracts(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='ximcs-contract-')
        self.addCleanup(self.temp.cleanup)
        self.course = Path(self.temp.name) / 'course'
        self.scripts = self.course / 'scripts'
        self.scripts.mkdir(parents=True)
        for name in ('_common.R', '11_geo_limma.R', '14_meta_analysis.R'):
            shutil.copy2(COURSE / 'scripts' / name, self.scripts / name)
        lib = COURSE / 'env/R-library'
        if lib.is_dir():
            (self.course / 'env').mkdir()
            (self.course / 'env/R-library').symlink_to(lib, target_is_directory=True)

    def run_script(self, name, *args):
        return subprocess.run(['Rscript', str(self.scripts / name), *map(str, args)],
                              capture_output=True, text=True, timeout=90)

    def geo_fixture(self):
        # Deterministic log2 values with a known paired mean/standard error.
        import math
        self.rows = [['feature_id', *[f'S{i}' for i in range(1, 7)]]]
        for g in range(30):
            control = [5 + math.sin(g + j) + j * .2 for j in range(3)]
            diff = [1 + .1 * g, 1.4 + .1 * g, .7 + .1 * g]
            self.rows.append([f'probe_{g}', *control,
                              *[c + d for c, d in zip(control, diff)]])
        self.expr = self.course / 'expression.csv'
        with self.expr.open('w', newline='') as f:
            csv.writer(f).writerows(self.rows)
        self.samples = [dict(sample_id=f'S{i+1}', condition='control' if i < 3 else 'treated',
                             donor_id=f'D{i % 3}') for i in range(6)]
        self.meta = self.course / 'samples.csv'

    def write_samples(self):
        with self.meta.open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(self.samples[0]))
            w.writeheader(); w.writerows(self.samples)

    def geo(self, mode='paired'):
        self.write_samples()
        return self.run_script('11_geo_limma.R', self.expr, self.meta,
                               'control', 'treated', mode, 'log2_verified')

    def test_paired_effect_and_sampling_se(self):
        import statistics
        self.geo_fixture()
        self.samples.reverse()  # Ordering must be based on IDs, not position.
        result = self.geo()
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.course / 'results/geo_limma/all_features.csv').open() as f:
            table = {r['feature_id']: r for r in csv.DictReader(f)}
        d = [1, 1.4, .7]
        self.assertAlmostEqual(float(table['probe_0']['logFC']), statistics.mean(d), places=10)
        self.assertAlmostEqual(float(table['probe_0']['SE_classical']),
                               statistics.stdev(d) / 3 ** .5, places=10)

    def test_independent_effect(self):
        self.geo_fixture()
        for i, row in enumerate(self.samples):
            row['donor_id'] = f'person_{i}'
        result = self.geo('independent')
        self.assertEqual(result.returncode, 0, result.stderr)
        with (self.course / 'results/geo_limma/all_features.csv').open() as f:
            table = {r['feature_id']: r for r in csv.DictReader(f)}
        self.assertAlmostEqual(float(table['probe_0']['logFC']), (1 + 1.4 + .7) / 3)

    def test_mismatched_sample_rejected(self):
        self.geo_fixture(); self.samples[0]['sample_id'] = 'not_in_expression'
        result = self.geo()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('one-to-one', result.stderr)

    def test_confounding_rejected(self):
        self.geo_fixture()
        for row in self.samples:
            row['batch'] = row['condition']
        result = self.geo()
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Nonidentifiable', result.stderr)

    def meta_input(self, mutate):
        with (COURSE / 'examples/meta_demo_effects.csv').open() as f:
            rows = list(csv.DictReader(f))
        mutate(rows)
        path = self.course / 'effects.csv'
        with path.open('w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
        return self.run_script('14_meta_analysis.R', path)

    def test_overlapping_cohort_rejected(self):
        result = self.meta_input(lambda r: r[1].update(cohort_id=r[0]['cohort_id']))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Repeated study/cohort', result.stderr)

    def test_incompatible_scale_rejected(self):
        result = self.meta_input(lambda r: r[0].update(effect_scale='SMD'))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('common scale', result.stderr)

    def test_unreviewed_overlap_rejected(self):
        result = self.meta_input(lambda r: r[0].update(independence_reviewed='FALSE'))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('Review participant overlap', result.stderr)

    def test_nonpositive_se_rejected(self):
        result = self.meta_input(lambda r: r[0].update(se='0'))
        self.assertNotEqual(result.returncode, 0)
        self.assertIn('positive sampling', result.stderr)


if __name__ == '__main__':
    unittest.main()
