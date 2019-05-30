import unittest, sys, os
from contextlib import contextmanager
from io import StringIO
from shutil import rmtree
from argparse import ArgumentError

from tapestry.assembly import Assembly
from tapestry.misc import get_args, get_weave_args
from tapestry._version import __version__

#https://stackoverflow.com/questions/4219717/how-to-assert-output-with-nosetest-unittest-in-python
@contextmanager
def captured_output():
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err

class TestArgs(unittest.TestCase):

    def test_noargs(self):
        with captured_output() as (out, err), self.assertRaises(SystemExit) as cm:
            get_args()
        self.assertEqual(cm.exception.code, None)
        output = out.getvalue().strip()
        self.assertEqual(output[:6], 'usage:')

    def test_help(self):
        with captured_output() as (out, err), self.assertRaises(SystemExit) as cm:
            get_args('-h')
        self.assertEqual(cm.exception.code, 2)
        output = err.getvalue().strip()
        self.assertEqual(output[:6], 'usage:')

    def test_version(self):
        with captured_output() as (out, err), self.assertRaises(SystemExit) as cm:
            get_args(['-v'])
        self.assertEqual(cm.exception.code, 0)
        output = out.getvalue().strip().split('\n')
        self.assertEqual(output[0], f'Tapestry version {__version__}')

    def test_weave_assembly_reads_missing(self):
        with captured_output() as (out, err), self.assertRaises(SystemExit) as cm:
            get_weave_args(['-c', '2'])
        self.assertEqual(cm.exception.code, 2)
        output = err.getvalue().strip()
        self.assertEqual(output[-35:], 'required: -a/--assembly, -r/--reads')

    def test_weave_reads_missing(self):
        with captured_output() as (out, err), self.assertRaises(SystemExit) as cm:
            get_weave_args(['-a', 'test/test_assembly.fasta', '-c', '2'])
        self.assertEqual(cm.exception.code, 2)
        output = err.getvalue().strip()
        self.assertEqual(output[-20:], 'required: -r/--reads')

    def test_weave_assembly_missing(self):
        with captured_output() as (out, err), self.assertRaises(SystemExit) as cm:
            get_weave_args(['-r', 'test/test_reads.fastq.gz', '-c', '2'])
        self.assertEqual(cm.exception.code, 2)
        output = err.getvalue().strip()
        self.assertEqual(output[-23:], 'required: -a/--assembly')

    def test_weave_zero_cores(self):
        with self.assertRaises(SystemExit) as scm, self.assertLogs() as lcm:
            get_weave_args(['-a', 'test/test_assembly.fasta', '-r', 'test/test_reads.fastq.gz', '-c', '0'])
        self.assertEqual(scm.exception.code, None)
        self.assertEqual(lcm.output, ['ERROR:root:Please specify at least one core'])

    def test_weave_assembly_file_missing(self):
        with self.assertRaises(SystemExit) as scm, self.assertLogs() as lcm:
            get_weave_args(['-a', 'missing.fasta', '-r', 'test/test_reads.fastq.gz'])
        self.assertEqual(scm.exception.code, None)
        self.assertEqual(lcm.output, ['ERROR:root:Assembly file missing.fasta does not exist'])

    def test_weave_read_file_missing(self):
        with self.assertRaises(SystemExit) as scm, self.assertLogs() as lcm:
            get_weave_args(['-a', 'test/test_assembly.fasta', '-r', 'missing.fastq.gz'])
        self.assertEqual(scm.exception.code, None)
        self.assertEqual(lcm.output, ['ERROR:root:Reads file missing.fastq.gz does not exist'])

    def test_weave_read_file_unzipped(self):
        reads_file = 'test/test_reads.fastq'
        with self.assertRaises(SystemExit) as scm, self.assertLogs() as lcm:
            get_weave_args(['-a', 'test/test_assembly.fasta', '-r', reads_file])
        self.assertEqual(scm.exception.code, None)
        self.assertEqual(lcm.output, [f'ERROR:root:Reads file {reads_file} must be gzipped'])

    def test_run(self):
        if os.path.exists('test_run'):
            rmtree('test_run')
        with captured_output() as (out, err):
            assembly = Assembly('test/test_assembly.fasta', 'test/test_reads.fastq.gz', None, 'test_run', 1, 0, 10000, 10000, False)
            assembly.report()
        self.assertTrue(os.path.exists("test_run/test_run.tapestry_report.html"))
        self.assertTrue(os.path.exists("test_run/contig_details.tsv"))
        with open("test_run/contig_details.tsv") as details:
            lines = details.readlines()
        self.assertEqual(lines[-1][:25], "1	test_run_contig1	4	50.0")
        if os.path.exists('test_run'):
            rmtree('test_run')

    def test_run_zipped_assembly(self):
        if os.path.exists('test_run'):
            rmtree('test_run')
        with captured_output() as (out, err):
            assembly = Assembly('test/test_assembly.fasta.gz', 'test/test_reads.fastq.gz', None, 'test_run', 1, 0, 10000, 10000, False)
            assembly.report()
        self.assertTrue(os.path.exists("test_run/test_run.tapestry_report.html"))
        self.assertTrue(os.path.exists("test_run/contig_details.tsv"))
        with open("test_run/contig_details.tsv") as details:
            lines = details.readlines()
        self.assertEqual(lines[-1][:25], "1	test_run_contig1	4	50.0")
        if os.path.exists('test_run'):
            rmtree('test_run')

if __name__ == '__main__':
    unittest.main()