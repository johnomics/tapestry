import unittest, sys, os
import logging as log
from contextlib import contextmanager
from io import StringIO

from Tapestry.misc import get_args
from Tapestry._version import __version__

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
            get_args(['-V'])
            self.assertEqual(cm.exception.code, None)
            output = out.getvalue().strip().split('\n')
            self.assertEqual(output[0], f'Tapestry version {__version__}')

    def test_assembly_missing(self):
        with captured_output() as (out, err), self.assertRaises(SystemExit) as cm:
            get_args(['-c', '15'])
            self.assertEqual(cm.exception.code, None)
            output = err.getvalue().strip()
            self.assertEqual(output[-16:], '(-a, --assembly)')

    def test_zero_cores(self):
        with captured_output() as (out, err), self.assertRaises(SystemExit) as cm:
            get_args(['-a', 'test/test_assembly.fasta', '-c', '0'])
            self.assertEqual(cm.exception.code, None)
            output = err.getvalue().strip()
            self.assertEqual(output[-30:], 'Must specify at least one core')

    def test_too_many_cores(self):
        with captured_output() as (out, err), self.assertRaises(SystemExit) as cm:
            get_args(['-a', 'test/test_assembly.fasta', '-c', str(len(os.sched_getaffinity(0))+1)])
            self.assertEqual(cm.exception.code, None)
            output = err.getvalue().strip()
            self.assertEqual(output[-19:], 'please reduce cores')

    def test_verbosity_warn(self):
        get_args(['-a', 'test/test_assembly.fasta'])
        self.assertEqual(log.getLogger().level, 30) # WARNING

    def test_verbosity_info(self):
        get_args(['-a', 'test/test_assembly.fasta', '-v'])
        self.assertEqual(log.getLogger().level, 20) # INFO

    def test_verbosity_debug(self):
        get_args(['-a', 'test/test_assembly.fasta', '-vv'])
        self.assertEqual(log.getLogger().level, 10) # DEBUG