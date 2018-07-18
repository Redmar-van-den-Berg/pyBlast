#!/usr/bin/env python3

import unittest
import os
import shutil
import tempfile

import pyBlast

class test_pyblast(unittest.TestCase):
    def test_unsupported_ncbi_command(self):
        self.verbose = False
        from Bio.Blast.Applications import NcbiblastnCommandline
        from Bio.Application import ApplicationError
        cmd = NcbiblastnCommandline(
            query = 'test/query.fasta',
            db = 'test/target.fasta',
            evalue = 0.001
        )
        cmd.program_name = 'bla_blast'
        self.assertRaises(
            ApplicationError,
            pyBlast.pyBlast.run_blast,
            self, cmd
        )

class testPyBlast(unittest.TestCase):
    def setUp(self):
        from Bio.Blast.Applications import NcbiblastnCommandline

        cmd = NcbiblastnCommandline(
            query = 'test/sul2_1_AF542061.fasta',
            db = 'test/102637-001-018_k64-contigs.fa',
            evalue = 0.001
        )
        with pyBlast.pyBlast(cmd, rm_tmp=False) as pb:
            self.sul2_1 = next(pb)

    def test_true(self):
        self.assertTrue(True)

    def test_min_max_two(self):
        self.assertEqual(pyBlast._minmax(3,2), (2,3))

    def test_min_max_ten_list(self):
        numbers = list(range(10))
        self.assertEqual(pyBlast._minmax(numbers), (0,9))

    def test_min_max_ten(self):
        self.assertEqual(pyBlast._minmax(0,1,2,3,4,9), (0,9))

    def test_hit_overlap(self):
        self.assertFalse(pyBlast._hit_overlap(*self.sul2_1.alignments[0].hsps))

    def test_srange_five(self):
        self.assertEqual(pyBlast._srange(0,5), set([0,1,2,3,4]))

    def test_srange_zero(self):
        self.assertEqual(pyBlast._srange(0,0), set())

    def test_srange_one(self):
        self.assertEqual(pyBlast._srange(0,1), {0})


class test_file_operations(unittest.TestCase):
    def setUp(self):
        self.temp_dir = tempfile.mkdtemp()
        self.target = shutil.copy('test/target.fasta', self.temp_dir)
        self.query = shutil.copy('test/query.fasta', self.temp_dir)
        self.verbose = False

    def tearDown(self):
        shutil.rmtree(self.temp_dir)

    def test_temp_folder(self):
        self.assertTrue(os.path.isdir(self.temp_dir))

    def test_files_exist(self):
        self.assertTrue(os.path.isfile(self.target))
        self.assertTrue(os.path.isfile(self.query))

    def test_makeblastdb(self):
        blastdb = os.path.join(self.temp_dir,'testdb')
        pyBlast.pyBlast.makeblastdb(
            self,
            target = self.target,
            dbtype = 'nucl',
            blastdb = blastdb
        )
        self.assertTrue(os.path.isfile(blastdb + '.nhr'))
        self.assertTrue(os.path.isfile(blastdb + '.nin'))
        self.assertTrue(os.path.isfile(blastdb + '.nsq'))

    def test_unsupported_makeblastdb_command(self):
        self.assertRaises(
            FileNotFoundError,
            pyBlast.pyBlast.makeblastdb,
            self, self.target, 'nucl', self.temp_dir, makeblastdb='fakecommand'
        )

    def test_unsupported_makeblastdb_dbtype(self):
        self.assertRaises(
            SyntaxError,
            pyBlast.pyBlast.makeblastdb,
            self, self.target, 'wrong_dbtype', self.temp_dir
        )

    def test_missing_file_makeblastdb(self):
        self.assertRaises(
            FileNotFoundError,
            pyBlast.pyBlast.makeblastdb,
            self,'no_file','nucl', self.temp_dir
        )

    def test_makeblastdb_target_is_folder(self):
        self.assertRaises(
            FileNotFoundError,
            pyBlast.pyBlast.makeblastdb,
            self, self.temp_dir, 'nucl', self.temp_dir
        )


if __name__ == '__main__':
    unittest.main()
