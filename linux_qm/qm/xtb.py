from uuid import uuid4
import logging
import cclib
import re
import os
import shutil
from time import time
import io
import subprocess
import numpy as np
from rdkit import Chem
from rdkit.Chem import rdchem
from linux_qm.qm.driver import set_positions


class XTBDriver:
    # default options
    options = {
        'method': 'gfn2',
        'opt_geom': False,
        'solvent': None,
        'geom_maxiter': 100,
        'n_jobs': 1,
    }
    jobname: str = 'task'

    def __init__(
            self,
            options: dict = None,
            show_progress: bool = False,
            tmp_root: str = '.xtb_tmp',
    ):
        self.options = {**self.options, **(options or {})}
        self.tmp_root = tmp_root

    # Main methods

    def single_point(self, conf: rdchem.Conformer):
        """
        Single point property calculation.
        On success updates properties of the conformer.
        :param conf: rdkit Conformer
        :return: None
        """
        logging.info(f"Method: {self.options['method']}")
        self.options['opt_geom'] = False

        data = self.run(conf)

        success = data.metadata['success']
        conf.SetBoolProp('success', success)

        if success:
            # self.update_properties(conf, data)
            logging.info(f'Success: {success}')

        return data

    def geometry_optimization(self, conf: rdchem.Conformer):
        """
        Runs geometry optimization using self.options. On
        success - updates coordinates and properties of conformer
        :param conf: rdkit Conformer
        """
        start = time()
        logging.info(f"Method: {self.options['method']}")
        self.options['opt_geom'] = True
        self.options['calc_npa'] = calc_npa

        data = self.run(conf)

        success = data.metadata['success']
        conf.SetBoolProp('success', success)

        if success:
            # self.update_properties(conf, data)
            self.update_geometry(conf, data)
            logging.info(f'Success: {success}')
            logging.info(f'Num Iter: {len(data.atomcoords)}')
            logging.info(f'Elapsed Time: {time() - start:.1f}s')

        return data

    # Additional methods

    def run(self, conf: rdchem.Conformer):
        """
        Runs isolated xtb calculation in a separate tmp dir which is cleared
        after the run.
        Creates tmp dir, runs the calculation, parses the data, cleans the tmp dir.
        :param conf: rdkit Conformer
        :return: cclib parsed data
        """
        self._start_routine()

        input_str = self.gen_input(conf)
        output = self.run_xtb(input_str)
        data = self.parse(output)

        self._end_routine()

        return data

    def run_xtb(self, input_str):
        """
        Run xtb command line
        :param input_str:
        :return:
        """
        self._write_input(input_str)

        _cmd = f"""xtb {xyz_name} --{method} --opt {crit} --cycles {maxiter} {"--input param.inp" if xtbinp else ""} -P {n_jobs}"""
        # print(_cmd)
        res = subprocess.run(
            _cmd.split(' '),
            capture_output=True,
            text=True
        )
        output, error = res.stdout, res.stderr

        if res.returncode != 0:
            raise Exception(f'Error running orca: \nINPUT {input_str}\nSTDOUT: {output}\nSTDERR:{error}')

        last_lines = '\n'.join(output.splitlines()[-10:])
        logging.debug(f"ORCA OUTPUT:\n{last_lines}\n")
        # logging.debug(output)

        return output

    def parse(self, output: str):
        """
        Parses ORCA output. Additional parsing of other
        properties to be added.
        :param output: output from the orca run.
        :return: cclib parsed data
        """
        # with open('orca_output', 'w') as f:
        #     f.write(output)

        data = cclib.io.ccread(io.StringIO(output))

        logging.debug(f"CALC_NPA: {self.options['calc_npa']}")

        if self.options['calc_npa']:
            data.atomcharges['npa'] = self.run_junpa(data)
        # print(self.parse_janpa(res.stdout))
        # logging.debug(repr(self.parse_janpa(res.stdout)))

        return data

    def update_geometry(self, conf, cclib_data):
        set_positions(conf, cclib_data.atomcoords[-1])

    def update_properties(self, conf, cclib_data):
        conf.SetDoubleProp('energy', cclib_data.scfenergies[-1])

    # Auxilary methods


    def _create_tmp_dir(self):
        uid = str(uuid4())
        self.tmp_dir = f"{self.tmp_root}/{uid}"
        os.makedirs(self.tmp_dir, exist_ok=True)

    def _start_routine(self):
        # save current dir
        self.saved_work_dir = os.getcwd()
        self._create_tmp_dir()
        os.chdir(self.tmp_dir)

    def _end_routine(self):
        os.chdir(self.saved_work_dir)
        shutil.rmtree(self.tmp_dir)

    def gen_input(self, conf: rdchem.Conformer):
        xyz = Chem.MolToXYZBlock(conf)


    def _write_input(self, content: str):
        with open(self.jobname + '.inp', 'w') as f:
            f.write(content)
    #
    # def write_output(self, content: str):
    #     fname = 'output'
    #     with open(fname, 'w') as f:
    #         f.write(content)
    #     return fname

