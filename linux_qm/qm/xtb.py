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
from linux_qm.src.util import set_atom_positions


COORD = 'coord.xyz'
INPUT = 'xtb.inp'
TMP_ROOT = '.xtb_tmp'

class XTBDriver:
    # default options
    options = {
        'method': 'gfn 2',
        'opt_geom': False,
        'opt_level': 'normal',
        'opt_engine': None,
        'opt_maxiter': None,
        'etemp': None,
        'input': None,
        'constrain': False,
        'solvent': None,
        'n_jobs': 1,
    }

    def __init__(
            self,
            options: dict = None,
            show_progress: bool = False,
    ):
        self.options = {**self.options, **(options or {})}

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
            self.update_properties(conf, data)
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

        data = self.run(conf)

        success = data.metadata['success']
        conf.SetBoolProp('success', success)

        if success:
            self.update_properties(conf, data)
            self.update_geometry(conf, data)

            num_iter = len(data.scfenergies)
            elapsed = time() - start
            time_per_iter = elapsed / num_iter

            logging.info(f'Success: {success}')
            logging.info(f'Num Iter: {num_iter}')
            logging.info(f'Elapsed Time: {elapsed:.1f}s')
            logging.info(f'Time per Iter: {time_per_iter:.1f}s')


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
        self._init()

        self._write_input(conf)
        output = self.run_xtb()
        data = self.parse(output)
        self._clean()

        return data

    def run_xtb(self):
        """
        Run xtb command line
        :param input_str:
        :return:
        """
        method = self.options['method']
        opt = self.options['opt_geom']
        level = self.options['opt_level']
        maxiter = self.options['opt_maxiter']
        inp = self.options['input']
        etemp = self.options['etemp']
        solv = self.options['solvent']
        n_jobs = self.options['n_jobs']

        _cmd = f"""xtb {COORD} --{method} \
        {"--input " + INPUT if inp else ""} \
        {"--opt " + level if opt else ""} \
        {"--cycles " + str(maxiter) if maxiter else ""} \
        {"--alpb " + solv if solv else ""} \
        {"--etemp " + str(etemp) if etemp else ""} \
        -P {n_jobs} \
        --verbose"""
        _cmd = re.sub(' +', ' ', _cmd)

        logging.debug(_cmd)

        res = subprocess.run(
            _cmd.split(' '),
            capture_output=True,
            text=True
        )

        if res.returncode != 0:
            raise Exception(f'Error running xtb: \nINPUT {_cmd}\nSTDOUT: {res.stdout}\nSTDERR:{res.stderr}')

        # logging.debug(res.stdout)

        return res.stdout

    def parse(self, output: str):
        """
        Parses ORCA output. Additional parsing of other
        properties to be added.
        :param output: output from the orca run.
        :return: cclib parsed data
        """

        data = cclib.io.ccread(io.StringIO(output))

        if self.options['opt_geom']:
            with open('xtbopt.xyz', 'r') as f:
                xyz = f.read()

        # custom parse before cclib fix double letter elements
        if xyz:
            lines = xyz.split('\n')
            match = re.search(r'energy:\s*(-?\d+\.\d+)', lines[1])
            if match:
                energy = float(match.group(1))
                data.scfenergies = [energy]
            atom_pos = []
            for line in lines[2:]:
                coords = re.findall('(-?\d+\.\d+)', line)
                coords = [float(x) for x in coords]
                if coords:
                    atom_pos.append(coords)
            data.atomcoords = [atom_pos]

        return data


    def update_geometry(self, conf, cclib_data):
        set_atom_positions(conf, cclib_data.atomcoords[-1])


    def update_properties(self, conf, cclib_data):
        conf.SetDoubleProp('energy', cclib_data.scfenergies[-1])

    # Auxilary methods


    def _create_tmp_dir(self):
        uid = str(uuid4())
        self.tmp_dir = f"{TMP_ROOT}/{uid}"
        os.makedirs(self.tmp_dir, exist_ok=True)

    def _init(self):
        # save current dir
        self.saved_work_dir = os.getcwd()
        self._create_tmp_dir()
        os.chdir(self.tmp_dir)

    def _clean(self):
        os.chdir(self.saved_work_dir)
        shutil.rmtree(self.tmp_dir)

    def _write_input(self, conf: rdchem.Conformer):
        conf_id = conf.GetId()
        mol = conf.GetOwningMol()
        xyz = Chem.MolToXYZBlock(mol, conf_id)

        with open(COORD, 'w') as f:
            f.write(xyz)

        if self.options['input']:
            with open(INPUT, 'w') as f:
                f.write(self.options['input'])
    #
    # def write_output(self, content: str):
    #     fname = 'output'
    #     with open(fname, 'w') as f:
    #         f.write(content)
    #     return fname

