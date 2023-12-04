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
from rdkit.Chem import rdchem
from linux_qm.qm.driver import set_positions

"""
Main ORCA driver.
The general idea is to have a rdkit Conformer as an input, to keep the atom connectivity information. After 
optimization or single point calculation coordinates and properties are updated. Molecular, atomic and bond
properties are set up for corresponding entity. 
```
properties are save as corresponding rdkit properties, ex: 
mol.SetDoubleProp('energy', value)
atom.SetDoubleProp('nmr_shielding', value)
bond.SetDoubleProp('ir_freq', value)
```
After each calculation for conformer specific boolean property is created `success`, which indicates reaching 
the scf-convergence and geometry optimization targets. Options are passed through the `options` dict argument and 
can be updated on the class instance. It is necessary to indicate the location of ORCA and OpenMPI installation folders.


Examples:

# Initialization
orca = OrcaDriver(
    orca_path='/opt/orca-5.0.3',
    omp_path='/opt/openmpi-4.1.1'
)

# read smiles and generate initial guess for 3D coordinates
mol = Chem.MolFromSmiles('COC')
mol = Chem.AddHs(mol)
AllChem.EmbedMolecule(mol)

# additional conformer generation pipeline may be added

conf = mol.GetConformer()
orca.geometry_optimization(conf)
conf.GetDoubleProp('energy') # get optimized energy
atom_pos = conf.GetPositions() # get optimized geometry
"""

class OrcaDriver:
    # default options
    options = {
        'method': 'B3LYP def2-SVP',
        'opt_geom': False,
        'calc_nmr': False,
        'calc_freq': False,
        'calc_npa': False,
        'solvent': None,
        'geom_maxiter': 100,
        'n_jobs': 1,
        'mem_per_core': 2000,  # in MB
    }
    jobname: str = 'task'

    def __init__(
            self,
            orca_path: str = '/opt/orca-5.0.3',
            omp_path: str = '/opt/openmpi-4.1.1',
            janpa_path: str = '/home/ergot/linux_qm/install/janpa',
            options: dict = None,
            show_progress: bool = False,
            tmp_root: str = '.orca_tmp'
    ):
        self.orca_path = orca_path
        self.omp_path = omp_path
        self.janpa_path = janpa_path
        self.options = {**self.options, **(options or {})}
        self.tmp_root = tmp_root
        self.set_orca_env()

    # Main methods

    def single_point(self, conf: rdchem.Conformer, calc_npa=False):
        """
        Single point property calculation.
        On success updates properties of the conformer.
        :param conf: rdkit Conformer
        :return: None
        """
        logging.info(f"Method: {self.options['method']}")
        self.options['opt_geom'] = False
        self.options['calc_npa'] = calc_npa

        data = self.run(conf)

        success = data.metadata['success']
        conf.SetBoolProp('success', success)

        if success:
            # self.update_properties(conf, data)
            logging.info(f'Success: {success}')

        return data

    def geometry_optimization(self, conf: rdchem.Conformer, calc_npa=False):
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
        Runs isolated orca calculation in a separate tmp dir which is cleared
        after the run.
        Creates tmp dir, runs the calculation, parses the data, cleans the tmp dir.
        :param conf: rdkit Conformer
        :return: cclib parsed data
        """
        self._start_routine()

        input_str = self.gen_input(conf)
        output = self.run_orca(input_str)
        data = self.parse(output)

        self._end_routine()

        return data

    def run_orca(self, input_str):
        """
        Run orca command line
        :param input_str:
        :return:
        """
        self._write_input(input_str)

        res = subprocess.run(
            [self.orca_path + '/orca', self.jobname + '.inp'],  #, '--use-hwthread-cpus --allow-run-as-root'], # '--use-hwthread-cpus --allow-run-as-root'
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

    def run_junpa(self, cclib_data):
        jobname = cclib_data.metadata['input_file_name'].strip('.inp')
        logging.debug(f'JANPA jobname:{jobname}')

        _cmd = f"orca_2mkl {jobname} -molden"
        logging.debug(_cmd)
        res = subprocess.run(
            _cmd.split(' '),
            capture_output=True,
            text=True
        )
        output, error = res.stdout, res.stderr

        if res.returncode != 0:
            raise Exception(f'Error running {_cmd}: \nSTDOUT: {output}\nSTDERR:{error}')
        if not os.path.exists(f'{jobname}.molden.input'):
            raise FileNotFoundError(f'Missing file {jobname}.molden.input: \nSTDOUT: {output}\nSTDERR:{error}')

        # logging.debug(str([f for f in os.listdir('.')]))
        # logging.debug(f'JANPA jobname:{jobname}')

        _cmd = f"java -jar {self.janpa_path}/molden2molden.jar -i {jobname}.molden.input -o joutput.PURE -fromorca3bf -orca3signs"

        logging.debug(_cmd)
        res = subprocess.run(
            _cmd.split(' '),
            capture_output=True,
            text=True
        )
        output, error = res.stdout, res.stderr

        if res.returncode != 0:
            raise Exception(f'Error running {_cmd}: \nSTDOUT: {output}\nSTDERR:{error}')
        if not os.path.exists('joutput.PURE'):
            raise FileNotFoundError(f'Missing file joutput.PURE: \nSTDOUT: {output}\nSTDERR:{error}')

        _cmd = f"java -jar {self.janpa_path}/janpa.jar -i joutput.PURE"

        logging.debug(_cmd)
        res = subprocess.run(
            _cmd.split(' '),
            capture_output=True,
            text=True
        )
        output, error = res.stdout, res.stderr

        if res.returncode != 0:
            raise Exception(f'Error running {_cmd}: \nSTDOUT: {output}\nSTDERR:{error}')

        npa_charges = self._parse_janpa(output)
        return npa_charges

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

    def _parse_janpa(self, output: str):
        lines = output.split('\n')
        lines = [l.strip() for l in lines if l.strip() != '']
        npa_lines = []
        npa_found = False
        for line in lines:
            if 'Final electron populations and NPA charges' in line:
                npa_found = True
            if 'Angular momentum contributions of the total atomic population' in line:
                break
            if npa_found:
                npa_lines.append(line)

        # logging.debug(npa_lines[1])
        # logging.debug(npa_lines[2])
        npa_charges = np.empty(len(npa_lines[3:]))
        for i, line in enumerate(npa_lines[3:]):
            logging.debug(line.split('\t')[4])
            charge = float(line.split('\t')[4])
            npa_charges[i] = round(charge, 6)
        return npa_charges

    def update_geometry(self, conf, cclib_data):
        set_positions(conf, cclib_data.atomcoords[-1])

    def update_properties(self, conf, cclib_data):
        conf.SetDoubleProp('energy', cclib_data.scfenergies[-1])

    # Auxilary methods

    def set_orca_env(self):
        """
        Sets up necessary env variables
        :return:
        """
        PATH = os.environ.get('PATH', '')
        LD_LIBRARY_PATH = os.environ.get('LD_LIBRARY_PATH', '')

        if self.orca_path not in PATH:
            PATH = f"{self.orca_path}:{PATH}"
            LD_LIBRARY_PATH = f"{self.orca_path}:{LD_LIBRARY_PATH}"
        if self.omp_path not in PATH:
            PATH = f"{self.omp_path}/bin:{PATH}"
            LD_LIBRARY_PATH = f"{self.omp_path}/lib:{LD_LIBRARY_PATH}"

        # # intel MKL
        # os.environ['LD_LIBRARY_PATH'] = f"/usr/lib/x86_64-linux-gnu:{os.environ['LD_LIBRARY_PATH']}"

        os.environ['PATH'] = PATH
        os.environ['LD_LIBRARY_PATH'] = LD_LIBRARY_PATH

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

    def _gen_xyz_block(self, conf: rdchem.Conformer):
        mol = conf.GetOwningMol()
        xyz_block = "*xyz 0 1\n"
        for i in range(conf.GetNumAtoms()):
            symbol = mol.GetAtomWithIdx(i).GetSymbol()
            x, y, z = conf.GetAtomPosition(i)
            xyz_block += f"{symbol:3}{x:20.8f}{y:20.8f}{z:20.8f}\n"
        xyz_block += "*\n"
        return xyz_block

    def gen_input(self, conf: rdchem.Conformer):
        # unpack options
        geom = 'OPT' if self.options['opt_geom'] else ''
        nmr = 'NMR' if self.options['calc_nmr'] else ''
        freq = 'FREQ' if self.options['calc_freq'] else ''
        method, solvent = self.options['method'], self.options['solvent']
        geom_maxiter = self.options['geom_maxiter']
        n_jobs, mem_per_core = self.options['n_jobs'], self.options['mem_per_core']

        input_str = ""

        input_str += f"!{method} {geom} {nmr} {freq}\n"
        input_str += f"%geom MaxIter {geom_maxiter} end\n"

        if solvent:
            input_str += (f"%cpcm\n"
                          f"  smd true\n"
                          f"  SMDsolvent \"{solvent}\"\n"
                          f"  end\n")

        # input_str += (f"%output\n "
        #     f"  PrintLevel Huge\n"
        #     # f"  Print[P_MOs] 1\n"
        #     # f"  Print[ P_gradient ] 1\n"
        #     f"  end\n")

        # parallel runs
        if n_jobs > 1:
            input_str += f"%pal nprocs {n_jobs} end\n"

        input_str += f"%maxcore {mem_per_core}\n"
        input_str += self._gen_xyz_block(conf)

        logging.debug(f"ORCA INPUT:\n{input_str}\n")
        return input_str

    def _write_input(self, content: str):
        with open(self.jobname + '.inp', 'w') as f:
            f.write(content)
    #
    # def write_output(self, content: str):
    #     fname = 'output'
    #     with open(fname, 'w') as f:
    #         f.write(content)
    #     return fname

