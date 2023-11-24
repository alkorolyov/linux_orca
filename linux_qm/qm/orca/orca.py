from uuid import uuid4
import logging
import cclib
import os
import shutil
from time import time
import io
import subprocess
from rdkit.Chem import rdchem
from linux_qm.src.util import SetPositions


class OrcaDriver:
    # orca options
    options = {
        'method': 'B3LYP def2-SVP',
        'opt_geom': False,
        'calc_nmr': False,
        'calc_freq': False,
        'solvent': None,
        'geom_maxiter': 100,
        'n_jobs': 16,
        'mem_per_core': 2000,  # in MB
    }

    def __init__(
            self,
            orca_path: str = '/opt/orca-5.0.3/',
            openmpi_path: str = '/opt/openmpi-4.1.1',
            options: dict = None,
            show_progress: bool = False,
            tmp_root: str = '.orca_tmp'
    ):
        self.orca_path = orca_path
        self.omp_path = openmpi_path
        self.options = {**self.options, **(options or {})}
        self.tmp_root = tmp_root
        self.set_orca_env()

    def set_orca_env(self):
        if self.orca_path not in os.environ['PATH']:
            os.environ['PATH'] = f"{self.orca_path}:{os.environ.get('PATH')}"
            os.environ['LD_LIBRARY_PATH'] = f"{self.orca_path}:{os.environ.get('LD_LIBRARY_PATH', '')}"
        if self.omp_path not in os.environ['PATH']:
            os.environ['PATH'] = f"{self.omp_path}/bin:{os.environ['PATH']}"
            os.environ['LD_LIBRARY_PATH'] = f"{self.omp_path}/lib:{os.environ['LD_LIBRARY_PATH']}"

        # # intel MKL
        # os.environ['LD_LIBRARY_PATH'] = f"/usr/lib/x86_64-linux-gnu:{os.environ['LD_LIBRARY_PATH']}"

    def _create_tmp_dir(self):
        uid = str(uuid4())
        self.tmp_dir = f"{self.tmp_root}/{uid}"
        os.makedirs(self.tmp_dir, exist_ok=True)

    def start_routine(self):
        # save current dir
        self.saved_work_dir = os.getcwd()
        self._create_tmp_dir()
        os.chdir(self.tmp_dir)

    def end_routine(self):
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

    def _gen_input(self, conf: rdchem.Conformer):
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

        input_str += f"%pal nprocs {n_jobs} end\n"
        input_str += f"%maxcore {mem_per_core}\n"
        input_str += self._gen_xyz_block(conf)

        logging.debug(f"ORCA INPUT:\n{input_str}\n")
        return input_str

    def write_input(self, content: str):
        fname = 'task.inp'
        with open(fname, 'w') as f:
            f.write(content)
        return fname

    def write_output(self, content: str):
        fname = 'output'
        with open(fname, 'w') as f:
            f.write(content)
        return fname

    def run(self, conf: rdchem.Conformer):
        self.start_routine()

        input_str = self._gen_input(conf)
        output = self._run_orca(input_str)
        data = self._parse(output)

        self.end_routine()

        return data

    def _run_orca(self, input_str):
        fname = self.write_input(input_str)

        res = subprocess.run(
            [self.orca_path + '/orca', fname, '--use-hwthread-cpus'],
            capture_output=True,
            text=True
        )

        output, error = res.stdout, res.stderr

        if res.returncode != 0:
            raise Exception(f'Error running orca: \nINPUT {input_str}\nSTDOUT: {output}\nSTDERR:{error}')

        logging.debug(f"ORCA OUTPUT:\n{output}\n")

        return output

    def _parse(self, output: str):
        return cclib.io.ccread(io.StringIO(output))

    def single_point(self, conf):
        self.options['opt_geom'] = False
        data = self.run(conf)
        success = data.metadata['success']
        conf.SetBoolProp('success', success)

        if success:
            self.update_properties(conf, data)
            logging.info(f'Success: {success}')

    def geometry_optimization(self, conf):
        start = time()
        logging.info(f"Geometry optimization")
        logging.info(f"Method: {self.options['method']}")
        self.options['opt_geom'] = True
        data = self.run(conf)

        success = data.metadata['success']
        conf.SetBoolProp('success', success)

        if success:
            self.update_properties(conf, data)
            self.update_geometry(conf, data)
            logging.info(f'Success: {success}')
            logging.info(f'Num Iter: {len(data.atomcoords)}')
            logging.info(f'Elapsed Time: {time() - start:.1f}s')


    def update_geometry(self, conf, cclib_data):
        SetPositions(conf, cclib_data.atomcoords[-1])

    def update_properties(self, conf, cclib_data):
        conf.SetDoubleProp('energy', cclib_data.scfenergies[-1])