import os
import shutil
from uuid import uuid4
import logging
from time import time
from rdkit.Chem import rdchem

from linux_qm.src.util import set_atom_positions


class _Driver:
    # default options
    options = {}
    jobname: str = 'task'

    def __init__(
            self,
            options: dict = None,
            show_progress: bool = False,
            tmp_root: str = None
    ):
        self.options = {**self.options, **(options or {})}
        self.tmp_root = tmp_root

    # Main methods
    def single_point(self, conf: rdchem.Conformer):
        """
        Single point property calculation.
        On success updates properties of the conformer.
        :param conf: rdkit Conformer
        :return: data: cclib parsed data
        """
        logging.info(f"Method: {self.options['method']}")
        self.options['opt_geom'] = False

        data = self.run(conf)

        success = data.metadata['success']
        conf.SetBoolProp('success', success)
        logging.info(f'Status: {success}')

        if success:
            # self.update_properties(conf, data)
            pass

        return data

    def geometry_optimization(self, conf: rdchem.Conformer):
        """
        Runs geometry optimization using self.options. On
        success - updates coordinates of conformer
        :param conf: rdkit Conformer
        :return: data: cclib parsed data
        """
        start = time()
        logging.info(f"Method: {self.options['method']}")
        self.options['opt_geom'] = True

        data = self.run(conf)

        success = data.metadata['success']
        conf.SetBoolProp('success', success)
        logging.info(f'Status: {success}')

        if success:
            # self.update_properties(conf, data)
            self.update_geometry(conf, data)
            logging.info(f'Num Iter: {len(data.atomcoords)}')
            logging.info(f'Elapsed Time: {time() - start:.1f}s')

        return data

    def run(self, conf: rdchem.Conformer):
        self._start_routine()

        input_str = self.gen_input(conf)
        output = self.run_cmd(input_str)
        data = self.parse(output)

        self._end_routine()

        return data

    def gen_input(self, conf: rdchem.Conformer):
        raise NotImplementedError

    def run_cmd(self, input_str):
        raise NotImplementedError

    def parse(self, output: str):
        raise NotImplementedError

    def update_geometry(self, conf, cclib_data):
        set_atom_positions(conf, cclib_data.atomcoords[-1])

    def _create_tmp_dir(self):
        uid = str(uuid4())
        self.tmp_dir = f"{self.tmp_root}/{uid}"
        os.makedirs(self.tmp_dir, exist_ok=True)

    def _start_routine(self):
        # save
        self.saved_work_dir = os.getcwd()

        self._create_tmp_dir()
        os.chdir(self.tmp_dir)

    def _end_routine(self):
        # restore
        os.chdir(self.saved_work_dir)

        shutil.rmtree(self.tmp_dir)


