import glob
import os
import shutil
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor

import argparse

parser = argparse.ArgumentParser(description='Execute notebooks for inclusion '
                                 'in documentation')
parser.add_argument('--kernel', metavar='kernel', type=str, nargs=1,
                    help='Kernel to execute the notebooks with. If none'
                    'will use the current kernel')

args = parser.parse_args()

excluded_notebooks = ['CB.ipynb',
                      'CH82 -- optimization.ipynb',
                      'DC-pyps comparison.ipynb',
                      'Example_MLL_Fit_GlyR_4patches.ipynb']
full_path = os.path.realpath(__file__)
mydir, a = os.path.split(full_path)

input_dir = os.path.join(mydir, '..', 'exploration')
output_dir = os.path.join(mydir, '..', 'documentation', 'source', 'notebooks')
notebooks = glob.glob(os.path.join(input_dir, '*.ipynb'))
if args.kernel:
    kernel_name = args.kernel[0]
else:
    kernel_name = ''
ep = ExecutePreprocessor(timeout=600, kernel_name=kernel_name)
for notebook in notebooks:
    _, filename = os.path.split(notebook)
    executed_nbname = os.path.join(output_dir, filename)
    if filename not in excluded_notebooks:
        print("Executing {}".format(notebook))
        with open(notebook) as f:
            nb = nbformat.read(f, as_version=4)
            ep.preprocess(nb, {'metadata': {'path': input_dir}})
        print("Saving to {}".format(executed_nbname))
        nbformat.write(nb, executed_nbname)
    else:
        print("{} excluded from execution".format(notebook))
        shutil.copyfile(notebook, executed_nbname)
