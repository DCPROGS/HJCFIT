# This script creates a virtual environment with all necessary packages
# to install HJCFIT and run its tests.

# Set things up to create the dcprogs environment:
# Note installation defaults to $HOME/.conda/envs, but it has to be created in the $WORK
# filesystem for the compute nodes to be able to access it.
# For that, we need to manually create a folder in $HOME:
mkdir -pv  $WORK/.conda/envs

# In other to make "activate" find the environment, we need to set this env var:
export CONDA_ENVS_PATH=$WORK/.conda/envs

# Create dcprogs virtual env in $WORK folder. Note it can only have Anaconda-compute's python 2
# since anaconda+python3 is not available in the MOM nodes.
conda create -p $CONDA_ENVS_PATH/dcprogs python numpy scipy six pip -y

# Behave needs to be installed with pip from inside the virtual env, so activate it
source activate dcprogs

# To pip install behave with --user in the right place, we need to set up another env var.
# Note $CONDA_ENV_PATH points to $CONDA_ENVS_PATH/<virtual_env_name> after activating.
export PYTHONUSERBASE=$CONDA_ENV_PATH
pip install --user behave

# Get out from the virtual environment
source deactivate
