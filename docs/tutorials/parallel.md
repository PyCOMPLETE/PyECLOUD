# Parallel execution

PyPARIS supports multiprocessing on one machine and MPI across processes that
may run on multiple machines. Choose process counts to match your available
resources and simulation configuration. The
[single-bunch tutorial](single-bunch.md) gives the corresponding launch commands.

## Local multiprocessing

The `PyPARIS.multiprocexec` launcher uses Python's multiprocessing support and
does not require an MPI installation. Its process count includes the master.

## MPI environments

MPI execution requires both an MPI runtime and the Python `mpi4py` package.
For a local Conda environment, one option is:

```bash
conda activate ecloud
conda install -c conda-forge mpi4py mpich
```

On a managed cluster, follow the site's instructions for loading its MPI module
and installing a matching `mpi4py`. Use the launcher belonging to that MPI
installation; mixing a Conda MPI runtime and a cluster launcher can cause
incompatibilities. Check what your active environment selects:

```bash
which mpiexec
python -c "from mpi4py import MPI; print(MPI.Get_library_version())"
```

The scheduler may require `srun` or another launch command instead of `mpiexec`.
PyPARIS MPI execution does not require parallel HDF5.

See the [mpi4py installation guide](https://mpi4py.readthedocs.io/en/stable/install.html)
for supported installation methods.

