# build container with latest R and sctk

* result: sctk_v2.18.0_R_v4.5.3.sif (latest sctk and R)
* and where possible: latest versions of other R tools we need but with patches
* build and test in a sandbox before finalising a sif, incorporating changes into the def file for rebuilding
* the reason we use a base R docker - so we can run various R tools later from this container and use BioC
* for dedicated py tools later - build a dedicated Py container 
* when sctk_v2.18.0_R_v4.5.3.sif is done - test locally on PBMC 1K

prebuild the local R v4.5.3 docker first in bash console to use as a local image - just a few mins
```bash
export APPTAINER_CACHEDIR="/global/scratch/$USER/test_cache"
export APPTAINER_TMPDIR="/global/scratch/$USER/test_tmp"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
apptainer pull local_rocker.sif docker://rocker/r-ver:4.5.3
rm -rf "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"
```
check the OS specs of rocker as it affects the python definitions in the def script
```bash
apptainer exec local_rocker.sif cat /etc/os-release

PRETTY_NAME="Ubuntu 24.04.4 LTS"
NAME="Ubuntu"
VERSION_ID="24.04"
VERSION="24.04.4 LTS (Noble Numbat)"
VERSION_CODENAME=noble
ID=ubuntu
ID_LIKE=debian
...
```
* for a robust def script: specify the deps used in the SCTK_runQC.R (v2.18.0)
* if you use sandbox building - keep track of changes for reverse engineering into the final def script
```bash
# reverse engineered from the finalised sandbox sif, testing in progress - build and use at own risk

cat << 'EOF' > sctk_delphine_22april.def
Bootstrap: localimage
From: ./local_rocker.sif

%post
# 1. OS dependencies & Python Environment
apt-get update && apt-get install -y \
tree curl wget ca-certificates git \
libcurl4-openssl-dev libssl-dev libxml2-dev zlib1g-dev \
libfontconfig1-dev libharfbuzz-dev libfribidi-dev \
libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
libicu-dev libssh2-1-dev libgit2-dev \
libhdf5-dev libbz2-dev liblzma-dev libglpk-dev \
make g++ gcc cmake \
python3-pip python3-dev libpython3-dev python3-full

# System-level R configuration
export MAKEFLAGS="-j8"
{
echo "MAKEFLAGS=-j8"
echo "PKG_SYSREQS=true"
echo "PAK_BUILD_MAX_CONCURRENCY=8"
echo "RETICULATE_PYTHON=/opt/sctk_env/bin/python"
echo "BASILISK_USE_EXTERNAL=1"
echo "RETICULATE_AUTOCREATE=FALSE"
} >> /usr/local/lib/R/etc/Renviron

# Isolated Python Env
python3 -m venv /opt/sctk_env
/opt/sctk_env/bin/pip install --no-cache-dir --upgrade pip setuptools wheel
/opt/sctk_env/bin/pip install --no-cache-dir "setuptools<70" "numpy<2" scrublet scanpy loompy anndata

# Install SCTK via pak
Rscript -e "install.packages('pak', repos = 'https://cloud.r-project.org')"
Rscript -e "pak::repo_add(BiocManager = '3.22')"
Rscript -e "pak::pkg_install(c('optparse', 'yaml', 'igraph', 'Rtsne', 'spam', 'MCMCprecision', 'celda', 'BiocParallel', 'GSEABase', 'SummarizedExperiment', 'compbiomed/singleCellTK@v2.18.0'))"

# custom sctk 2.18.0 script
wget -O /usr/local/bin/sctk_mito3_newnet_23april.R https://raw.githubusercontent.com/GITC2025/integration_benchmark/main/sctk_mito3_newnet_23april.R
chmod +x /usr/local/bin/sctk_mito3_newnet_23april.R

# --- GLOBAL RPROFILE LOCKDOWN ---
cat << 'RPROF' > /usr/local/lib/R/etc/Rprofile.site
local({
options(download.file.method = 'unsupported')
setHook(packageEvent("reticulate", "onLoad"), function(...) {
if (exists("pyenv_bootstrap_unix", envir = asNamespace("reticulate"))) {
message(">>> [HPC LOCKDOWN] Neutralizing reticulate bootstrap")
unlockBinding("pyenv_bootstrap_unix", asNamespace("reticulate"))
assignInNamespace("pyenv_bootstrap_unix", function(...) { return(NULL) }, ns="reticulate")
lockBinding("pyenv_bootstrap_unix", asNamespace("reticulate"))
}
})
})
RPROF

%environment
export LC_ALL=C
export LANG=C.UTF-8
export PATH=/usr/local/bin:$PATH
export RETICULATE_PYTHON=/opt/sctk_env/bin/python
export BASILISK_USE_EXTERNAL=1
EOF
```

## build container from local image
* temp files must be pointed to local node SSD for speed - otherwise it takes forever to write to scratch space
* local node SSD extraction time to temp about 10mins, on scratch space temp - many hours
```bash
#!/bin/bash
#SBATCH --job-name=build_sctk_v2
#SBATCH --time=04:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=128G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=delphine.girona@gmail.com
#SBATCH --output=/global/scratch/%u/sctk_build_%j.out
#SBATCH --error=/global/scratch/%u/sctk_build_%j.err

set -euo pipefail

# resource limit adjustments
ulimit -n $(ulimit -Hn)
ulimit -u $(ulimit -Hu)

module load apptainer/1.4.5

# local node SSD for temp files to maximize R compilation speed!!
export APPTAINER_TMPDIR="$SLURM_TMPDIR"

# cache on scratch to avoid filling up node RAM/disk
export APPTAINER_CACHEDIR="/global/scratch/$USER/apptainer_cache"

mkdir -p "$APPTAINER_CACHEDIR"

cd /global/scratch/$USER

apptainer build --fakeroot sctk_v2.18.0_R_v4.5.3.sif sctk_delphine_22april.def

echo -e "\n[ALERT] Build complete."
```
check disk usage in the $SLURM_TMPDIR - it should be growing rapidly up to 4GB+
```bash
srun --jobid=8223292 --overlap bash -c 'du -sh $SLURM_TMPDIR'
```
apptainer sif files can build even when there are errors
* on completion, parse the logs to check for specific errors
* try it out on a small sample in console 

## output
sif build should complete < 15 mins (writing to local node was critical for speed!)
```bash
apptainer exec sctk_v2.18.0_R_v4.5.3.sif Rscript -e "cat(sprintf('R Version: %s\nSCTK Version: %s\nBioC Version: %s\n', R.version.string, packageVersion('singleCellTK'), BiocManager::version()))"

R Version: R version 4.5.3 (2026-03-11)
SCTK Version: 2.18.0
BioC Version: 3.22
```

* specs, debug, shimmy, test and comparison [here](https://github.com/GITC2025/integration_benchmark/blob/main/stck_container_test.md)
* when using the container, CPU architecture may need to be [specified](https://github.com/GITC2025/integration_benchmark/blob/main/stck_container_test.md#cpu-architecture-specs)

## sandbox notes
entering an apptainer sandbox with HPC internet connection to install deps
```bash
apptainer shell --fakeroot --writable \
--env http_proxy=$http_proxy \
--env https_proxy=$https_proxy \
sctk_sandbox/
```
when building from sandbox
```bash
apptainer build sctk_v2.18.0_R_v4.5.3.sif sctk_sandbox/
```
