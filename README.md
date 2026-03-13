# sci-scripts

A collection of scripts for quantum chemistry workflows: HPC job submission, Gaussian input/output processing, and FCclasses3 integration.

---

## Repository structure

```
sci-scripts/
├── Cluster/        # HPC job management and monitoring
├── Gaussian/       # Gaussian input/output utilities
└── FCclasses/      # FCclasses3 utilities
```

---

## Cluster

Scripts for submitting and monitoring jobs on a PBS-based HPC cluster. Designed for a two-hop connection (local machine → university machine → cluster), with SSH host support throughout.

### `gaussian-qsub-umu.sh`

PBS submission script for running multiple Gaussian calculations in parallel. Self-submitting: run it directly and it calls `qsub` on itself.

After each calculation, converts `.chk` → `.fchk` (via `formchk`) and compresses to `.fchk.gz`. Progress is logged to `progress.log`.

```bash
./gaussian-qsub-umu.sh [options]
```

| Flag | Description | Default |
|------|-------------|---------|
| `-N` | PBS job name | `numders-fw` |
| `-n` | Node name | `nv81` |
| `-p` | Cores per node (ppn) | `32` |
| `-f` | File list of `.com` inputs | `FILES` |
| `-c` | Cores per Gaussian calculation | `8` |
| `-m` | Max parallel calculations | `4` |
| `-h` | Show help | |

The `FILES` file should contain one `.com` filename per line.

---

### `gaussian-qsub-umu-remote.py`

Submit Gaussian calculations **from a local machine** to the cluster. Handles file transfer and job submission in one step.

Reads the `.com` files listed in `FILES`, uploads them along with `gaussian-qsub-umu.sh` to the remote directory, and submits the job.

```bash
gaussian-qsub-umu-remote.py -H HOST [options]
```

| Flag | Description | Default |
|------|-------------|---------|
| `-H` / `--host` | SSH host **(required)** | |
| `-d` / `--dir` | Remote working directory | `$HOME/remote-calculations-tmp` |
| `-f` / `--files` | Local file list | `FILES` |
| `-N` | PBS job name | `numders-fw` |
| `-n` | Cluster node | `nv81` |
| `-p` | Cores per node | `32` |
| `-c` | Cores per Gaussian calculation | `8` |
| `-m` | Max parallel calculations | `4` |

After job completion, retrieve results with `folder-sync.py --down`.

**Requires** `gaussian-qsub-umu.sh` in the same directory.

---

### `folder-sync.py`

Sync a directory between local machine and remote host via `rsync`. Only transfers files that have changed.

```bash
folder-sync.py -H HOST -ld LOCAL_DIR -hd HOST_DIR (--up | --down) [options]
```

| Flag | Description |
|------|-------------|
| `-H` / `--host` | SSH host **(required)** |
| `-ld` / `--local-dir` | Local directory **(required)** |
| `-hd` / `--host-dir` | Remote directory **(required)** |
| `--up` / `--upload` | Upload: local → remote |
| `--down` / `--download` | Download: remote → local |
| `--dry-run` | Preview without transferring |
| `--delete` | Delete destination files absent from source |

```bash
# Upload
folder-sync.py -H cluster -ld ./my_calcs -hd /scratch/jose/my_calcs --up

# Download results
folder-sync.py -H cluster -ld ./my_calcs -hd /scratch/jose/my_calcs --down

# Preview
folder-sync.py -H cluster -ld ./my_calcs -hd /scratch/jose/my_calcs --up --dry-run
```

---

### `qstat-info.py`

Monitor PBS job progress. Shows job state, walltime, host, calculation list, and a progress bar (running / completed / failed) parsed from `progress.log`.

```bash
qstat-info.py                        # all jobs for $USER
qstat-info.py -u alice               # jobs for another user
qstat-info.py 12345.scc              # specific job ID
qstat-info.py --host cluster         # query via SSH
qstat-info.py --watch 60             # auto-refresh every 60 s
```

---

### `pbsnodes-info.py`

Show PBS node status: free / busy / down, cores available per node, and running job IDs.

```bash
pbsnodes-info.py                     # all nodes
pbsnodes-info.py -f                  # free nodes only
pbsnodes-info.py --host cluster      # query via SSH
pbsnodes-info.py -f --watch 30       # auto-refresh every 30 s
```

---

## Gaussian

### `extract-optimized-geom-gaussian.py`

Extract the final optimized geometry from a Gaussian `.log` file and write it in `.xyz` format.

```bash
extract-optimized-geom-gaussian.py file.log
```

---

### `gen-findiff-numder-displaced-input-files.py`

Generate displaced-geometry Gaussian input files for numerical differentiation of transition dipole moments (TDMs). Supports Cartesian (`-der X`) and normal mode (`-der Q`) displacements.

```bash
gen-findiff-numder-displaced-input-files.py -f data.fcc -der X -disp 0.001 [options]
```

Key flags: `-f`, `-der X|Q`, `-disp`, `-outdir`, `-gauhead`. See `-h` for full list.

---

## FCclasses

### `gen-findiff-numder-displaced-input-files.py`

Same as the Gaussian version but with a different default header (`B3LYP/6-31G(d)`). Generates displaced geometry inputs for FCclasses3 workflows.

---

### `gen-d2num-findiff-dipfiles.py`

Compute second-derivative TDMs (d²μ/dxᵢdxⱼ) from displaced `.fchk` files using finite differences. Supports velocity and length gauges, 5-point differentiation, and optional symmetrization.

```bash
gen-d2num-findiff-dipfiles.py -f ref.fchk -disp_dir displaced/ -gauge vel [options]
```

Key flags: `-f`, `-f0`, `-gauge`, `-delta`, `-disp_dir`, `-extract`. See `-h` for full list.

---

## SSH setup

All scripts with `--host` run commands on the remote via SSH. For passwordless access from the university machine to the cluster, add a key:

```bash
ssh-keygen -t ed25519 -f ~/.ssh/id_ed25519_cdc   # on the university machine
ssh-copy-id -i ~/.ssh/id_ed25519_cdc user@cluster
```

Then set `~/.ssh/config`:

```
Host cluster
    HostName cluster.hostname
    User your_username
    IdentityFile ~/.ssh/id_ed25519_cdc
```
