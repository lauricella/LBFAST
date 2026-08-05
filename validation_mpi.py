"""Shared MPI command-line support for LBFAST validation drivers."""

from pathlib import Path


def add_mpi_arguments(parser):
    parser.add_argument("--mpi-procs", type=int, default=1,
                        help="number of MPI ranks (default: 1, sequential run)")
    parser.add_argument("--decomposition", type=int, nargs=3,
                        metavar=("PX", "PY", "PZ"),
                        help="MPI Cartesian decomposition; product must equal --mpi-procs")
    parser.add_argument("--launcher", choices=("mpirun", "srun"), default="mpirun",
                        help="MPI process launcher")


def validate_mpi_arguments(args):
    if args.mpi_procs < 1:
        raise SystemExit("--mpi-procs must be positive")
    mpi_target = "mpi" in args.make_target.lower()
    if args.mpi_procs == 1:
        if args.decomposition not in (None, [1, 1, 1]):
            raise SystemExit("a sequential run only accepts decomposition 1 1 1")
        if mpi_target:
            raise SystemExit("an MPI Make target requires --mpi-procs greater than one")
        return
    if not mpi_target:
        raise SystemExit("--mpi-procs greater than one requires an MPI Make target")
    if args.decomposition is None:
        raise SystemExit("MPI runs require --decomposition PX PY PZ")
    if any(value < 1 for value in args.decomposition):
        raise SystemExit("MPI decomposition values must be positive")
    if args.decomposition[0] * args.decomposition[1] * args.decomposition[2] != args.mpi_procs:
        raise SystemExit("PX*PY*PZ must equal --mpi-procs")


def validate_local_tiles(input_path, decomposition, tile=(8, 8, 8)):
    if decomposition is None:
        return
    text = Path(input_path).read_text()
    import re
    sizes = []
    for name in ("lx", "ly", "lz"):
        match = re.search(r"(?m)^\s*{}\s*=\s*([0-9]+)".format(name), text)
        if not match:
            raise SystemExit("could not read {} from {}".format(name, input_path))
        sizes.append(int(match.group(1)))
    local = []
    for size, processes, tile_size, axis in zip(sizes, decomposition, tile, "xyz"):
        if size % processes:
            raise SystemExit("global {} size {} is not divisible by {} MPI ranks".format(
                axis, size, processes))
        local_size = size // processes
        if local_size % tile_size:
            raise SystemExit("local {} size {} is not a multiple of tile {}".format(
                axis, local_size, tile_size))
        local.append(local_size)
    return tuple(local)


def solver_command(args, binary, input_name):
    if args.mpi_procs == 1:
        return [str(binary), str(input_name)]
    px, py, pz = args.decomposition
    count_flag = "-np" if args.launcher == "mpirun" else "-n"
    return [args.launcher, count_flag, str(args.mpi_procs), str(binary),
            str(px), str(py), str(pz), str(input_name)]
