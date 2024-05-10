import zepyros as zp
import pandas as pd
import numpy as np
from pathlib import Path
import argparse


class ComputeZernike:
    def __init__(self, surface_file, input_path, output_path, verso):
        self.surface = surface_file
        self.input_path = input_path
        self.output_path = output_path
        self.verso = verso

    def get_invariants(self):
        surf = pd.read_csv(Path(self.input_path).joinpath(self.surface))

        coeff, disk_data, _, _ = zp.get_zernike(
            surf[['x', 'y', 'z', 'nx', 'ny', 'nz']], 
            radius=6.0, ndx=0, order=20, verso=int(self.verso))
        if int(self.verso) == 1:
            verso = 'up'
        elif int(self.verso) == -1:
            verso = 'down'
        else:
            raise ValueError('Wrong value for `verso`. It must be 1 or -1')
        zp.plot_disk(disk_data, Path(self.output_path).joinpath(f'{str(self.surface)[:-4]}_{verso}_disk.png'))
        zp.plot_coeff(coeff, Path(self.output_path).joinpath(f'{str(self.surface)[:-4]}_{verso}_coeff.png'))

        np.savetxt(Path(self.output_path).joinpath(f'{str(self.surface)[:-4]}_{verso}_coeff.csv'), 
                   coeff, fmt="%f", delimiter="\n")


def build_cli():
    des = ("Given as input a `.pdb` file obtains a molecular surface using dms software")
    parser = argparse.ArgumentParser(description=des, add_help=False)
    req_grp = parser.add_argument_group(title="positional arguments")
    req_grp.add_argument("-sf", "--surface",
                         required=True,
                         help="surface file; the csv file name (ex: 1a1u_A.csv)",
                         metavar="")
    req_grp.add_argument("-i", "--input",
                         required=True,
                         help="input path; folder with surface file",
                         metavar="")
    req_grp.add_argument("-o", "--output",
                         required=True,
                         help="output path; destination folder of output files",
                         metavar="")
    req_grp.add_argument("-v", "--verso",
                         required=True,
                         help="verso; direction of the patch with respect to the z axis; 1 positive, -1 negative",
                         metavar="")
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-h", "--help",
                          action="help",
                          help="show this help message and exit")
    return parser.parse_args()


def main():
    args = build_cli()
    cz = ComputeZernike(args.surface, args.input, args.output, args.verso)
    cz.get_invariants()


if __name__ == "__main__":
    main()