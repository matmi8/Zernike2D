import zepyros as zp
import pandas as pd
import numpy as np
from pathlib import Path
from docs import build_cli_zernike_invariants

np.seterr(divide='ignore', invalid='ignore')


class ZernikeInvariants:
    def __init__(self, surface_file, input_path, output_path, verso):
        self.surface = surface_file
        self.input_path = input_path
        self.output_path = output_path
        self.verso = verso

    def get_invariants(self):
        """
        Compute Zernike invariants and plot results
        """
        if int(self.verso) == 1:
            verso = 'up'
        elif int(self.verso) == -1:
            verso = 'down'
        else:
            raise ValueError('Wrong value for `verso`. It must be 1 or -1')
        
        surf = pd.read_csv(Path(self.input_path).joinpath(self.surface))

        coeff, disk_data, _, _ = zp.get_zernike(
            surf[['x', 'y', 'z', 'nx', 'ny', 'nz']], 
            radius=6.0, ndx=0, order=20, verso=int(self.verso))
        
        zp.plot_disk(disk_data, Path(self.output_path).joinpath(f'{str(self.surface)[:-4]}_{verso}_disk.png'))
        zp.plot_coeff(coeff, Path(self.output_path).joinpath(f'{str(self.surface)[:-4]}_{verso}_coeff.png'))

        np.savetxt(Path(self.output_path).joinpath(f'{str(self.surface)[:-4]}_{verso}_coeff.csv'), 
                   coeff, fmt="%f", delimiter="\n")


def main():
    args = build_cli_zernike_invariants()
    cz = ZernikeInvariants(args.surface, args.input, args.output, args.verso)
    cz.get_invariants()


if __name__ == "__main__":
    main()
