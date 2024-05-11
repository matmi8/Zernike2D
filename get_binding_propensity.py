import zepyros as zp
import pandas as pd
import numpy as np
from pathlib import Path
import tqdm
from docs import build_cli_binding_propensity

np.seterr(divide='ignore', invalid='ignore')


class BindingPropensity:
    def __init__(self, surface_file1, surface_file2, output_path):
        self.surface1 = surface_file1
        self.surface2 = surface_file2
        self.output_path = output_path

    @staticmethod
    def get_all_invariants(surface, verso):
        surf = pd.read_csv(surface)
        nrow = len(surf.index)
        coeff_array = np.zeros((nrow, 121))

        for i in tqdm.tqdm(range(nrow)):
            coeff, _, _, _ = zp.get_zernike(
                surf[['x', 'y', 'z', 'nx', 'ny', 'nz']], 6.0, i, 20, int(verso)
            )
            coeff_array[i, :] = coeff

        return surf, coeff_array
    
    @staticmethod
    def get_bp_mask(df, point):
        mask = np.where((
            df['x'].isin(point[0]) & 
            df['y'].isin(point[1]) & 
            df['z'].isin(point[2])), 1, 0)
        return mask

    
    def get_binding_propensity(self):
        df_surf1, mat_coeff1 = self.get_all_invariants(self.surface1, verso=1)
        df_surf2, mat_coeff2 = self.get_all_invariants(self.surface2, verso=-1)

        patch1, patch2 = zp.contact_points(
            df_surf1[['x', 'y', 'z']].to_numpy(), 
            df_surf2[['x', 'y', 'z']].to_numpy(), 
            1
        )

        dist = np.sqrt(((mat_coeff1 - mat_coeff2[:, None])**2).sum(2))

        color_1 = np.min(dist, axis=0)
        color_2 = np.min(dist, axis=1)
        bp1 = self.get_bp_mask(df_surf1, patch1)
        bp2 = self.get_bp_mask(df_surf2, patch2)

        df_bp1 = pd.DataFrame({
            'x': df_surf1['x'],
            'y': df_surf1['y'],
            'z': df_surf1['z'],
            'color': color_1,
            'bp': bp1
        })

        df_bp2 = pd.DataFrame({
            'x': df_surf2['x'],
            'y': df_surf2['y'],
            'z': df_surf2['z'],
            'color': color_2,
            'bp': bp2
        })

        file_name1 = Path(self.surface1).stem
        file_name2 = Path(self.surface2).stem

        df_bp1.to_csv(Path(self.output_path).joinpath(f'{file_name1}_bp.csv'))
        df_bp2.to_csv(Path(self.output_path).joinpath(f'{file_name2}_bp.csv'))


def main():
    args = build_cli_binding_propensity()
    bp = BindingPropensity(args.surface1, args.surface2, args.output)
    bp.get_binding_propensity()


if __name__ == "__main__":
    main()
