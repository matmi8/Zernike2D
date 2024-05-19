import zepyros as zp
import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
from pathlib import Path
import tqdm
from docs import build_cli_binding_propensity

np.seterr(divide='ignore', invalid='ignore')


class BindingPropensity:
    def __init__(self, surface_file1, surface_file2, output_path):
        self.surface1 = pd.read_csv(surface_file1)
        self.surface2 = pd.read_csv(surface_file2)
        self.file_name1 = Path(surface_file1).stem
        self.file_name2 = Path(surface_file2).stem
        self.output_path = output_path

    @staticmethod
    def get_all_invariants(surface, verso):
        """
        Compute Zernike invariants for each point of a surface
        Parameters
        ----------
        - `surface`: full path of the surface .csv file
        - `verso`: +1 or -1
        """
        surf = surface
        nrow = len(surf.index)
        sampled = range(0, nrow, 10)
        coeff_array = np.zeros((len(sampled), 121))

        for i, ndx in enumerate(tqdm.tqdm(sampled)):
            coeff, _, _, _ = zp.get_zernike(
                surf[['x', 'y', 'z', 'nx', 'ny', 'nz']], 6.0, ndx, 20, int(verso)
            )
            coeff_array[i, :] = coeff

        return surf.iloc[sampled], coeff_array
    
    @staticmethod
    def get_bp_mask(df, points):
        """
        Filter patch points indexes from reduced surface dataframe
        
        Paramters
        ---------
        - `df`: filtered surface dataframe
        - `points`: patch points

        Return
        ------
        - `mask`: indexes of the patch point in filtered surface dataframe
        """
        mask = np.where((
            df['x'].isin(points[:, 0]) & 
            df['y'].isin(points[:, 1]) & 
            df['z'].isin(points[:, 2])), 1, 0)
        return mask

    
    def get_binding_propensity(self):
        """
        Compute binding propensity
        """
        print(f'compute Zernike invariants for {self.file_name1}')
        df_surf1, mat_coeff1 = self.get_all_invariants(self.surface1, verso=1)
        print(f'compute Zernike invariants for {self.file_name2}')
        df_surf2, mat_coeff2 = self.get_all_invariants(self.surface2, verso=-1)

        patch1, patch2 = zp.contact_points(
            df_surf1[['x', 'y', 'z']].to_numpy(), 
            df_surf2[['x', 'y', 'z']].to_numpy(), 
            3
        )

        dist = cdist(mat_coeff1, mat_coeff2)
        
        color_1 = np.min(dist, axis=1)
        color_2 = np.min(dist, axis=0)
        bp1 = self.get_bp_mask(df_surf1, patch1)
        bp2 = self.get_bp_mask(df_surf2, patch2)

        df_bp1 = pd.DataFrame({
            'x': df_surf1['x'].to_numpy(),
            'y': df_surf1['y'].to_numpy(),
            'z': df_surf1['z'].to_numpy(),
            'color': color_1,
            'bp': bp1
        })

        df_bp2 = pd.DataFrame({
            'x': df_surf2['x'].to_numpy(),
            'y': df_surf2['y'].to_numpy(),
            'z': df_surf2['z'].to_numpy(),
            'color': color_2,
            'bp': bp2
        })

        df_bp1.to_csv(Path(self.output_path).joinpath(f'{self.file_name1}_bp.csv'), index=False)
        df_bp2.to_csv(Path(self.output_path).joinpath(f'{self.file_name2}_bp.csv'), index=False)


def main():
    args = build_cli_binding_propensity()
    bp = BindingPropensity(args.surface1, args.surface2, args.output)
    bp.get_binding_propensity()


if __name__ == "__main__":
    main()
