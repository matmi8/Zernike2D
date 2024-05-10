import subprocess
import argparse
from pathlib import Path
import pandas as pd
import numpy as np


class GetMolecularSurface:
    def __init__(self, pdb, input_path, output_path, all_atom):
        self.pdb = pdb
        self.input_path = input_path
        self.output_path = output_path
        self.all_atom = all_atom
        self.dms_output = None

    def get_dms_surface(self):
        """
        Calculates the molecular surface of a pdb using dms
        """
        dms_input = Path(self.input_path).joinpath(self.pdb)
        self.dms_output = Path(self.output_path).joinpath(f'{self.pdb[:-4]}.dms')
        
        if self.all_atom is True:
            run_dms = f'dms {dms_input} -n -a -o {self.dms_output}'
        else:
            run_dms = f'dms {dms_input} -n -o {self.dms_output}'    
        subprocess.call(run_dms, shell=True)

    def edit_dms_file(self):
        """
        Adds a character to column 14 in `.dms` file to prevent two adjacent columns 
        from interfering. This happens in some cases because of the x coordinates
        """
        with open(self.dms_output, 'r') as file:
            file_lines = [f'{x[:14]} {x[14:]}' for x in file.readlines()]

        with open(self.dms_output, 'w') as file:
            file.writelines(file_lines) 

    @staticmethod
    def range_char(start, stop):
        return (chr(n) for n in range(ord(start), ord(stop) + 1))
    
    def get_surface_csv(self):
        """
        Gets a `.csv` file with columns: 
        ["res", "x", "y", "z", "nx", "ny", "nz"]
        from a `.dms` file
        """
        colnames = [character for character in self.range_char("A", "K")]
        dms_surf = pd.read_csv(self.dms_output, names=colnames, header=None, delimiter=r'\s+')

        df = (dms_surf
            .replace(r'[^\w\s]|_', '', regex=True)
            .assign(res = np.where(
                dms_surf["B"].astype(str).str[-1].str.isnumeric(), 
                dms_surf["A"].astype(str) + "_" + 
                dms_surf["B"].astype(str) + "_" + 
                dms_surf["C"].astype(str), 
                dms_surf["A"].astype(str) + "_" + 
                (dms_surf["B"].astype(str).str
                 .extract(r'(\d+\.?\d*)([A-Za-z]*)', expand = True)
                 .agg('_'.join, axis=1)) + "_" + 
                dms_surf["C"])
            )
            .query('G.str[0] == "S"')
            .filter(items = ["res", "D", "E", "F", "I", "J", "K"])
            .set_axis(["res", "x", "y", "z", "nx", "ny", "nz"], axis = 1)
            )
        df.to_csv(f'{str(self.dms_output)[:-4]}.csv', index=False)


def build_cli():
    des = ("Given as input a `.pdb` file obtains a molecular surface using dms software")
    parser = argparse.ArgumentParser(description=des, add_help=False)
    req_grp = parser.add_argument_group(title="positional arguments")
    req_grp.add_argument("-pdb", "--pdb",
                         required=True,
                         help="pdb file; the pdb file name (ex: 1a1u_A.pdb)",
                         metavar="")
    req_grp.add_argument("-i", "--input",
                         required=True,
                         help="input path; folder with pdb files",
                         metavar="")
    req_grp.add_argument("-o", "--output",
                         required=True,
                         help="output path; destination folder of output files",
                         metavar="")
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('-a', 
                          action='store_true', 
                          help="boolean; use all atom in .pdb file to calculate dms surface")
    optional.add_argument("-h", "--help",
                          action="help",
                          help="show this help message and exit")
    return parser.parse_args()


def main():
    args = build_cli()
    ms = GetMolecularSurface(args.pdb, args.input, args.output, args.a)
    ms.get_dms_surface()
    ms.edit_dms_file()
    ms.get_surface_csv()


if __name__ == "__main__":
    main()