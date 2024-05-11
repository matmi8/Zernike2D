import argparse


def build_cli_molecular_surface():
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


def build_cli_zernike_invariants():
    des = (f"Given a surface of which coordinates and versors of each point are defined, "
           f"calculate the 2D Zernike invariants")
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


def build_cli_binding_propensity():
    des = ("Given two surfaces to compare, compute the binding propensity")
    parser = argparse.ArgumentParser(description=des, add_help=False)
    req_grp = parser.add_argument_group(title="positional arguments")
    req_grp.add_argument("-sf1", "--surface1",
                         required=True,
                         help="surface1 file; full path of the surface1 csv file",
                         metavar="")
    req_grp.add_argument("-sf2", "--surface2",
                         required=True,
                         help="surface2 file; full path of the surface2 csv file",
                         metavar="")
    req_grp.add_argument("-o", "--output",
                         required=True,
                         help="output path; destination folder of output files",
                         metavar="")
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-h", "--help",
                          action="help",
                          help="show this help message and exit")
    return parser.parse_args()