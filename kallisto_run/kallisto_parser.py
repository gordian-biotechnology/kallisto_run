from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
import sys

def get_parser():
    """Parser command line args."""
    from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter
    parser = ArgumentParser(description=__doc__,
                            formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("-i", "--input",
                        dest="input",
                        type=str,
                        required =True,
                        help="Path to directory with fastqs")
    parser.add_argument("-o", "--out",
                        dest="out",
                        type=str,
                        required =True,
                        help="Path to output directory.")
    parser.add_argument("-dr", "--dry_run",
                        dest="dry_run",
                        default=False,
                        action="store_true",
                        help="Don't run just print out commands.")
    parser.add_argument("-ref", "--ref_dir",
                        dest="ref_dir",
                        type=str,
                        help="Path to directory with references, should contain cdna.all.fa files. If not found will download by species, if .idx isn't also found it will build it.")
    parser.add_argument("-s", "--species",
                        dest="species",
                        default="mouse",
                        type=str,
                        help="Species of run, default is mouse.")
    parser.add_argument("-n", "--name",
                        dest="name",
                        default="",
                        type=str,
                        help="Name of to append to files.")
    parser.add_argument("-rr", "--rerun",
                        dest="rerun",
                        default=False,
                        action="store_true",)
    parser.add_argument("-v", "--version",
                        dest="version",
                        default='10xv3',
                        type=str,
                        help="Version of technology to provide to kallisto bus.")

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args=parser.parse_args()
    return args
