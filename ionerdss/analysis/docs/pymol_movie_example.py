"""
Example: export a PyMOL movie from a NERDSS PDB trajectory.

Run from a simulation directory, or update `SIM_DIR` below:

    pymol -cq /absolute/path/to/pymol_movie_example.py
"""

from pathlib import Path

from ionerdss.analysis.visualization import export_pymol_pdb_movie


# Point this at the simulation root that contains:
# - `parms.inp`
# - `PDB/1.pdb`, `PDB/2.pdb`, ...
# - optionally exported `*_system.json` or `*detailed_summary*.txt`
SIM_DIR = Path("/Users/yueying/Downloads/microtubule/35")


def main() -> None:
    frames_dir = export_pymol_pdb_movie(
        base_dir=SIM_DIR,
        pdb_dir=SIM_DIR / "PDB",
        add_timestamp_overlay=True,
        timestamp_font_size=96,
        movie_fps=10,
        out_dir_name="frames",
    )


if __name__ == "__main__":
    main()
