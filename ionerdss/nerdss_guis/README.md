# Install lib

```bash
conda create --name myenv python=3.10
conda activate myenv
conda install pip
pip install numpy pandas matplotlib seaborn tqdm PyQt6 pyqtgraph PyOpenGL biopython
pip install pyqt6-tools
```

From the repository root:

```bash
python -m ionerdss.nerdss_guis.nerdss
```

# For developer

generate the .ui file using `designer`

convert the .ui to its module in `gen/` (e.g. the main window):
```bash
pyuic6 -x gen/mainwindow.ui -o gen/mainwindow.py
```
