# cell2location
## install
https://github.com/BayraktarLab/cell2location
```bash
conda create -y -n cell2loc_env python=3.10

conda activate cell2loc_env
pip install cell2location[tutorials]
```
#如果要加入到jupyter-notebook中调用，请使用下述代码
```bash
conda activate cell2loc_env
python -m ipykernel install --user --name=cell2loc_env --display-name='Environment (cell2loc_env)'
```
