# 简单能带&态密度

`easydft`大部分功能为计算数据后处理所写。对于普通的能带态密度后处理，考虑到`pymatgen`，`vaspvis`以及`pyprocar`等包已经可以直接读取`VASP`输出文件直接出图。`easydft`目前主要用于提取出可绘图的csv数据供二次直接绘图。从数据储存来说，相比于保存`VASP`的输出文件以供随时后处理，csv数据文件大小更友好，适合大量长期保存。

## 态密度

从`easydft`中导入`DosAnalyzer`模块：

```python
from easydft.toolkit.electronic_properties import DosAnalyzer
```

只需要传入`vasprun.xml`文件或完整的态密度计算文件夹，同时指定数据文件保存的目录即可。

```python
from easydft.toolkit.electronic_properties import DosAnalyzer

analyzer = DosAnalyzer(work_dir='./Si', output_path='./')
# analyzer = DosAnalyzer(vasprun_file='./Si/vasprun.xml', output_path='./')
```

### 元素态密度

```python
from easydft.toolkit.electronic_properties import DosAnalyzer

analyzer = DosAnalyzer(work_dir='./Si', output_path='./')

analyzer.parse_element_dos()
```

当前路径下会生成`element_dos.csv`文件：

```text
Energy,tdos,Si
-13.809566,0.0,0.0
-13.796866,0.0,0.0
-13.784165999999999,0.0,0.0
-13.771466,0.0,0.0
-13.758766,0.0,0.0
-13.746065999999999,0.0,0.0
-13.733366,0.0,0.0
-13.720666,0.0,0.0
-13.707866,0.0,0.0
-13.695166,0.0,0.0
-13.682466,0.0,0.0
-13.669766,0.0,0.0
-13.657066,0.0,0.0
......
```

### 轨道态密度

```python
from easydft.toolkit.electronic_properties import DosAnalyzer

analyzer = DosAnalyzer(work_dir='./Si', output_path='./')

analyzer.parse_element_spd_dos()
```
当前路径下会生成`element_spd_dos.csv`文件。


### 亚轨道态密度

处理亚轨道态密度要求`LORBIT >= 11`

```python
from easydft.toolkit.electronic_properties import DosAnalyzer

analyzer = DosAnalyzer(work_dir='./Si', output_path='./')

analyzer.parse_element_spd_dos()
```

当前路径下会生成`orbital_spd_dos.csv`文件。

----

后续可根据绘图习惯自由选择软件进行绘图~~~

## 能带

用法与态密度处理一致：

```python
from easydft.toolkit.electronic_properties import BandStructureAnalyzer

analyzer = BandStructureAnalyzer(vasprun_file='./Si11/11/vasprun.xml', output_path='./')

analyzer.parse_normal_bandstructure()
```

当前路径下会生成`band.csv`文件：

```text
band_index,k_idx,kpoint,kpoint_distance,energy
0,0,GAMMA,0.0,-12.653966
0,1,"0.0263,0.0000,0.0263",0.03721614712256227,-12.641166
0,2,"0.0526,0.0000,0.0526",0.07443229424512454,-12.602666
0,3,"0.0789,0.0000,0.0789",0.11164844136768681,-12.538566
0,4,"0.1053,0.0000,0.1053",0.14886458849024908,-12.448865999999999
0,5,"0.1316,0.0000,0.1316",0.18608073561281135,-12.333966
0,6,"0.1579,0.0000,0.1579",0.22329688273537363,-12.193566
```

绘图时选择`kpoints_distance`为x，`energy`为y即可。

band_index表示自旋方向，`easydft`自动识别`ISPIN`参数，当`ISPIN=1`时，band_index无意义；当`ISPIN=2`时，band_index为0表示自旋向上，为1表示自旋向下。