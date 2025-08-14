# 生成应变结构

`easydft`支持对结构进行**晶格常数**，**基矢角度**，**体积**三个维度批量生成应变结构: 

## 晶格常数应变

以bulk-Si为例，对应的模块为`easydft.toolkit.BulkModify`

```python
from easydft.toolkit.structure_modify import BulkModify
from pymatgen.core import Structure
import os

work_dir = os.getcwd()
Si = Structure.from_file(os.path.join(work_dir, 'Si.vasp'))
lattice_strain_list = [-0.1, -0.05, 0, 0.05, 0.1]
lattice_strain_structures = BulkModify.add_lattice_strain(
                                    structure=Si, 
                                    strain_list=lattice_strain_list, 
                                    mode='all'
                                    )
# save structures
for strain, structure in lattice_strain_structures.items(): 
    # show lattice
    print(f"应变值: {strain}")
    print(f"晶格常数: ")
    print(structure.lattice)
    filename = f"Si_{strain}.cif"
    structure.to(filename=os.path.join(work_dir, filename))
    print(f"应变值为{strain}的结构已保存为: {filename} \n")
```

Note: `mode`参数控制应变模式，`all`为a,b,c轴等比应变，单独`a`/`b`/`c`为单轴应变

输出: 

```text
应变值: -0.1
晶格常数: 
4.784400 0.000000 0.000000
-0.000000 4.784400 0.000000
0.000000 0.000000 4.784400
应变值为-0.1的结构已保存为: Si_-0.1.cif 

应变值: -0.05
晶格常数: 
5.050200 0.000000 0.000000
-0.000000 5.050200 0.000000
0.000000 0.000000 5.050200
应变值为-0.05的结构已保存为: Si_-0.05.cif 

应变值: 0
晶格常数: 
5.316000 0.000000 0.000000
-0.000000 5.316000 0.000000
0.000000 0.000000 5.316000
应变值为0的结构已保存为: Si_0.cif 

应变值: 0.05
晶格常数: 
5.581800 0.000000 0.000000
-0.000000 5.581800 0.000000
0.000000 0.000000 5.581800
应变值为0.05的结构已保存为: Si_0.05.cif 

应变值: 0.1
晶格常数: 
5.847600 0.000000 0.000000
-0.000000 5.847600 0.000000
0.000000 0.000000 5.847600
应变值为0.1的结构已保存为: Si_0.1.cif
```

## 基矢角度应变

```python
from easydft.toolkit.structure_modify import BulkModify
from pymatgen.core import Structure
import os

work_dir = os.getcwd()
Si = Structure.from_file(os.path.join(work_dir, 'Si.vasp'))
lattice_strain_list = [-0.1, -0.05, 0, 0.05, 0.1]
lattice_strain_structures = BulkModify.add_angle_strain(
                                    structure=Si, 
                                    strain_list=lattice_strain_list, 
                                    mode='alpha'
                                    )
# save structures
for strain, structure in lattice_strain_structures.items(): 
    # show angle
    print(f"应变值: {strain}")
    print(f"alpha: ")
    print(structure.lattice.alpha)
    filename = f"Si_{strain}.cif"
    structure.to(filename=os.path.join(work_dir, filename))
    print(f"应变值为{strain}的结构已保存为: {filename} \n")
```

Note: `mode`支持分别对$\alpha /\beta /\gamma $进行操作

## 体积应变

```python
from easydft.toolkit.structure_modify import BulkModify
from pymatgen.core import Structure
import os

work_dir = os.getcwd()

Si = Structure.from_file(os.path.join(work_dir, 'Si.vasp'))

lattice_strain_list = [-0.1, -0.05, 0, 0.05, 0.1]

lattice_strain_structures = BulkModify.add_volume_strain(
                                    structure=Si, 
                                    strain_list=lattice_strain_list, 
                                    )
# save structures
for strain, structure in lattice_strain_structures.items(): 
    # show volume and lattice
    print(f"应变值: {strain}")
    print(f"volume: ")
    print(structure.volume)
    print(f"lattice: ")
    print(structure.lattice)
    filename = f"Si_{strain}.cif"
    structure.to(filename=os.path.join(work_dir, filename))
    print(f"应变值为{strain}的结构已保存为: {filename} \n")
```

Note: `volume`应变是通过对abc等比缩放实现。

输出: 

```python
应变值: -0.1
volume: 
135.20645504639998
lattice: 
5.132542 0.000000 0.000000
-0.000000 5.132542 0.000000
0.000000 0.000000 5.132542
应变值为-0.1的结构已保存为: Si_-0.1.cif 

应变值: -0.05
volume: 
142.71792477120002
lattice: 
5.225881 0.000000 0.000000
-0.000000 5.225881 0.000000
0.000000 0.000000 5.225881
应变值为-0.05的结构已保存为: Si_-0.05.cif 

应变值: 0
volume: 
150.229394496
lattice: 
5.316000 0.000000 0.000000
-0.000000 5.316000 0.000000
0.000000 0.000000 5.316000
应变值为0的结构已保存为: Si_0.cif 

应变值: 0.05
volume: 
157.7408642208
lattice: 
5.403163 0.000000 0.000000
-0.000000 5.403163 0.000000
0.000000 0.000000 5.403163
应变值为0.05的结构已保存为: Si_0.05.cif 

应变值: 0.1
volume: 
165.25233394560001
lattice: 
5.487601 0.000000 0.000000
-0.000000 5.487601 0.000000
0.000000 0.000000 5.487601
应变值为0.1的结构已保存为: Si_0.1.cif
```

---

Slab模型的操作类似，不赘述~~~~