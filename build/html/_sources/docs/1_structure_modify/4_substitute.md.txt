# 掺杂结构操作

## 替换掺杂

`easydft`支持快速生成指定数量的**不重复**掺杂结构

```python
from easydft.toolkit.structure_modify import gen_substitute_structures
from pymatgen.core import Structure

Si = Structure.from_file('./Si.vasp')

supercell = Si.make_supercell([2,2,2])

doped_structures = gen_substitute_structures(
    structure=supercell,
    substitute_element='Si',
    dopant='Ge',
    dopant_num=5,
    max_structures=5,
    save=False
)

for i, struct in enumerate(doped_structures):
    print(f"No.{i}")
    print(struct)
```

输出：

```text
No.0
Full Formula (Si59 Ge5)
Reduced Formula: Si59Ge5
abc   :  10.632000  10.632000  10.632000
angles:  90.000000  90.000000  90.000000
pbc   :       True       True       True
Sites (64)
......
......
```

## 空位结构

生成空位结构，只需要在上述代码基础上，`dopant`参数传入`vac`即可。

```python
from easydft.toolkit.structure_modify import gen_substitute_structures
from pymatgen.core import Structure

Si = Structure.from_file('./Si.vasp')

supercell = Si.make_supercell([2,2,2])

doped_structures = gen_substitute_structures(
    structure=supercell,
    substitute_element='Si',
    dopant='vac',
    dopant_num=5,
    max_structures=5,
    save=False
)

for i, struct in enumerate(doped_structures):
    print(f"No.{i}")
    print(struct)
```

输出：
```text
No.0
Full Formula (Si59)
Reduced Formula: Si
abc   :  10.632000  10.632000  10.632000
angles:  90.000000  90.000000  90.000000
pbc   :       True       True       True
Sites (59)
......
......
```

Note：生成的结构为**不重复**结构，而非**对称性去重**结构，更加适用于多原子掺杂或对称性较低的使用场景。