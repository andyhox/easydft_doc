# Slab结构操作

## 指定hkl生成所有不重复暴露面表面

```python
from easydft.toolkit.structure_modify import SlabModify
from pymatgen.core import Structure

Si = Structure.from_file('./Si.vasp')

slabs = SlabModify.gen_all_slabs(Si, miller_index=(1,1,1), min_slab_size=4, min_vacuum_size=4)

for i, slab in enumerate(slabs):
    print(f"No.{i+1}")
    print(slab)
```

输出：

```text
No.1
Slab Summary (Si8)
Reduced Formula: Si
Miller index: (1, 1, 1)
Shift: 0.3750, Scale Factor: [[-1  1  0]
 [-1  0  1]
 [ 1  0  0]]
abc   :   3.758980   3.758980  24.649278
angles:  85.626993  85.626993  60.000000
Sites (8)
1 Si     0.015625     0.015625     0.703125
2 Si     0.359375     0.359375     0.671875
3 Si     0.390625     0.390625     0.578125
4 Si     0.734375     0.734375     0.546875
5 Si     0.765625     0.765625     0.453125
6 Si     0.109375     0.109375     0.421875
7 Si     0.140625     0.140625     0.328125
8 Si     0.484375     0.484375     0.296875
No.2
Slab Summary (Si8)
Reduced Formula: Si
Miller index: (1, 1, 1)
Shift: 0.8750, Scale Factor: [[-1  1  0]
 [-1  0  1]
 [ 1  0  0]]
abc   :   3.758980   3.758980  24.649278
angles:  85.626993  85.626993  60.000000
Sites (8)
1 Si     0.703125     0.703125     0.640625
2 Si     0.671875     0.671875     0.734375
3 Si     0.078125     0.078125     0.515625
4 Si     0.046875     0.046875     0.609375
5 Si     0.453125     0.453125     0.390625
6 Si     0.421875     0.421875     0.484375
7 Si     0.828125     0.828125     0.265625
8 Si     0.796875     0.796875     0.359375
```

## 生成唯一表面

```python
from easydft.toolkit.structure_modify import SlabModify
from pymatgen.core import Structure

Si = Structure.from_file('./Si.vasp')

slab = SlabModify.gen_slab(Si, miller_index=(1,1,1), min_slab_size=4, min_vacuum_size=4, shift=0.5)

print(slab)
```

通过调控`shift`得到指定的暴露面。

`shift`: The termination coordinate along the lattice c direction in fractional coordinates.