# 异质结建模操作

异质结建模关键参数除了`film`和`substrate`结构以及对应的`miller index`之外，还需要读取指定的`film`和`substrate`的暴露面，即`termination`参数。

因此在生成异质结之间，需要先选择合适的`termination`参数，可以通过下列代码把所有对应的`termination`列出来：

```python
from easydft.toolkit.structure_modify import HeterojunctionModify
from pymatgen.core import Structure

film_structure = Structure.from_file('./LiF.cif')
substrate_structure = Structure.from_file('./Si.vasp')

hetergen = HeterojunctionModify(
    film_structure=film_structure,
    substrate_structure=substrate_structure,
    film_miller_index=(0,0,1),
    substrate_miller_index=(1,1,1))

terminations = hetergen.terminations

print(terminations)
```

输出结果返回一个列表：
```text
[('LiF_P4/mmm_2', 'Si_P4/mmm_1')]
```

选择对应的`termination`进行建模.

```python
from easydft.toolkit.structure_modify import HeterojunctionModify
from pymatgen.core import Structure

film_structure = Structure.from_file('./LiF.cif')
substrate_structure = Structure.from_file('./Si.vasp')

hetergen = HeterojunctionModify(
    film_structure=film_structure,
    substrate_structure=substrate_structure,
    film_miller_index=(0,0,1),
    substrate_miller_index=(0,0,1))

terminations = hetergen.terminations

heter_structures = hetergen.gen_heter_interface(termination=terminations[0])
i = 0
for i, heter in enumerate(heter_structures):
    i += 1
print(f"total structures: {i}")
```

输出：total structures: 446

为什么会生成这么多结构呢？

实际上代码内部进行晶面匹配的时候会设置相应的cutoff：

- max_area_ratio_tol=0.09<p></p>

- max_area=400<p></p>

- max_length_tol=0.03<p></p>

- max_angle_tol=0.01<p></p>

代码会输出所有面积小于400的异质结结构。

当然，在实际使用中，不需要每个结构都输出，446个结构是按照由简到繁，原子数由少到多排列的，一般情况下我们只需要用第一个结构即可。

```python
from easydft.toolkit.structure_modify import HeterojunctionModify
from pymatgen.core import Structure

film_structure = Structure.from_file('./LiF.cif')
substrate_structure = Structure.from_file('./Si.vasp')

hetergen = HeterojunctionModify(
    film_structure=film_structure,
    substrate_structure=substrate_structure,
    film_miller_index=(0,0,1),
    substrate_miller_index=(0,0,1))

terminations = hetergen.terminations

heter_structures = hetergen.gen_heter_interface(termination=terminations[0])

heter_lists = list(heter_structures)
print(heter_lists[0])
```

输出：
```text
Full Formula (Li16 Si20 F16)
Reduced Formula: Li4Si5F4
abc   :   8.405334   8.405334  28.028714
angles:  90.000000  90.000000  90.000000
pbc   :       True       True       True
Sites (52)
  #  SP       a      b         c    bulk_equivalent  bulk_wyckoff    interface_label
---  ----  ----  -----  --------  -----------------  --------------  -----------------
  0  Li+   0      0     0.60523                   0  a               film
  1  Li+   0      0.5   0.60523                   0  a               film
  2  Li+   0.25   0.25  0.60523                   0  a               film
  3  Li+   0.5    0     0.60523                   0  a               film
  4  Li+   0.25   0.75  0.60523                   0  a               film
  5  Li+   0.5    0.5   0.60523                   0  a               film
  6  Li+   0.75   0.25  0.60523                   0  a               film
  7  Li+   0.75   0.75  0.60523                   0  a               film
  8  Li+   0.25   0     0.532386                  0  a               film
  9  Li+   0.25   0.5   0.532386                  0  a               film
 10  Li+   0.5    0.25  0.532386                  0  a               film
 11  Li+   0.75   0     0.532386                  0  a               film
 12  Li+   0.5    0.75  0.532386                  0  a               film
 13  Li+   0.75   0.5   0.532386                  0  a               film
 14  Li+   0      0.25  0.532386                  0  a               film
 15  Li+   1      0.75  0.532386                  0  a               film
 16  Si    0.7    0.6   0.461031                  0  a               substrate
 17  Si    0.1    0.8   0.461031                  0  a               substrate
 18  Si    0.5   -0     0.461031                  0  a               substrate
 19  Si    0.9    0.2   0.461031                  0  a               substrate
 20  Si    0.3    0.4   0.461031                  0  a               substrate
 21  Si    0.2    0.6   0.318784                  0  a               substrate
 22  Si    0.6    0.8   0.318784                  0  a               substrate
 23  Si    0     -0     0.318784                  0  a               substrate
 24  Si    0.4    0.2   0.318784                  0  a               substrate
 25  Si    0.8    0.4   0.318784                  0  a               substrate
 26  Si    0.6    0.3   0.366199                  0  a               substrate
 27  Si    0      0.5   0.366199                  0  a               substrate
 28  Si    0.4    0.7   0.366199                  0  a               substrate
 29  Si    0.8    0.9   0.366199                  0  a               substrate
 30  Si    0.2    0.1   0.366199                  0  a               substrate
 31  Si    0.5    0.5   0.413615                  0  a               substrate
 32  Si    0.9    0.7   0.413615                  0  a               substrate
 33  Si    0.3    0.9   0.413615                  0  a               substrate
 34  Si    0.7    0.1   0.413615                  0  a               substrate
 35  Si    0.1    0.3   0.413615                  0  a               substrate
 36  F-    0      0     0.532386                  4  b               film
 37  F-    0      0.5   0.532386                  4  b               film
 38  F-    0.25   0.25  0.532386                  4  b               film
 39  F-    0.5    0     0.532386                  4  b               film
 40  F-    0.25   0.75  0.532386                  4  b               film
 41  F-    0.5    0.5   0.532386                  4  b               film
 42  F-    0.75   0.25  0.532386                  4  b               film
 43  F-    0.75   0.75  0.532386                  4  b               film
 44  F-    0.25   0     0.60523                   4  b               film
 45  F-    0.25   0.5   0.60523                   4  b               film
 46  F-    0.5    0.25  0.60523                   4  b               film
 47  F-    0.75   0     0.60523                   4  b               film
 48  F-    0.5    0.75  0.60523                   4  b               film
 49  F-    0.75   0.5   0.60523                   4  b               film
 50  F-    0      0.25  0.60523                   4  b               film
 51  F-    1      0.75  0.60523                   4  b               film
 ```