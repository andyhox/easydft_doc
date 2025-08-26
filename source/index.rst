Introduction
============

**easydft** 是一个用于材料科学和计算化学中 DFT 计算自动化的工具包。它通过自动任务生成、批量提交、结果收集与分析等功能简化了工作流程。
easydft 适用于所有用户，能够集成主流 DFT 软件，并基于 `pymatgen`_、 `jobflow-remote`_、 `ASE`_ 以及多个 `open-source VASP_script`_ 提供高效的 DFT 自动化平台。

.. _pymatgen: https://pymatgen.org/
.. _jobflow-remote: https://matgenix.github.io/jobflow-remote/
.. _ASE: https://wiki.fysik.dtu.dk/ase/
.. _open-source VASP_script: https://github.com/tamaswells/VASP_script/

主要功能包括：

- 自动生成和管理 DFT 计算任务
  
- 支持多种主流 DFT 软件包
  
- 便捷的批量作业提交与远程作业管理
  
- 自动结果收集与数据分析工具
  
- 易于扩展并集成到自定义工作流中

在使用 easydft 之前，请确保你已熟悉 `pymatgen <https://pymatgen.org/>`_ 的基本用法。
更重要的是，对 VASP 的使用有基本的认识，因为大部分代码读取 VASP 计算的输入和输出，请确保你了解相关内容。

如果你对工作流感兴趣，建议先阅读以下官方教程：

- `atomate2  <https://materialsproject.github.io/atomate2/>`_
- `jobflow <https://materialsproject.github.io/jobflow/tutorials.html>`_
- `jobflow-remote-quickstart <https://matgenix.github.io/jobflow-remote/user/quickstart.html>`_

欢迎联系老司机加入交流群

.. image:: /images/wechat.jpg
   :width: 200px
   :align: left

.. toctree::
   :maxdepth: 3
   :caption: 目录
   :hidden:
   
   api/modules
   docs/content
