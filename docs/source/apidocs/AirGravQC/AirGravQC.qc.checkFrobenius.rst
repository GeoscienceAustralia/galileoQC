:py:mod:`AirGravQC.qc.checkFrobenius`
=====================================

.. py:module:: AirGravQC.qc.checkFrobenius

.. autodoc2-docstring:: AirGravQC.qc.checkFrobenius
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`checkFrobenius <AirGravQC.qc.checkFrobenius.checkFrobenius>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkFrobenius.checkFrobenius
          :summary:
   * - :py:obj:`_FTGeigen <AirGravQC.qc.checkFrobenius._FTGeigen>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkFrobenius._FTGeigen
          :summary:
   * - :py:obj:`_FTGTransform <AirGravQC.qc.checkFrobenius._FTGTransform>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkFrobenius._FTGTransform
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`groupName <AirGravQC.qc.checkFrobenius.groupName>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkFrobenius.groupName
          :summary:

API
~~~

.. py:data:: groupName
   :canonical: AirGravQC.qc.checkFrobenius.groupName
   :value: None

   .. autodoc2-docstring:: AirGravQC.qc.checkFrobenius.groupName

.. py:function:: checkFrobenius(whizzFile, lines=[], il1='Inline1_raw', il2='Inline2_raw', il3='Inline3_raw', cr1='Cross1_raw', cr2='Cross2_raw', cr3='Cross3_raw', noiselimit=30.0, verbose=True, plot_flag=False)
   :canonical: AirGravQC.qc.checkFrobenius.checkFrobenius

   .. autodoc2-docstring:: AirGravQC.qc.checkFrobenius.checkFrobenius

.. py:function:: _FTGeigen(Txx, Txy, Txz, Tyy, Tyz, Tzz, line='', noiselimit=30.0, plot_flag=False)
   :canonical: AirGravQC.qc.checkFrobenius._FTGeigen

   .. autodoc2-docstring:: AirGravQC.qc.checkFrobenius._FTGeigen

.. py:function:: _FTGTransform(il1, il2, il3, cr1, cr2, cr3)
   :canonical: AirGravQC.qc.checkFrobenius._FTGTransform

   .. autodoc2-docstring:: AirGravQC.qc.checkFrobenius._FTGTransform
