:py:mod:`AirGravQC.qc.checkHighFreq`
====================================

.. py:module:: AirGravQC.qc.checkHighFreq

.. autodoc2-docstring:: AirGravQC.qc.checkHighFreq
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`checkHighFreq <AirGravQC.qc.checkHighFreq.checkHighFreq>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkHighFreq.checkHighFreq
          :summary:
   * - :py:obj:`_subplot_hiF_analysis <AirGravQC.qc.checkHighFreq._subplot_hiF_analysis>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkHighFreq._subplot_hiF_analysis
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`groupName <AirGravQC.qc.checkHighFreq.groupName>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkHighFreq.groupName
          :summary:

API
~~~

.. py:data:: groupName
   :canonical: AirGravQC.qc.checkHighFreq.groupName
   :value: None

   .. autodoc2-docstring:: AirGravQC.qc.checkHighFreq.groupName

.. py:function:: checkHighFreq(whizzFile, lines=[], noiseLimit=50, channels=[], cutoffs=[0.15, 3.6], tChannel='', vertaccel='', vertvelocity='', vertdispl='', verbose=False, plot_flag=False)
   :canonical: AirGravQC.qc.checkHighFreq.checkHighFreq

   .. autodoc2-docstring:: AirGravQC.qc.checkHighFreq.checkHighFreq

.. py:function:: _subplot_hiF_analysis(fig, num_subplots, plotIdx, plotTitle, x1, y1, x2=np.array([]), y2=np.array([]), bounds=[])
   :canonical: AirGravQC.qc.checkHighFreq._subplot_hiF_analysis

   .. autodoc2-docstring:: AirGravQC.qc.checkHighFreq._subplot_hiF_analysis
