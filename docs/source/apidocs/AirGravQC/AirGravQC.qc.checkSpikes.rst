:py:mod:`AirGravQC.qc.checkSpikes`
==================================

.. py:module:: AirGravQC.qc.checkSpikes

.. autodoc2-docstring:: AirGravQC.qc.checkSpikes
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`checkSpikes <AirGravQC.qc.checkSpikes.checkSpikes>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkSpikes.checkSpikes
          :summary:
   * - :py:obj:`spike_test2 <AirGravQC.qc.checkSpikes.spike_test2>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkSpikes.spike_test2
          :summary:
   * - :py:obj:`spike_test1 <AirGravQC.qc.checkSpikes.spike_test1>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkSpikes.spike_test1
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`groupName <AirGravQC.qc.checkSpikes.groupName>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkSpikes.groupName
          :summary:

API
~~~

.. py:data:: groupName
   :canonical: AirGravQC.qc.checkSpikes.groupName
   :value: None

   .. autodoc2-docstring:: AirGravQC.qc.checkSpikes.groupName

.. py:function:: checkSpikes(whizzFile, channels=[], lines=[], numStd=8.0, window=0, verbose=False)
   :canonical: AirGravQC.qc.checkSpikes.checkSpikes

   .. autodoc2-docstring:: AirGravQC.qc.checkSpikes.checkSpikes

.. py:function:: spike_test2(data, threshold, report_start, window=0, verbose=False)
   :canonical: AirGravQC.qc.checkSpikes.spike_test2

   .. autodoc2-docstring:: AirGravQC.qc.checkSpikes.spike_test2

.. py:function:: spike_test1(data, threshold, report_start, verbose=False)
   :canonical: AirGravQC.qc.checkSpikes.spike_test1

   .. autodoc2-docstring:: AirGravQC.qc.checkSpikes.spike_test1
