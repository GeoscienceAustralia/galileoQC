:py:mod:`AirGravQC.qc.qualityAnalysis`
======================================

.. py:module:: AirGravQC.qc.qualityAnalysis

.. autodoc2-docstring:: AirGravQC.qc.qualityAnalysis
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`calcDrift <AirGravQC.qc.qualityAnalysis.calcDrift>`
     - .. autodoc2-docstring:: AirGravQC.qc.qualityAnalysis.calcDrift
          :summary:
   * - :py:obj:`checkVertAcc <AirGravQC.qc.qualityAnalysis.checkVertAcc>`
     - .. autodoc2-docstring:: AirGravQC.qc.qualityAnalysis.checkVertAcc
          :summary:
   * - :py:obj:`checkVertAccStats <AirGravQC.qc.qualityAnalysis.checkVertAccStats>`
     - .. autodoc2-docstring:: AirGravQC.qc.qualityAnalysis.checkVertAccStats
          :summary:
   * - :py:obj:`checkStatcor <AirGravQC.qc.qualityAnalysis.checkStatcor>`
     - .. autodoc2-docstring:: AirGravQC.qc.qualityAnalysis.checkStatcor
          :summary:
   * - :py:obj:`lineStats <AirGravQC.qc.qualityAnalysis.lineStats>`
     - .. autodoc2-docstring:: AirGravQC.qc.qualityAnalysis.lineStats
          :summary:
   * - :py:obj:`statsChannelDiff <AirGravQC.qc.qualityAnalysis.statsChannelDiff>`
     - .. autodoc2-docstring:: AirGravQC.qc.qualityAnalysis.statsChannelDiff
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`groupName <AirGravQC.qc.qualityAnalysis.groupName>`
     - .. autodoc2-docstring:: AirGravQC.qc.qualityAnalysis.groupName
          :summary:

API
~~~

.. py:data:: groupName
   :canonical: AirGravQC.qc.qualityAnalysis.groupName
   :value: None

   .. autodoc2-docstring:: AirGravQC.qc.qualityAnalysis.groupName

.. py:function:: calcDrift(whizzFile, time, gradient)
   :canonical: AirGravQC.qc.qualityAnalysis.calcDrift

   .. autodoc2-docstring:: AirGravQC.qc.qualityAnalysis.calcDrift

.. py:function:: checkVertAcc(whizzFile, vertvelocity)
   :canonical: AirGravQC.qc.qualityAnalysis.checkVertAcc

   .. autodoc2-docstring:: AirGravQC.qc.qualityAnalysis.checkVertAcc

.. py:function:: checkVertAccStats(whizzFile)
   :canonical: AirGravQC.qc.qualityAnalysis.checkVertAccStats

   .. autodoc2-docstring:: AirGravQC.qc.qualityAnalysis.checkVertAccStats

.. py:function:: checkStatcor(whizzFile, statcor, flight='')
   :canonical: AirGravQC.qc.qualityAnalysis.checkStatcor

   .. autodoc2-docstring:: AirGravQC.qc.qualityAnalysis.checkStatcor

.. py:function:: lineStats(whizzFile, channel)
   :canonical: AirGravQC.qc.qualityAnalysis.lineStats

   .. autodoc2-docstring:: AirGravQC.qc.qualityAnalysis.lineStats

.. py:function:: statsChannelDiff(whizzFile, channel1, channel2, flightLines=[])
   :canonical: AirGravQC.qc.qualityAnalysis.statsChannelDiff

   .. autodoc2-docstring:: AirGravQC.qc.qualityAnalysis.statsChannelDiff
