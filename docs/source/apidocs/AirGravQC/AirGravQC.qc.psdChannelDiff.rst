:py:mod:`AirGravQC.qc.psdChannelDiff`
=====================================

.. py:module:: AirGravQC.qc.psdChannelDiff

.. autodoc2-docstring:: AirGravQC.qc.psdChannelDiff
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`psdChannelDiff <AirGravQC.qc.psdChannelDiff.psdChannelDiff>`
     - .. autodoc2-docstring:: AirGravQC.qc.psdChannelDiff.psdChannelDiff
          :summary:
   * - :py:obj:`psdChannelGain <AirGravQC.qc.psdChannelDiff.psdChannelGain>`
     - .. autodoc2-docstring:: AirGravQC.qc.psdChannelDiff.psdChannelGain
          :summary:
   * - :py:obj:`_period_to_dist <AirGravQC.qc.psdChannelDiff._period_to_dist>`
     - .. autodoc2-docstring:: AirGravQC.qc.psdChannelDiff._period_to_dist
          :summary:
   * - :py:obj:`_dist_to_period <AirGravQC.qc.psdChannelDiff._dist_to_period>`
     - .. autodoc2-docstring:: AirGravQC.qc.psdChannelDiff._dist_to_period
          :summary:
   * - :py:obj:`_time_frequency <AirGravQC.qc.psdChannelDiff._time_frequency>`
     - .. autodoc2-docstring:: AirGravQC.qc.psdChannelDiff._time_frequency
          :summary:
   * - :py:obj:`_mean_line_speed <AirGravQC.qc.psdChannelDiff._mean_line_speed>`
     - .. autodoc2-docstring:: AirGravQC.qc.psdChannelDiff._mean_line_speed
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`groupName <AirGravQC.qc.psdChannelDiff.groupName>`
     - .. autodoc2-docstring:: AirGravQC.qc.psdChannelDiff.groupName
          :summary:

API
~~~

.. py:data:: groupName
   :canonical: AirGravQC.qc.psdChannelDiff.groupName
   :value: None

   .. autodoc2-docstring:: AirGravQC.qc.psdChannelDiff.groupName

.. py:function:: psdChannelDiff(whizzFile, channel1, channel2, flightLines=[])
   :canonical: AirGravQC.qc.psdChannelDiff.psdChannelDiff

   .. autodoc2-docstring:: AirGravQC.qc.psdChannelDiff.psdChannelDiff

.. py:function:: psdChannelGain(whizzFile, rawchan, filchan, flightLines=[], nominalPeriod=0.0, verbose=False)
   :canonical: AirGravQC.qc.psdChannelDiff.psdChannelGain

   .. autodoc2-docstring:: AirGravQC.qc.psdChannelDiff.psdChannelGain

.. py:function:: _period_to_dist(p)
   :canonical: AirGravQC.qc.psdChannelDiff._period_to_dist

   .. autodoc2-docstring:: AirGravQC.qc.psdChannelDiff._period_to_dist

.. py:function:: _dist_to_period(x)
   :canonical: AirGravQC.qc.psdChannelDiff._dist_to_period

   .. autodoc2-docstring:: AirGravQC.qc.psdChannelDiff._dist_to_period

.. py:function:: _time_frequency(group)
   :canonical: AirGravQC.qc.psdChannelDiff._time_frequency

   .. autodoc2-docstring:: AirGravQC.qc.psdChannelDiff._time_frequency

.. py:function:: _mean_line_speed(group, line)
   :canonical: AirGravQC.qc.psdChannelDiff._mean_line_speed

   .. autodoc2-docstring:: AirGravQC.qc.psdChannelDiff._mean_line_speed
