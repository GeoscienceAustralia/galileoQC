:py:mod:`AirGravQC.qc.checkSpeeds`
==================================

.. py:module:: AirGravQC.qc.checkSpeeds

.. autodoc2-docstring:: AirGravQC.qc.checkSpeeds
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`checkSpeeds <AirGravQC.qc.checkSpeeds.checkSpeeds>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkSpeeds.checkSpeeds
          :summary:
   * - :py:obj:`_reportSpeeds <AirGravQC.qc.checkSpeeds._reportSpeeds>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkSpeeds._reportSpeeds
          :summary:
   * - :py:obj:`_get_speeddata <AirGravQC.qc.checkSpeeds._get_speeddata>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkSpeeds._get_speeddata
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`groupName <AirGravQC.qc.checkSpeeds.groupName>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkSpeeds.groupName
          :summary:

API
~~~

.. py:data:: groupName
   :canonical: AirGravQC.qc.checkSpeeds.groupName
   :value: None

   .. autodoc2-docstring:: AirGravQC.qc.checkSpeeds.groupName

.. py:function:: checkSpeeds(whizzFile, xChannel='', yChannel='', tChannel='', vel_north='', vel_east='', nominalSpeed=60.0, maxDuration=0.0, maxDistance=0.0, allowance=0.1, allowed_range=[], minSafeSpeed=42.0, known='', plot_flag=False)
   :canonical: AirGravQC.qc.checkSpeeds.checkSpeeds

   .. autodoc2-docstring:: AirGravQC.qc.checkSpeeds.checkSpeeds

.. py:function:: _reportSpeeds(group, maxDuration=0.0, maxDistance=0.0, xChannel='X', yChannel='Y', tChannel='time', vel_north='', vel_east='', nominalSpeed=60.0, allowance=0.1, allowed_range=[], minSafeSpeed=42.0, title_str='', known='', plot_flag=False)
   :canonical: AirGravQC.qc.checkSpeeds._reportSpeeds

   .. autodoc2-docstring:: AirGravQC.qc.checkSpeeds._reportSpeeds

.. py:function:: _get_speeddata(line_group, xChannel, yChannel, tChannel, vel_north, vel_east)
   :canonical: AirGravQC.qc.checkSpeeds._get_speeddata

   .. autodoc2-docstring:: AirGravQC.qc.checkSpeeds._get_speeddata
