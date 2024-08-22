:py:mod:`AirGravQC.qc.checkRepeatLines`
=======================================

.. py:module:: AirGravQC.qc.checkRepeatLines

.. autodoc2-docstring:: AirGravQC.qc.checkRepeatLines
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`checkRepeatLines <AirGravQC.qc.checkRepeatLines.checkRepeatLines>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkRepeatLines.checkRepeatLines
          :summary:
   * - :py:obj:`reportTrackDirection <AirGravQC.qc.checkRepeatLines.reportTrackDirection>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkRepeatLines.reportTrackDirection
          :summary:
   * - :py:obj:`_xBaseInterpolant <AirGravQC.qc.checkRepeatLines._xBaseInterpolant>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkRepeatLines._xBaseInterpolant
          :summary:
   * - :py:obj:`_plotRepeatAnalysis <AirGravQC.qc.checkRepeatLines._plotRepeatAnalysis>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkRepeatLines._plotRepeatAnalysis
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`groupName <AirGravQC.qc.checkRepeatLines.groupName>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkRepeatLines.groupName
          :summary:

API
~~~

.. py:data:: groupName
   :canonical: AirGravQC.qc.checkRepeatLines.groupName
   :value: None

   .. autodoc2-docstring:: AirGravQC.qc.checkRepeatLines.groupName

.. py:function:: checkRepeatLines(whizzFiles, channel, repeatLines, x='', z='', xOffset=True, verbose=False)
   :canonical: AirGravQC.qc.checkRepeatLines.checkRepeatLines

   .. autodoc2-docstring:: AirGravQC.qc.checkRepeatLines.checkRepeatLines

.. py:function:: reportTrackDirection(surveygroup, line, east='', north='')
   :canonical: AirGravQC.qc.checkRepeatLines.reportTrackDirection

   .. autodoc2-docstring:: AirGravQC.qc.checkRepeatLines.reportTrackDirection

.. py:function:: _xBaseInterpolant(whizzFiles, channel, repeatLines, x='', z='', verbose=False)
   :canonical: AirGravQC.qc.checkRepeatLines._xBaseInterpolant

   .. autodoc2-docstring:: AirGravQC.qc.checkRepeatLines._xBaseInterpolant

.. py:function:: _plotRepeatAnalysis(xBase, xOffset, nLines, xData, yData, zData, channel, flightLines, z, chan_z_units, chan_y_label, chan_y_units, baseLine=-1.0)
   :canonical: AirGravQC.qc.checkRepeatLines._plotRepeatAnalysis

   .. autodoc2-docstring:: AirGravQC.qc.checkRepeatLines._plotRepeatAnalysis
