:py:mod:`AirGravQC.whizzPlots.whizzPlot`
========================================

.. py:module:: AirGravQC.whizzPlots.whizzPlot

.. autodoc2-docstring:: AirGravQC.whizzPlots.whizzPlot
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`plotWsLineChannel <AirGravQC.whizzPlots.whizzPlot.plotWsLineChannel>`
     - .. autodoc2-docstring:: AirGravQC.whizzPlots.whizzPlot.plotWsLineChannel
          :summary:
   * - :py:obj:`_subplotCompare <AirGravQC.whizzPlots.whizzPlot._subplotCompare>`
     - .. autodoc2-docstring:: AirGravQC.whizzPlots.whizzPlot._subplotCompare
          :summary:
   * - :py:obj:`_get_data <AirGravQC.whizzPlots.whizzPlot._get_data>`
     - .. autodoc2-docstring:: AirGravQC.whizzPlots.whizzPlot._get_data
          :summary:
   * - :py:obj:`plot_xcohere <AirGravQC.whizzPlots.whizzPlot.plot_xcohere>`
     - .. autodoc2-docstring:: AirGravQC.whizzPlots.whizzPlot.plot_xcohere
          :summary:
   * - :py:obj:`make_plot_title <AirGravQC.whizzPlots.whizzPlot.make_plot_title>`
     - .. autodoc2-docstring:: AirGravQC.whizzPlots.whizzPlot.make_plot_title
          :summary:
   * - :py:obj:`specificLinesMap <AirGravQC.whizzPlots.whizzPlot.specificLinesMap>`
     - .. autodoc2-docstring:: AirGravQC.whizzPlots.whizzPlot.specificLinesMap
          :summary:
   * - :py:obj:`statusMap <AirGravQC.whizzPlots.whizzPlot.statusMap>`
     - .. autodoc2-docstring:: AirGravQC.whizzPlots.whizzPlot.statusMap
          :summary:
   * - :py:obj:`_plot_speed <AirGravQC.whizzPlots.whizzPlot._plot_speed>`
     - .. autodoc2-docstring:: AirGravQC.whizzPlots.whizzPlot._plot_speed
          :summary:
   * - :py:obj:`_plotcheckSafeClearance <AirGravQC.whizzPlots.whizzPlot._plotcheckSafeClearance>`
     - .. autodoc2-docstring:: AirGravQC.whizzPlots.whizzPlot._plotcheckSafeClearance
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`groupName <AirGravQC.whizzPlots.whizzPlot.groupName>`
     - .. autodoc2-docstring:: AirGravQC.whizzPlots.whizzPlot.groupName
          :summary:

API
~~~

.. py:data:: groupName
   :canonical: AirGravQC.whizzPlots.whizzPlot.groupName
   :value: None

   .. autodoc2-docstring:: AirGravQC.whizzPlots.whizzPlot.groupName

.. py:function:: plotWsLineChannel(whizzFile1, flightLine1, channel1, whizzFile2, flightLine2, channel2, x1='', x2='', y1='', y2='', h1='', h2='', plotTitle='', xOffset=False)
   :canonical: AirGravQC.whizzPlots.whizzPlot.plotWsLineChannel

   .. autodoc2-docstring:: AirGravQC.whizzPlots.whizzPlot.plotWsLineChannel

.. py:function:: _subplotCompare(ax, x1, y1, x2, y2, c1, c2, xLabel, yLabel, plotTitle)
   :canonical: AirGravQC.whizzPlots.whizzPlot._subplotCompare

   .. autodoc2-docstring:: AirGravQC.whizzPlots.whizzPlot._subplotCompare

.. py:function:: _get_data(whizzFile, flightLine, channel, x='', y='', h='')
   :canonical: AirGravQC.whizzPlots.whizzPlot._get_data

   .. autodoc2-docstring:: AirGravQC.whizzPlots.whizzPlot._get_data

.. py:function:: plot_xcohere(whizzFile, flightLine, xchannel, ychannel)
   :canonical: AirGravQC.whizzPlots.whizzPlot.plot_xcohere

   .. autodoc2-docstring:: AirGravQC.whizzPlots.whizzPlot.plot_xcohere

.. py:function:: make_plot_title(group)
   :canonical: AirGravQC.whizzPlots.whizzPlot.make_plot_title

   .. autodoc2-docstring:: AirGravQC.whizzPlots.whizzPlot.make_plot_title

.. py:function:: specificLinesMap(whizzFile, lines, easting='', northing='')
   :canonical: AirGravQC.whizzPlots.whizzPlot.specificLinesMap

   .. autodoc2-docstring:: AirGravQC.whizzPlots.whizzPlot.specificLinesMap

.. py:function:: statusMap(planFile='', planEast='', planNorth='', plotTitle='')
   :canonical: AirGravQC.whizzPlots.whizzPlot.statusMap

   .. autodoc2-docstring:: AirGravQC.whizzPlots.whizzPlot.statusMap

.. py:function:: _plot_speed(t, t_label, speed, min_speed=54, max_speed=66, plot_title='')
   :canonical: AirGravQC.whizzPlots.whizzPlot._plot_speed

   .. autodoc2-docstring:: AirGravQC.whizzPlots.whizzPlot._plot_speed

.. py:function:: _plotcheckSafeClearance(projName, line, distance, clearance_chan='', altitude_chan='', terrain_chan='', alt=[], dtm=[], clearance=[])
   :canonical: AirGravQC.whizzPlots.whizzPlot._plotcheckSafeClearance

   .. autodoc2-docstring:: AirGravQC.whizzPlots.whizzPlot._plotcheckSafeClearance
