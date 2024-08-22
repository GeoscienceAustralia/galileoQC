:py:mod:`AirGravQC.qc.checkIntersection`
========================================

.. py:module:: AirGravQC.qc.checkIntersection

.. autodoc2-docstring:: AirGravQC.qc.checkIntersection
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`checkIntersection <AirGravQC.qc.checkIntersection.checkIntersection>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkIntersection.checkIntersection
          :summary:
   * - :py:obj:`_ccw <AirGravQC.qc.checkIntersection._ccw>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkIntersection._ccw
          :summary:
   * - :py:obj:`_intersect <AirGravQC.qc.checkIntersection._intersect>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkIntersection._intersect
          :summary:
   * - :py:obj:`_intersection_height <AirGravQC.qc.checkIntersection._intersection_height>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkIntersection._intersection_height
          :summary:
   * - :py:obj:`_lines_cross <AirGravQC.qc.checkIntersection._lines_cross>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkIntersection._lines_cross
          :summary:
   * - :py:obj:`_calc_bearing <AirGravQC.qc.checkIntersection._calc_bearing>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkIntersection._calc_bearing
          :summary:
   * - :py:obj:`_mean_1std <AirGravQC.qc.checkIntersection._mean_1std>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkIntersection._mean_1std
          :summary:
   * - :py:obj:`controls_lessthan_1000 <AirGravQC.qc.checkIntersection.controls_lessthan_1000>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkIntersection.controls_lessthan_1000
          :summary:
   * - :py:obj:`controls_nineteen <AirGravQC.qc.checkIntersection.controls_nineteen>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkIntersection.controls_nineteen
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`groupName <AirGravQC.qc.checkIntersection.groupName>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkIntersection.groupName
          :summary:

API
~~~

.. py:data:: groupName
   :canonical: AirGravQC.qc.checkIntersection.groupName
   :value: None

   .. autodoc2-docstring:: AirGravQC.qc.checkIntersection.groupName

.. py:function:: checkIntersection(whizzFile, controls=[], travs=[], xChannel='', yChannel='', zChannel='', max_allowed_deltaZ=10.0, plot_flag=False)
   :canonical: AirGravQC.qc.checkIntersection.checkIntersection

   .. autodoc2-docstring:: AirGravQC.qc.checkIntersection.checkIntersection

.. py:function:: _ccw(x1, y1, x2, y2, x3, y3)
   :canonical: AirGravQC.qc.checkIntersection._ccw

   .. autodoc2-docstring:: AirGravQC.qc.checkIntersection._ccw

.. py:function:: _intersect(cx1, cy1, cx2, cy2, tx1, ty1, tx2, ty2)
   :canonical: AirGravQC.qc.checkIntersection._intersect

   .. autodoc2-docstring:: AirGravQC.qc.checkIntersection._intersect

.. py:function:: _intersection_height(x_trav, y_trav, z_trav, x_ctrl, y_ctrl, z_ctrl, bearingc)
   :canonical: AirGravQC.qc.checkIntersection._intersection_height

   .. autodoc2-docstring:: AirGravQC.qc.checkIntersection._intersection_height

.. py:function:: _lines_cross(x_ctrl, y_ctrl, x_trav, y_trav)
   :canonical: AirGravQC.qc.checkIntersection._lines_cross

   .. autodoc2-docstring:: AirGravQC.qc.checkIntersection._lines_cross

.. py:function:: _calc_bearing(x, y)
   :canonical: AirGravQC.qc.checkIntersection._calc_bearing

   .. autodoc2-docstring:: AirGravQC.qc.checkIntersection._calc_bearing

.. py:function:: _mean_1std(x)
   :canonical: AirGravQC.qc.checkIntersection._mean_1std

   .. autodoc2-docstring:: AirGravQC.qc.checkIntersection._mean_1std

.. py:function:: controls_lessthan_1000(all_lines)
   :canonical: AirGravQC.qc.checkIntersection.controls_lessthan_1000

   .. autodoc2-docstring:: AirGravQC.qc.checkIntersection.controls_lessthan_1000

.. py:function:: controls_nineteen(all_lines)
   :canonical: AirGravQC.qc.checkIntersection.controls_nineteen

   .. autodoc2-docstring:: AirGravQC.qc.checkIntersection.controls_nineteen
