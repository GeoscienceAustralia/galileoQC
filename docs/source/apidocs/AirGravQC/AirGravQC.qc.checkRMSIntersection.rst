:py:mod:`AirGravQC.qc.checkRMSIntersection`
===========================================

.. py:module:: AirGravQC.qc.checkRMSIntersection

.. autodoc2-docstring:: AirGravQC.qc.checkRMSIntersection
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`checkRMSIntersection <AirGravQC.qc.checkRMSIntersection.checkRMSIntersection>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkRMSIntersection.checkRMSIntersection
          :summary:
   * - :py:obj:`_ccw <AirGravQC.qc.checkRMSIntersection._ccw>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkRMSIntersection._ccw
          :summary:
   * - :py:obj:`_intersect <AirGravQC.qc.checkRMSIntersection._intersect>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkRMSIntersection._intersect
          :summary:
   * - :py:obj:`_intersection_height <AirGravQC.qc.checkRMSIntersection._intersection_height>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkRMSIntersection._intersection_height
          :summary:
   * - :py:obj:`_lines_cross <AirGravQC.qc.checkRMSIntersection._lines_cross>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkRMSIntersection._lines_cross
          :summary:
   * - :py:obj:`_calc_bearing <AirGravQC.qc.checkRMSIntersection._calc_bearing>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkRMSIntersection._calc_bearing
          :summary:
   * - :py:obj:`_mean_1std <AirGravQC.qc.checkRMSIntersection._mean_1std>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkRMSIntersection._mean_1std
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`groupName <AirGravQC.qc.checkRMSIntersection.groupName>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkRMSIntersection.groupName
          :summary:

API
~~~

.. py:data:: groupName
   :canonical: AirGravQC.qc.checkRMSIntersection.groupName
   :value: None

   .. autodoc2-docstring:: AirGravQC.qc.checkRMSIntersection.groupName

.. py:function:: checkRMSIntersection(whizzFile, controls=[], travs=[], xChannel='', yChannel='', zChannel='', max_allowed_deltaZ=10.0)
   :canonical: AirGravQC.qc.checkRMSIntersection.checkRMSIntersection

   .. autodoc2-docstring:: AirGravQC.qc.checkRMSIntersection.checkRMSIntersection

.. py:function:: _ccw(x1, y1, x2, y2, x3, y3)
   :canonical: AirGravQC.qc.checkRMSIntersection._ccw

   .. autodoc2-docstring:: AirGravQC.qc.checkRMSIntersection._ccw

.. py:function:: _intersect(cx1, cy1, cx2, cy2, tx1, ty1, tx2, ty2)
   :canonical: AirGravQC.qc.checkRMSIntersection._intersect

   .. autodoc2-docstring:: AirGravQC.qc.checkRMSIntersection._intersect

.. py:function:: _intersection_height(x_trav, y_trav, z_trav, x_ctrl, y_ctrl, z_ctrl, bearingc)
   :canonical: AirGravQC.qc.checkRMSIntersection._intersection_height

   .. autodoc2-docstring:: AirGravQC.qc.checkRMSIntersection._intersection_height

.. py:function:: _lines_cross(x_ctrl, y_ctrl, x_trav, y_trav)
   :canonical: AirGravQC.qc.checkRMSIntersection._lines_cross

   .. autodoc2-docstring:: AirGravQC.qc.checkRMSIntersection._lines_cross

.. py:function:: _calc_bearing(x, y)
   :canonical: AirGravQC.qc.checkRMSIntersection._calc_bearing

   .. autodoc2-docstring:: AirGravQC.qc.checkRMSIntersection._calc_bearing

.. py:function:: _mean_1std(x)
   :canonical: AirGravQC.qc.checkRMSIntersection._mean_1std

   .. autodoc2-docstring:: AirGravQC.qc.checkRMSIntersection._mean_1std
