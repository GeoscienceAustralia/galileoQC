:py:mod:`AirGravQC.qc.checkXYPlan`
==================================

.. py:module:: AirGravQC.qc.checkXYPlan

.. autodoc2-docstring:: AirGravQC.qc.checkXYPlan
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`checkXYPlan <AirGravQC.qc.checkXYPlan.checkXYPlan>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkXYPlan.checkXYPlan
          :summary:
   * - :py:obj:`_plot_exceeding_line <AirGravQC.qc.checkXYPlan._plot_exceeding_line>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkXYPlan._plot_exceeding_line
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`groupName <AirGravQC.qc.checkXYPlan.groupName>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkXYPlan.groupName
          :summary:

API
~~~

.. py:data:: groupName
   :canonical: AirGravQC.qc.checkXYPlan.groupName
   :value: None

   .. autodoc2-docstring:: AirGravQC.qc.checkXYPlan.groupName

.. py:function:: checkXYPlan(planPath, measPath, lines=[], planX='', planY='', measX='', measY='', allowance=200.0, maxCounter=0, maxDistance=0, known='', plot_flag=False, verbose=False)
   :canonical: AirGravQC.qc.checkXYPlan.checkXYPlan

   .. autodoc2-docstring:: AirGravQC.qc.checkXYPlan.checkXYPlan

.. py:function:: _plot_exceeding_line(x, y, xP, yP, xM, yM, measX, measY, allowance, line, planLine, dirn)
   :canonical: AirGravQC.qc.checkXYPlan._plot_exceeding_line

   .. autodoc2-docstring:: AirGravQC.qc.checkXYPlan._plot_exceeding_line
