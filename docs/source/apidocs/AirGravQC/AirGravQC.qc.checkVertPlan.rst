:py:mod:`AirGravQC.qc.checkVertPlan`
====================================

.. py:module:: AirGravQC.qc.checkVertPlan

.. autodoc2-docstring:: AirGravQC.qc.checkVertPlan
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`checkVertPlan <AirGravQC.qc.checkVertPlan.checkVertPlan>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkVertPlan.checkVertPlan
          :summary:
   * - :py:obj:`_plot_vert_exceedance <AirGravQC.qc.checkVertPlan._plot_vert_exceedance>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkVertPlan._plot_vert_exceedance
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`groupName <AirGravQC.qc.checkVertPlan.groupName>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkVertPlan.groupName
          :summary:

API
~~~

.. py:data:: groupName
   :canonical: AirGravQC.qc.checkVertPlan.groupName
   :value: None

   .. autodoc2-docstring:: AirGravQC.qc.checkVertPlan.groupName

.. py:function:: checkVertPlan(planPath, measPath, lines=[], planX='', planY='', planZ='', measX='', measY='', measZ='', allowance=30.0, maxCounter=13, maxDistance=0.0, known='', plot_flag=False)
   :canonical: AirGravQC.qc.checkVertPlan.checkVertPlan

   .. autodoc2-docstring:: AirGravQC.qc.checkVertPlan.checkVertPlan

.. py:function:: _plot_vert_exceedance(xm, z_dev, xP, zP, xM, zM, measX, measZ, allowance, line, planLine, dirn)
   :canonical: AirGravQC.qc.checkVertPlan._plot_vert_exceedance

   .. autodoc2-docstring:: AirGravQC.qc.checkVertPlan._plot_vert_exceedance
