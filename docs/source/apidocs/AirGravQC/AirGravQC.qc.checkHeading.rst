:py:mod:`AirGravQC.qc.checkHeading`
===================================

.. py:module:: AirGravQC.qc.checkHeading

.. autodoc2-docstring:: AirGravQC.qc.checkHeading
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`checkHeading <AirGravQC.qc.checkHeading.checkHeading>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkHeading.checkHeading
          :summary:
   * - :py:obj:`angle_in_range <AirGravQC.qc.checkHeading.angle_in_range>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkHeading.angle_in_range
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`groupName <AirGravQC.qc.checkHeading.groupName>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkHeading.groupName
          :summary:

API
~~~

.. py:data:: groupName
   :canonical: AirGravQC.qc.checkHeading.groupName
   :value: None

   .. autodoc2-docstring:: AirGravQC.qc.checkHeading.groupName

.. py:function:: checkHeading(whizzFile, nominalHeadings, lines=[], headingchan='', x='', y='', tolerance=10.0, known='', plot_flag=False)
   :canonical: AirGravQC.qc.checkHeading.checkHeading

   .. autodoc2-docstring:: AirGravQC.qc.checkHeading.checkHeading

.. py:function:: angle_in_range(alpha, nominal, tolerance)
   :canonical: AirGravQC.qc.checkHeading.angle_in_range

   .. autodoc2-docstring:: AirGravQC.qc.checkHeading.angle_in_range
