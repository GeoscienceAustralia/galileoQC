:py:mod:`AirGravQC.qc.checkClearance`
=====================================

.. py:module:: AirGravQC.qc.checkClearance

.. autodoc2-docstring:: AirGravQC.qc.checkClearance
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`checkSafeClearance <AirGravQC.qc.checkClearance.checkSafeClearance>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkClearance.checkSafeClearance
          :summary:
   * - :py:obj:`checkClearance <AirGravQC.qc.checkClearance.checkClearance>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkClearance.checkClearance
          :summary:
   * - :py:obj:`checkDrape <AirGravQC.qc.checkClearance.checkDrape>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkClearance.checkDrape
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`groupName <AirGravQC.qc.checkClearance.groupName>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkClearance.groupName
          :summary:

API
~~~

.. py:data:: groupName
   :canonical: AirGravQC.qc.checkClearance.groupName
   :value: None

   .. autodoc2-docstring:: AirGravQC.qc.checkClearance.groupName

.. py:function:: checkSafeClearance(whizzFile, minimumAllowedClearance, clearance_chan='', altitude_chan='', terrain_chan='', xChannel='', yChannel='', plot_flag=False)
   :canonical: AirGravQC.qc.checkClearance.checkSafeClearance

   .. autodoc2-docstring:: AirGravQC.qc.checkClearance.checkSafeClearance

.. py:function:: checkClearance(whizzFile, nominalClearance, clearance_chan='', altitude_chan='', terrain_chan='', allowance=20.0, maxDistance=1000.0, xChannel='', yChannel='')
   :canonical: AirGravQC.qc.checkClearance.checkClearance

   .. autodoc2-docstring:: AirGravQC.qc.checkClearance.checkClearance

.. py:function:: checkDrape(whizzFile, altitude, drapeChannel, warningClearance=20.0, xChannel='', yChannel='', plot_flag=True)
   :canonical: AirGravQC.qc.checkClearance.checkDrape

   .. autodoc2-docstring:: AirGravQC.qc.checkClearance.checkDrape
