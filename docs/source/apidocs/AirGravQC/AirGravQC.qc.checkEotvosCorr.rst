:py:mod:`AirGravQC.qc.checkEotvosCorr`
======================================

.. py:module:: AirGravQC.qc.checkEotvosCorr

.. autodoc2-docstring:: AirGravQC.qc.checkEotvosCorr
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`checkEotvosCorr <AirGravQC.qc.checkEotvosCorr.checkEotvosCorr>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkEotvosCorr.checkEotvosCorr
          :summary:
   * - :py:obj:`_calc_speed <AirGravQC.qc.checkEotvosCorr._calc_speed>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkEotvosCorr._calc_speed
          :summary:
   * - :py:obj:`_eotvosCorrection <AirGravQC.qc.checkEotvosCorr._eotvosCorrection>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkEotvosCorr._eotvosCorrection
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`groupName <AirGravQC.qc.checkEotvosCorr.groupName>`
     - .. autodoc2-docstring:: AirGravQC.qc.checkEotvosCorr.groupName
          :summary:

API
~~~

.. py:data:: groupName
   :canonical: AirGravQC.qc.checkEotvosCorr.groupName
   :value: None

   .. autodoc2-docstring:: AirGravQC.qc.checkEotvosCorr.groupName

.. py:function:: checkEotvosCorr(whizzFile, eotCorr, latitude='', x='', y='', GRS80_height='', time='', east_vel='', north_vel='', plot_flag=False)
   :canonical: AirGravQC.qc.checkEotvosCorr.checkEotvosCorr

   .. autodoc2-docstring:: AirGravQC.qc.checkEotvosCorr.checkEotvosCorr

.. py:function:: _calc_speed(e, n, t)
   :canonical: AirGravQC.qc.checkEotvosCorr._calc_speed

   .. autodoc2-docstring:: AirGravQC.qc.checkEotvosCorr._calc_speed

.. py:function:: _eotvosCorrection(eSpeed, nSpeed, latitude, height=0)
   :canonical: AirGravQC.qc.checkEotvosCorr._eotvosCorrection

   .. autodoc2-docstring:: AirGravQC.qc.checkEotvosCorr._eotvosCorrection
