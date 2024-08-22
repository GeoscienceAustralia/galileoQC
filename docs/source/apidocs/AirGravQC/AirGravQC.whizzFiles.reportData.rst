:py:mod:`AirGravQC.whizzFiles.reportData`
=========================================

.. py:module:: AirGravQC.whizzFiles.reportData

.. autodoc2-docstring:: AirGravQC.whizzFiles.reportData
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`reportWhizz <AirGravQC.whizzFiles.reportData.reportWhizz>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.reportData.reportWhizz
          :summary:
   * - :py:obj:`reportFlights <AirGravQC.whizzFiles.reportData.reportFlights>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.reportData.reportFlights
          :summary:
   * - :py:obj:`reportSampling <AirGravQC.whizzFiles.reportData.reportSampling>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.reportData.reportSampling
          :summary:
   * - :py:obj:`_distanceFlown <AirGravQC.whizzFiles.reportData._distanceFlown>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.reportData._distanceFlown
          :summary:
   * - :py:obj:`_lineLength <AirGravQC.whizzFiles.reportData._lineLength>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.reportData._lineLength
          :summary:
   * - :py:obj:`whizzAttrExists <AirGravQC.whizzFiles.reportData.whizzAttrExists>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.reportData.whizzAttrExists
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`groupName <AirGravQC.whizzFiles.reportData.groupName>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.reportData.groupName
          :summary:
   * - :py:obj:`projectName <AirGravQC.whizzFiles.reportData.projectName>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.reportData.projectName
          :summary:

API
~~~

.. py:data:: groupName
   :canonical: AirGravQC.whizzFiles.reportData.groupName
   :value: None

   .. autodoc2-docstring:: AirGravQC.whizzFiles.reportData.groupName

.. py:data:: projectName
   :canonical: AirGravQC.whizzFiles.reportData.projectName
   :value: None

   .. autodoc2-docstring:: AirGravQC.whizzFiles.reportData.projectName

.. py:function:: reportWhizz(whizzFile, line='', channel='')
   :canonical: AirGravQC.whizzFiles.reportData.reportWhizz

   .. autodoc2-docstring:: AirGravQC.whizzFiles.reportData.reportWhizz

.. py:function:: reportFlights(whizzFile, flightChannel='', lines=[], detailed=False)
   :canonical: AirGravQC.whizzFiles.reportData.reportFlights

   .. autodoc2-docstring:: AirGravQC.whizzFiles.reportData.reportFlights

.. py:function:: reportSampling(whizzFile, timeChannel='', xChannel='', yChannel='')
   :canonical: AirGravQC.whizzFiles.reportData.reportSampling

   .. autodoc2-docstring:: AirGravQC.whizzFiles.reportData.reportSampling

.. py:function:: _distanceFlown(whizzFile, x='', y='', lines=[])
   :canonical: AirGravQC.whizzFiles.reportData._distanceFlown

   .. autodoc2-docstring:: AirGravQC.whizzFiles.reportData._distanceFlown

.. py:function:: _lineLength(x, y)
   :canonical: AirGravQC.whizzFiles.reportData._lineLength

   .. autodoc2-docstring:: AirGravQC.whizzFiles.reportData._lineLength

.. py:function:: whizzAttrExists(group, my_attr)
   :canonical: AirGravQC.whizzFiles.reportData.whizzAttrExists

   .. autodoc2-docstring:: AirGravQC.whizzFiles.reportData.whizzAttrExists
