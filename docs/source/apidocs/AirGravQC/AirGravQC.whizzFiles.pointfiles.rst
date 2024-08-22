:py:mod:`AirGravQC.whizzFiles.pointfiles`
=========================================

.. py:module:: AirGravQC.whizzFiles.pointfiles

.. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles
   :allowtitles:

Module Contents
---------------

Functions
~~~~~~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`_translate_date <AirGravQC.whizzFiles.pointfiles._translate_date>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles._translate_date
          :summary:
   * - :py:obj:`updateChannelAttributes <AirGravQC.whizzFiles.pointfiles.updateChannelAttributes>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.updateChannelAttributes
          :summary:
   * - :py:obj:`interpolateGridOntoLine <AirGravQC.whizzFiles.pointfiles.interpolateGridOntoLine>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.interpolateGridOntoLine
          :summary:
   * - :py:obj:`interpolateLine <AirGravQC.whizzFiles.pointfiles.interpolateLine>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.interpolateLine
          :summary:
   * - :py:obj:`_trim_monotonic <AirGravQC.whizzFiles.pointfiles._trim_monotonic>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles._trim_monotonic
          :summary:
   * - :py:obj:`_invDist <AirGravQC.whizzFiles.pointfiles._invDist>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles._invDist
          :summary:
   * - :py:obj:`_weightedAverage <AirGravQC.whizzFiles.pointfiles._weightedAverage>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles._weightedAverage
          :summary:
   * - :py:obj:`updateProject <AirGravQC.whizzFiles.pointfiles.updateProject>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.updateProject
          :summary:
   * - :py:obj:`updateCoordFrame <AirGravQC.whizzFiles.pointfiles.updateCoordFrame>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.updateCoordFrame
          :summary:
   * - :py:obj:`renameChannels <AirGravQC.whizzFiles.pointfiles.renameChannels>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.renameChannels
          :summary:
   * - :py:obj:`addWhizzToWhizz <AirGravQC.whizzFiles.pointfiles.addWhizzToWhizz>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.addWhizzToWhizz
          :summary:
   * - :py:obj:`addLineToWhizz <AirGravQC.whizzFiles.pointfiles.addLineToWhizz>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.addLineToWhizz
          :summary:
   * - :py:obj:`doyToISO8601 <AirGravQC.whizzFiles.pointfiles.doyToISO8601>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.doyToISO8601
          :summary:

Data
~~~~

.. list-table::
   :class: autosummary longtable
   :align: left

   * - :py:obj:`groupName <AirGravQC.whizzFiles.pointfiles.groupName>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.groupName
          :summary:
   * - :py:obj:`projectName <AirGravQC.whizzFiles.pointfiles.projectName>`
     - .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.projectName
          :summary:

API
~~~

.. py:data:: groupName
   :canonical: AirGravQC.whizzFiles.pointfiles.groupName
   :value: None

   .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.groupName

.. py:data:: projectName
   :canonical: AirGravQC.whizzFiles.pointfiles.projectName
   :value: None

   .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.projectName

.. py:function:: _translate_date(decimal_year)
   :canonical: AirGravQC.whizzFiles.pointfiles._translate_date

   .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles._translate_date

.. py:function:: updateChannelAttributes(whizzFile, channel, name='', units='', alias='', description='', chan_precision=-1)
   :canonical: AirGravQC.whizzFiles.pointfiles.updateChannelAttributes

   .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.updateChannelAttributes

.. py:function:: interpolateGridOntoLine(gridPath, hdfPath, lines=[])
   :canonical: AirGravQC.whizzFiles.pointfiles.interpolateGridOntoLine

   .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.interpolateGridOntoLine

.. py:function:: interpolateLine(timeIn, dataIn, timeOut, spare=[], plot_flag=False)
   :canonical: AirGravQC.whizzFiles.pointfiles.interpolateLine

   .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.interpolateLine

.. py:function:: _trim_monotonic(data_in, sync=[])
   :canonical: AirGravQC.whizzFiles.pointfiles._trim_monotonic

   .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles._trim_monotonic

.. py:function:: _invDist(x, y)
   :canonical: AirGravQC.whizzFiles.pointfiles._invDist

   .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles._invDist

.. py:function:: _weightedAverage(x, y, z, xo, yo)
   :canonical: AirGravQC.whizzFiles.pointfiles._weightedAverage

   .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles._weightedAverage

.. py:function:: updateProject(whizzFile, projectName='', blockID='', acquirer='', acquirerProjectID='', reportName='')
   :canonical: AirGravQC.whizzFiles.pointfiles.updateProject

   .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.updateProject

.. py:function:: updateCoordFrame(whizzFile, lat='', lon='', geoDatum='', alt='', htDatum='', x='', y='', projection='', utmz='', time='', timeDatum='', fid='')
   :canonical: AirGravQC.whizzFiles.pointfiles.updateCoordFrame

   .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.updateCoordFrame

.. py:function:: renameChannels(nchannels, chanNames)
   :canonical: AirGravQC.whizzFiles.pointfiles.renameChannels

   .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.renameChannels

.. py:function:: addWhizzToWhizz(inputWhizzFile, outputWhizzFile, lines=[])
   :canonical: AirGravQC.whizzFiles.pointfiles.addWhizzToWhizz

   .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.addWhizzToWhizz

.. py:function:: addLineToWhizz(hdf5FileId, databaseName, lineNo, lineType, plannedX=[], plannedY=[], plannedZ=[], distanceUnits='')
   :canonical: AirGravQC.whizzFiles.pointfiles.addLineToWhizz

   .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.addLineToWhizz

.. py:function:: doyToISO8601(doy)
   :canonical: AirGravQC.whizzFiles.pointfiles.doyToISO8601

   .. autodoc2-docstring:: AirGravQC.whizzFiles.pointfiles.doyToISO8601
