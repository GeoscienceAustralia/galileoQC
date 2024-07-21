///////////////////////////////////////////////////
// Atlas Control/Workspace File
// # or // for commment
///////////////////////////////////////////////////

CONTROL_BEGIN

  PROGRAM  = LoadASCtoGDB
  VERSION  = 1.0.0


InputFile           = \\Mac\Home\Documents\GitHub\AirGravQC\examples\AirFTG\VintonDome.xyz

ImportTemplateFile  = 

OutputChannel       = X           ,  Yes ,  double[  1] ,       *
OutputChannel       = Y           ,  Yes ,  double[  1] ,       *
OutputChannel       = Lon         ,  Yes ,  double[  1] ,       *
OutputChannel       = Lat         ,  Yes ,  double[  1] ,       *
OutputChannel       = Altitude    ,  Yes ,   float[  1] ,       *
OutputChannel       = Drape       ,  Yes ,   float[  1] ,       *
OutputChannel       = Terrain     ,  Yes ,   float[  1] ,       *
OutputChannel       = Cross1_raw  ,  Yes ,   float[  1] ,       *
OutputChannel       = Cross2_raw  ,  Yes ,   float[  1] ,       *
OutputChannel       = Cross3_raw  ,  Yes ,   float[  1] ,       *
OutputChannel       = Inline1_raw ,  Yes ,   float[  1] ,       *
OutputChannel       = Inline2_raw ,  Yes ,   float[  1] ,       *
OutputChannel       = Inline3_raw ,  Yes ,   float[  1] ,       *
OutputChannel       = TC_Txx_100  ,  Yes ,   float[  1] ,       *
OutputChannel       = TC_Txz_100  ,  Yes ,   float[  1] ,       *
OutputChannel       = TC_Tyx_100  ,  Yes ,   float[  1] ,       *
OutputChannel       = TC_Tyy_100  ,  Yes ,   float[  1] ,       *
OutputChannel       = TC_Tyz_100  ,  Yes ,   float[  1] ,       *
OutputChannel       = TC_Tzz_100  ,  Yes ,   float[  1] ,       *
OutputChannel       = Txx_slv     ,  Yes ,   float[  1] ,       *
OutputChannel       = Txz_slv     ,  Yes ,   float[  1] ,       *
OutputChannel       = Tyx_slv     ,  Yes ,   float[  1] ,       *
OutputChannel       = Tyy_slv     ,  Yes ,   float[  1] ,       *
OutputChannel       = Tyz_slv     ,  Yes ,   float[  1] ,       *
OutputChannel       = Tzz_slv     ,  Yes ,   float[  1] ,       *
OutputChannel       = HHMMSS      ,  Yes ,  double[  1] ,       *
OutputChannel       = Time        ,  Yes ,  double[  1] ,       *
OutputChannel       = YYMMDD      ,  Yes ,  double[  1] ,       *

OutputGDB           = \\Mac\Home\Documents\GitHub\AirGravQC\examples\AirFTG\VintonDome.gdb

  GDBCompression = 1
  GDBMaxChannels = 256
  GDBMaxLines    = 4095
  GDBMaxBlobs    = 4851

TextDelimiters      = @20,	


InputFileType       = GsfXYZ
DummyValue          = *

LinesToSkip         = 10
ChannelHeaderLine   = 6

LineChannelColumn   = 0
FlightChannelColumn = 0
FidChannelColumn    = 0
DateChannelColumn   = 0

DataSyncOption      = NoSync
  SyncChannelName   = 
  SyncChannelColumn = 0
  SyncOffset        = 
  SyncTolerance     = 
  SyncAdjustFracFid = No

  XChannelColumn    = 0
  YChannelColumn    = 0

  XYSearchDistance   = 500
  XYSearchLineExt    = 100
  AccurateLineSearch = No


CONTROL_END

