# -*- coding: utf-8 -*-
"""
Created on Sun Jan 31 14:14:24 2016

@author: William
"""
import os
import sys
import numpy as np
import scipy as sp
from PIL import Image

try:
    from osgeo import ogr, osr, gdal
except:
    sys.exit('ERROR: cannot find GDAL/OGR modules')


#==============================================================================
# GDAL Error handler
#==============================================================================
# Enable GDAL/OGR exceptions
#gdal.UseExceptions()

# example GDAL error handler function
def gdal_error_handler(err_class, err_num, err_msg):
    errtype = {
            gdal.CE_None:'None',
            gdal.CE_Debug:'Debug',
            gdal.CE_Warning:'Warning',
            gdal.CE_Failure:'Failure',
            gdal.CE_Fatal:'Fatal'
    }
    err_msg = err_msg.replace('\n',' ')
    err_class = errtype.get(err_class, 'None')
    print 'Error Number: %s' % (err_num)
    print 'Error Type: %s' % (err_class)
    print 'Error Message: %s' % (err_msg)


        # install error handler
        #gdal.PushErrorHandler(gdal_error_handler)
        
        # Raise a dummy error
        #gdal.Error(1, 2, 'test error')
    
        #uninstall error handler
        #gdal.PopErrorHandler()

#==============================================================================
# test datatype
#==============================================================================

def TestDataType(data):
    if isinstance(data,gdal.Dataset):
        img = data        
    elif isinstance(data,str):
        img = gdal.Open(data,gdal.GA_Update)        
    else:
        print "error: wrong datatype"
    return img
    
#==============================================================================
#     info functions
#==============================================================================

def GetImageInfo(filename):
    img = gdal.Open(filename,gdal.GA_Update)
    
    return img
    
def GetImageSize(filename):
    img = TestDataType(filename)
    nvar = img.RasterCount
    ncol = img.RasterXSize
    nrow = img.RasterYSize
    Size = (nrow,ncol,nvar)
    return Size
    
#==============================================================================
#     Create file
#==============================================================================

def CreateFileAs(filename,img,**warg):
    img = TestDataType(img)
    driver = img.GetDriver()
    projection = img.GetProjection()
    geoTransform = img.GetGeoTransform()
    nrow, ncol, nvar = GetImageSize(img)
    
    if "Nband" in warg:
        nvar = warg["Nband"]
    
    #if file exist delete it first
    if os.path.isfile(filename):
        os.remove(filename)
     
    #create a new
    Outdata = driver.Create(filename,ncol,nrow,nvar,gdal.GDT_Float32)  
    Outdata.SetProjection(projection)
    Outdata.SetGeoTransform(geoTransform)

    if "data" in warg:
        data = warg["data"]
    else:
        data = np.zeros((nrow,ncol,nvar))    
        
    WriteToFile(Outdata,data)
        
    return Outdata
    
#==============================================================================
#     Write functions 
#==============================================================================
    
def WriteToFile(filename,data):
    OutImg = TestDataType(filename)
    nrow, ncol, nvar = GetImageSize(OutImg)
    
    # Define NoData value of new raster
    NoData_value = -9999

    if data.shape == GetImageSize(OutImg):
        if nvar == 1:
            if len(np.shape(data))==3:
                OutImg.GetRasterBand(1).WriteArray(data[:,:,0],0,0)
                #set no data value
                OutImg.GetRasterBand(1).SetNoDataValue(NoData_value)
    
            
            elif len(np.shape(data))==2:
                OutImg.GetRasterBand(1).WriteArray(data,0,0)
                 #set no data value
                OutImg.GetRasterBand(1).SetNoDataValue(NoData_value)
    
            else:
                print "size of data is wrong"    
        else:
            for i in range(1, nvar+1):
                OutImg.GetRasterBand(i).WriteArray(data[:,:,i-1],0,0)
                #set no data value
                OutImg.GetRasterBand(i).SetNoDataValue(NoData_value)
    
    else:
        print "Error: image size mismatch" 
    
    OutImg.FlushCache()
    
    return None
    
def WriteImgRow(filename,data,RowNr):
    img = TestDataType(filename)
    Size = data.shape
    if (Size[0] == 1) & (len(Size) == 2):
        img.GetRasterBand(1).WriteArray(data,0,RowNr)
        
    elif Size[1] == GetImageSize(img)[1]:
        for i in range(1,Size[0]+1):
            img.GetRasterBand(i).WriteArray(data[i-1,np.newaxis],0,RowNr)
    else:
        print "Error: image size mismatch" 
    
    img.FlushCache()
    
    return None
    
#==============================================================================
#     read functions
#==============================================================================
    
def ReadData(filename):
    img = TestDataType(filename)
    data1 = np.array(img.GetRasterBand(1).ReadAsArray())
    data = data1[:,:,np.newaxis]
    if img.RasterCount > 1:
        for i in range(2, img.RasterCount + 1):
            data_array = np.array(img.GetRasterBand(i).ReadAsArray())
            data = np.concatenate((data,data_array[:,:,np.newaxis]),axis = 2)

    return data
    
def ReadImgRow(filename,RowNr):
    img = TestDataType(filename)
    nrow, ncol, nvar = GetImageSize(img)
    
    if nvar == 1:
        dataarray = np.array(img.GetRasterBand(1).ReadAsArray(0,RowNr,ncol,1))
    else:
        dataarray = np.empty((0,ncol))
        for i in range(1, nvar+1):
            data = np.array(img.GetRasterBand(i).ReadAsArray(0,RowNr,ncol,1))
            dataarray = np.vstack((dataarray,data))
 
    return dataarray
    
def ReadImgRect(filename,SquareInfo):
    Sq = SquareInfo
    img = TestDataType(filename)
    nrow, ncol, nvar = GetImageSize(img)
    
    if nvar == 1:
        dataarray = np.array(img.GetRasterBand(1).ReadAsArray(Sq[0],Sq[1],Sq[2],Sq[3]))
    else:
        dataarray = np.empty((0,ncol))
        for i in range(1, nvar+1):
            data = np.array(img.GetRasterBand(i).ReadAsArray(Sq[0],Sq[1],Sq[2],Sq[3]))
            dataarray = np.dstack((dataarray,data))
 
    return dataarray


        
        
#==============================================================================
# Polygon to binary
#==============================================================================

def polygon2Binary(inputfn,asraster):
    raster = TestDataType(asraster)
    geoTransform = raster.GetGeoTransform()
    [nrow,ncol,nvar] = GetImageSize(raster)
    
    # Open the data source
    source_ds = ogr.Open(inputfn)
    source_layer = source_ds.GetLayer()
    
    # Create the destination data source
    target_ds = gdal.GetDriverByName('MEM').Create('', ncol, nrow, gdal.GDT_Byte)
    target_ds.SetGeoTransform(geoTransform)
    band = target_ds.GetRasterBand(1)
    
    # Define NoData value of new raster
    NoData_value = -9999
    band.SetNoDataValue(NoData_value)
    
    # Rasterize
    gdal.RasterizeLayer(target_ds, [1], source_layer,None, burn_values=[1])
    
    # Read as array
    binary = band.ReadAsArray()
    
    return binary
    


#==============================================================================
# Convert Shape to another projection
#==============================================================================

def ShapeProjConvert(inputfn,projname,outputfn):
    
    if projname == 'UTM32':
        EPSGcode = 25832
    elif projname == 'WGS84':
        EPSGcode = 4326
    elif projname == 'ETRS89':
        EPSGcode = 3044
    else:
        print 'wrong datatype'
        
    #read input file and layer
    inDriver = ogr.GetDriverByName("ESRI Shapefile")
    inDataSource = inDriver.Open(inputfn, 0)
    inLayer = inDataSource.GetLayer()
    
    # Create the output Layer
    outDriver = ogr.GetDriverByName("ESRI Shapefile")
    
    # Remove output shapefile if it already exists
    if os.path.exists(outputfn):
        outDriver.DeleteDataSource(outputfn)
        
    #get datatype 
    lyr = inLayer.GetFeature(0)
    geom = lyr.GetGeometryRef()
    geomType = geom.GetGeometryType()
            

    #calculate the coordinate transformation
    sourceSR = inLayer.GetSpatialRef()
    targetSR = osr.SpatialReference()
    targetSR.ImportFromEPSG(EPSGcode)
    coordTrans = osr.CoordinateTransformation(sourceSR,targetSR)
    
    # Create the output shapefile
    outDataSource = outDriver.CreateDataSource(outputfn)
    outLayer = outDataSource.CreateLayer("myLayer",targetSR, geom_type=geomType)
    

    # Get the output Layer's Feature Definition
    outLayerDefn = outLayer.GetLayerDefn()
    

    
    # Add features to the ouput Layer
    for i in range(0, inLayer.GetFeatureCount()):
        # Get the input Feature
        inFeature = inLayer.GetFeature(i)
        # Create output Feature
        
        #get geometry from input file and do coordinate transformation
        geom = inFeature.GetGeometryRef()
        geom.Transform(coordTrans)
        
        #create feature
        outFeature = ogr.Feature(outLayerDefn)
        #add geometry to feature
        outFeature.SetGeometry(geom)
        # Add new feature to output Layer
        outLayer.CreateFeature(outFeature)
    
    # Close DataSources
    inDataSource.Destroy()
    outDataSource.Destroy()

    return 0
    
  
#==============================================================================
# Geometry projection conversion
#==============================================================================

def GeomProjConvert(geom,target,source):
    if target == 'UTM32':
        EPSGcodeTarget = 32632
    elif target == 'WGS84':
        EPSGcodeTarget = 4326
    else:
        print 'wrong datatype'
    
    if source == 'UTM32':
        EPSGcodeSource = 32632
    elif source == 'WGS84':
        EPSGcodeSource = 4326
    else:
        print 'wrong datatype'
        
    targetSR = osr.SpatialReference()
    sourceSR = osr.SpatialReference()
    targetSR.ImportFromEPSG(EPSGcodeTarget)
    sourceSR.ImportFromEPSG(EPSGcodeSource)
    
    coordTrans = osr.CoordinateTransformation(targetSR,sourceSR)
    geom.Transform(coordTrans)
    
    return geom

#==============================================================================
# Save geometry to shape file
#==============================================================================
    
def SaveGeomToShape(shapefn,inGeom):
    # Create the output Layer driver
    outDriver = ogr.GetDriverByName("ESRI Shapefile")
    
    # Remove output shapefile if it already exists
    if os.path.exists(shapefn):
        outDriver.DeleteDataSource(shapefn)
    
    if isinstance(inGeom,ogr.Geometry):
        geotest = inGeom
    elif isinstance(inGeom,list):
        Geotest = inGeom[0]
        
    geomtype = geotest.GetGeometryType()
 

    # Create the output shapefile
    outDataSource = outDriver.CreateDataSource(shapefn)
    outLayer = outDataSource.CreateLayer("mylayer", geom_type=geomtype)
    outLayerDefn =  outLayer.GetLayerDefn()
    
    if isinstance(inGeom,ogr.Geometry):
        outFeature = ogr.Feature(outLayerDefn)        
        outFeature.SetGeometry(inGeom)
        outLayer.CreateFeature(outFeature)
    elif isinstance(inGeom,list):
        # Add features to the ouput Layer
        for i in range(0, len(inGeom)):
            # Create output Feature
            outFeature = ogr.Feature(outLayerDefn)        
            outFeature.SetGeometry(inGeom[i])
            # Add new feature to output Layer
            outLayer.CreateFeature(outFeature)
            
    # Close DataSources
    outDataSource.Destroy()
    return 0
    
#==============================================================================
# CreateBuffer
#==============================================================================

def createBuffer(inputfn, outputBufferfn, bufferDist):
    inputds = ogr.Open(inputfn)
    inputlyr = inputds.GetLayer()

    shpdriver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(outputBufferfn):
        shpdriver.DeleteDataSource(outputBufferfn)
    outputBufferds = shpdriver.CreateDataSource(outputBufferfn)
    bufferlyr = outputBufferds.CreateLayer(outputBufferfn, geom_type=ogr.wkbPolygon)
    featureDefn = bufferlyr.GetLayerDefn()

    for feature in inputlyr:
        ingeom = feature.GetGeometryRef()
        geomBuffer = ingeom.Buffer(bufferDist)

        outFeature = ogr.Feature(featureDefn)
        outFeature.SetGeometry(geomBuffer)
        bufferlyr.CreateFeature(outFeature)
        
    return 0
   
#==============================================================================
# Maximum filter of raster data
#==============================================================================
def MaxFilter(Rasterin,Rasterout,WinSize):
    #WinSize must be odd number
    Rasterin = TestDataType(Rasterin)
    nrow,ncol,nvar = GetImageSize(Rasterin)
    Mrow = np.round((WinSize-1)/2) #the middle row number
    
    #create circle buffer at the size of window
    ax = range(WinSize)-Mrow
    
    ax1 = np.tile(np.mat(ax).T,(1,WinSize))
    ax2 = np.tile(ax,(WinSize,1))
    
    circ = np.sqrt(np.array(ax1)**2+ax2**2)
    
    blob = Mrow>=circ
    
        
    datarow = np.zeros((WinSize,ncol+2*Mrow))
    OutRow = np.zeros((1,ncol))
    
    for irow in range(nrow):
        print irow       
        
        for d in range(WinSize):
            ReadRow = irow + d - Mrow
            if 0 <= ReadRow | ReadRow < nrow: #
                datarow[d,Mrow:-Mrow] = ReadImgRow(Rasterin,ReadRow) 
            else:
                datarow[d,:] = 0
            
        for icol in range(ncol):
            win = datarow[0:WinSize,icol:icol+WinSize]
            OutRow[:,icol] = np.max(win*blob)
            
        

        WriteImgRow(Rasterout,OutRow,irow)
    
    
    return 0
    
    
    
def MaxGaussFilter(Rasterin,Rasterout,Sigma):

    WinSize = 1 + Sigma*4
    
    #WinSize must be odd number
    Rasterin = TestDataType(Rasterin)
    nrow,ncol,nvar = GetImageSize(Rasterin)
    Mrow = np.round((WinSize-1)/2) #the middle row number
    
    x = range(0,WinSize)
    Gaus = sp.stats.norm.pdf(x, loc=Mrow, scale=Sigma)
    
    G = np.mat(Gaus*Gaus)
    W = G.T*G
    W = np.sqrt(np.array(W))
    W /= W[Mrow,Mrow]    
    
    datarow = np.zeros((WinSize,ncol+2*Mrow))
    OutRow = np.zeros((1,ncol))
    
    for irow in range(nrow):
        print irow       
        
        for d in range(WinSize):
            ReadRow = irow + d - Mrow
            if 0 <= ReadRow | ReadRow < nrow: #
                datarow[d,Mrow:-Mrow] = ReadImgRow(Rasterin,ReadRow) 
            else:
                datarow[d,:] = 0
            
        for icol in range(ncol):
            win = datarow[0:WinSize,icol:icol+WinSize]
            OutRow[:,icol] = np.max(win*W)

        WriteImgRow(Rasterout,OutRow,irow)
    
    
    return 0
    
#==============================================================================
# Get no data value
#==============================================================================

def GetNoDataValue(filename):
    raster = TestDataType(filename)    
    NoDataValue = raster.GetRasterBand(1).GetNoDataValue()
    
    return NoDataValue

#==============================================================================
# raster to polygon
#==============================================================================
   
def RasterToPolygon(Raster,outShapefile):
    gdal.UseExceptions()
    sourceRaster = TestDataType(Raster)
    band = sourceRaster.GetRasterBand(1)

    driver = ogr.GetDriverByName("ESRI Shapefile")
    if os.path.exists(outShapefile):
        driver.DeleteDataSource(outShapefile)
    
    #UTM32
    EPSGcode = 25832
    targetSR = osr.SpatialReference()
    targetSR.ImportFromEPSG(EPSGcode)
    
    outDatasource = driver.CreateDataSource(outShapefile)
    outLayer = outDatasource.CreateLayer("polygonized",targetSR)#, srs=None)
    
    newField = ogr.FieldDefn('MYFLD', ogr.OFTInteger)
    outLayer.CreateField(newField)  
    
    gdal.Polygonize( band, None, outLayer, 0, [], callback=None )
    outDatasource.Destroy()
    sourceRaster = None
    
    return outShapefile
    
#==============================================================================
# Minimum distance to object
#==============================================================================

def MinDistanceToSHP(obj,shape,bufferDist):
    
    driver = ogr.GetDriverByName('ESRI Shapefile')
    
    #open skel shpfile
    Skeldata = driver.Open(shape,0)       
    SkelLayer = Skeldata.GetLayer()
    
    #open change shp file

    if isinstance(obj,ogr.Geometry):
        ChngData = obj   
    elif isinstance(obj,str):
        ChngData = driver.Open(obj, 0) # 0 means read-only. 1 means writeable.        
    else:
        print "error: wrong datatype"    
    
    ChngLayer = ChngData.GetLayer()
    ChngFeatCount = ChngLayer.GetFeatureCount()
    
    SkelNearfn = 'shape/skel/skelNear.shp'
    if os.path.exists(SkelNearfn):
        driver.DeleteDataSource(SkelNearfn)
        
    SkelNearData = driver.CreateDataSource(SkelNearfn)
    SkelNearLayer = SkelNearData.CreateLayer(SkelNearfn,geom_type = ogr.wkbLineString)
    
    dist = np.zeros(ChngFeatCount)+9999999
    index = np.zeros(ChngFeatCount,dtype=int)
    
    #loop for each Change feature
    for i in range(ChngFeatCount):    
        ChngFeat = ChngLayer.GetFeature(i)
        ChngGeom = ChngFeat.geometry()
        
        if bufferDist is not None:
            SkelLayer.SetSpatialFilter(ChngGeom.Buffer(bufferDist))
        
        SkelFeatCount = SkelLayer.GetFeatureCount()
        
        #find distance to lines within buffer dist
        for j in range(SkelFeatCount):
            SkelFeat = SkelLayer.GetFeature(j)
            SkelGeom = SkelFeat.geometry()        
            distj = ogr.Geometry.Distance(SkelGeom,ChngGeom)
            #save the smallest one        
            if dist[i] > distj:        
                dist[i] = distj
                index[i] = j
    
        SkelFeat = SkelLayer.GetFeature(index[i])
        
        SkelNearFeat = ogr.Feature(SkelNearLayer.GetLayerDefn())
        SkelNearFeat.SetGeometry(SkelFeat.GetGeometryRef().Clone())
        SkelNearLayer.CreateFeature(SkelNearFeat)
        SkelNearLayer.SyncToDisk()
                

        
    return dist
    
    
#==============================================================================
#     Max Area of polygon Removes objects with area smaller than AreaMin
#==============================================================================
def RemoveMinAreaPolygon(inshapefn,AreaMin):
    
    driver = ogr.GetDriverByName('ESRI Shapefile')
    
    #open skel shpfile
    inputdata = driver.Open(inshapefn,1) # 0 means read-only. 1 means writeable.      
    inputlayer = inputdata.GetLayer()
    
    ChngFeatCount = inputlayer.GetFeatureCount()
    
    #loop for each Change feature
    for i in range(ChngFeatCount):    
        inputfeat = inputlayer.GetFeature(i)
        ChngGeom = inputfeat.geometry()
        Area = ogr.Geometry.Area(ChngGeom)
        if (Area <= AreaMin) | (Area >= 1000):
            featId = inputfeat.GetFID()
            inputlayer.DeleteFeature(featId)
            
    return 0

#==============================================================================
# other functions
#==============================================================================
def Scale2uint8(I):
    I_out = (((I - I.min()) / (I.max() - I.min())) * 255.9)
    return I_out.astype(np.uint8)


#==============================================================================
# Select a list of object from shape using index
#==============================================================================
def SelectVectorObjects(index,fromShpfn,outputfn):
    
    driver = ogr.GetDriverByName('ESRI Shapefile')

    geometryType = ogr.wkbPolygon

    
    if os.path.exists(outputfn):
        driver.DeleteDataSource(outputfn)
    
    out_ds = driver.CreateDataSource(outputfn)
    out_layer = out_ds.CreateLayer(outputfn, geom_type=geometryType)
    
    #open inputfile
    data = driver.Open(fromShpfn,0)
    Layer = data.GetLayer()
    
    for i in range(Layer.GetFeatureCount()):
        if (index == i).any:
            feat = Layer.GetFeature(i)
            out_feat = ogr.Feature(out_layer.GetLayerDefn())
            out_feat.SetGeometry(feat.GetGeometryRef().Clone())
            out_layer.CreateFeature(out_feat)
            out_layer.SyncToDisk()
    return 0

       