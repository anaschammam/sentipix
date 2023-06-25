import os
from qgis.core import QgsProject, QgsRasterLayer
from scipy.ndimage import zoom
from osgeo import gdal
import numpy as np
import math as m
def calculate_ndvi(data_source_path, red_band_path, nir_band_path):
    # Open the red band raster
    red_dataset = gdal.Open(red_band_path)
    red_band = red_dataset.GetRasterBand(1)
    red_width = red_dataset.RasterXSize
    red_height = red_dataset.RasterYSize

    # Open the NIR band raster
    nir_dataset = gdal.Open(nir_band_path)
    nir_band = nir_dataset.GetRasterBand(1)
    nir_width = nir_dataset.RasterXSize
    nir_height = nir_dataset.RasterYSize

    # Determine the output dimensions
    output_width = min(red_width, nir_width)
    output_height = min(red_height, nir_height)

    if red_width != nir_width or red_height != nir_height:
        # Calculate the scaling factors in x and y dimensions
        # Resample all input arrays to the same size using the nearest-neighbor method
        if red_width > nir_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            red_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            red_dataset_resampled.SetGeoTransform(red_dataset.GetGeoTransform())
            red_dataset_resampled.SetProjection(red_dataset.GetProjection())
            red_dataset_resampled.GetRasterBand(1).Fill(0)
            red_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(red_dataset, red_dataset_resampled, None, None, resample_method)
            red_band = red_dataset_resampled.GetRasterBand(1)
        else:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            nir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1, gdal.GDT_Float32)
            nir_dataset_resampled.SetGeoTransform(nir_dataset.GetGeoTransform())
            nir_dataset_resampled.SetProjection(nir_dataset.GetProjection())
            nir_dataset_resampled.GetRasterBand(1).Fill(0)
            nir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(nir_dataset, nir_dataset_resampled, None, None, resample_method)
            nir_band = nir_dataset_resampled.GetRasterBand(1)


    # Read the band arrays as NumPy arrays
    red_array = red_band.ReadAsArray().astype(np.float32)
    nir_array = nir_band.ReadAsArray().astype(np.float32)

    # Calculate NDVI
    ndvi = (nir_array - red_array) / (nir_array + red_array)

    # Create the output file
    output_path = data_source_path
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_path, output_width, output_height, 1, gdal.GDT_Float32)
    output_dataset.GetRasterBand(1).WriteArray(ndvi)
    output_dataset.SetGeoTransform(red_dataset.GetGeoTransform())
    output_dataset.SetProjection(red_dataset.GetProjection())
    output_dataset.FlushCache()
    output_dataset = None
    red_dataset = None
    nir_dataset = None

    # Create the output file path
    output_folder = QgsProject.instance().homePath()  # Get the QGIS project's home folder
    output_file_path = os.path.join(output_folder, output_path)

    # Add the output raster to QGIS
    output_layer = QgsRasterLayer(output_file_path, "NDVI")
    QgsProject.instance().addMapLayer(output_layer)
    print("NDVI calculated and added to QGIS.")

def calculate_bai(data_source_path, red_band_path, nir_band_path):
    # Open the red band raster
    red_dataset = gdal.Open(red_band_path)
    red_band = red_dataset.GetRasterBand(1)
    red_width = red_dataset.RasterXSize
    red_height = red_dataset.RasterYSize

    # Open the NIR band raster
    nir_dataset = gdal.Open(nir_band_path)
    nir_band = nir_dataset.GetRasterBand(1)
    nir_width = nir_dataset.RasterXSize
    nir_height = nir_dataset.RasterYSize

    # Determine the output dimensions
    output_width = min(red_width, nir_width)
    output_height = min(red_height, nir_height)

    if red_width != nir_width or red_height != nir_height:
        # Calculate the scaling factors in x and y dimensions
        # Resample all input arrays to the same size using the nearest-neighbor method
        if red_width > nir_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            red_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            red_dataset_resampled.SetGeoTransform(red_dataset.GetGeoTransform())
            red_dataset_resampled.SetProjection(red_dataset.GetProjection())
            red_dataset_resampled.GetRasterBand(1).Fill(0)
            red_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(red_dataset, red_dataset_resampled, None, None, resample_method)
            red_band = red_dataset_resampled.GetRasterBand(1)
        else:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            nir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1, gdal.GDT_Float32)
            nir_dataset_resampled.SetGeoTransform(nir_dataset.GetGeoTransform())
            nir_dataset_resampled.SetProjection(nir_dataset.GetProjection())
            nir_dataset_resampled.GetRasterBand(1).Fill(0)
            nir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(nir_dataset, nir_dataset_resampled, None, None, resample_method)
            nir_band = nir_dataset_resampled.GetRasterBand(1)


    # Read the band arrays as NumPy arrays
    red_array = red_band.ReadAsArray().astype(np.float32)
    nir_array = nir_band.ReadAsArray().astype(np.float32)

    # Calculate BAI = 1/((0.1-RED)²+(0.06-NIR)²)
    bai =1 / ( (0.1-red_array)**2 + (0.06-nir_array)**2 )

    # Create the output file
    output_path = data_source_path
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_path, output_width, output_height, 1, gdal.GDT_Float32)
    output_dataset.GetRasterBand(1).WriteArray(bai)
    output_dataset.SetGeoTransform(red_dataset.GetGeoTransform())
    output_dataset.SetProjection(red_dataset.GetProjection())
    output_dataset.FlushCache()
    output_dataset = None
    red_dataset = None
    nir_dataset = None

    # Create the output file path
    output_folder = QgsProject.instance().homePath()  # Get the QGIS project's home folder
    output_file_path = os.path.join(output_folder, output_path)

    # Add the output raster to QGIS
    output_layer = QgsRasterLayer(output_file_path, "BAI")
    QgsProject.instance().addMapLayer(output_layer)
    print("BAI calculated and added to QGIS.")

def ndvi_sentinel_three(raster_path):
    output_path=raster_path+'NDVI.nc'
    # Create the output file path
    output_folder = QgsProject.instance().homePath()  # Get the QGIS project's home folder
    output_file_path = os.path.join(output_folder, output_path)
    # Add the output raster to QGIS
    output_layer = QgsRasterLayer(output_file_path, "NDVI")
    QgsProject.instance().addMapLayer(output_layer)
    print("NDVI calculated and added to QGIS.")

#ndmi = (nir_band - mir_band) / (nir_band + mir_band)
def calculate_ndmi(data_source_path, nir_band_path, mir_band_path):
    # Open the NIR band raster
    nir_dataset = gdal.Open(nir_band_path)
    nir_band = nir_dataset.GetRasterBand(1)
    nir_width = nir_dataset.RasterXSize
    nir_height = nir_dataset.RasterYSize

    # Open the red band raster
    mir_dataset = gdal.Open(mir_band_path)
    mir_band = mir_dataset.GetRasterBand(1)
    mir_width = mir_dataset.RasterXSize
    mir_height = mir_dataset.RasterYSize
    # Determine the output dimensions
    output_width = min(mir_width, nir_width)
    output_height = min(mir_height, nir_height)

    if mir_width != nir_width or mir_height != nir_height:
        # Calculate the scaling factors in x and y dimensions
        # Resample all input arrays to the same size using the nearest-neighbor method
        if mir_width > nir_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            mir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            mir_dataset_resampled.SetGeoTransform(mir_dataset.GetGeoTransform())
            mir_dataset_resampled.SetProjection(mir_dataset.GetProjection())
            mir_dataset_resampled.GetRasterBand(1).Fill(0)
            mir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(mir_dataset, mir_dataset_resampled, None, None, resample_method)
            mir_band = mir_dataset_resampled.GetRasterBand(1)
        else:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            nir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1, gdal.GDT_Float32)
            nir_dataset_resampled.SetGeoTransform(nir_dataset.GetGeoTransform())
            nir_dataset_resampled.SetProjection(nir_dataset.GetProjection())
            nir_dataset_resampled.GetRasterBand(1).Fill(0)
            nir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(nir_dataset, nir_dataset_resampled, None, None, resample_method)
            nir_band = nir_dataset_resampled.GetRasterBand(1)


    # Read the band arrays as NumPy arrays
    mir_array = mir_band.ReadAsArray().astype(np.float32)
    nir_array = nir_band.ReadAsArray().astype(np.float32)

    NDMI = (nir_array - mir_array) / (nir_array + mir_array)

    # Create the output file
    output_path = data_source_path
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_path, output_width, output_height, 1, gdal.GDT_Float32)
    output_dataset.GetRasterBand(1).WriteArray(NDMI)
    output_dataset.SetGeoTransform(nir_dataset.GetGeoTransform())
    output_dataset.SetProjection(nir_dataset.GetProjection())
    output_dataset.FlushCache()
    output_dataset = None
    mir_dataset = None
    nir_dataset = None

    # Create the output file path
    output_folder = QgsProject.instance().homePath()  # Get the QGIS project's home folder
    output_file_path = os.path.join(output_folder, output_path)

    # Add the output raster to QGIS
    output_layer = QgsRasterLayer(output_file_path, "NDMI")
    QgsProject.instance().addMapLayer(output_layer)

    print("NDMI calculated and added to QGIS.")




#MSAVI = ((2 * NIR + 1) - sqrt((2 * NIR + 1)^2 - 8 * (NIR - Red))) / 2
def calculate_msavi(data_source_path, red_band_path, nir_band_path):
    # Open the red band raster
    red_dataset = gdal.Open(red_band_path)
    red_band = red_dataset.GetRasterBand(1)
    red_width = red_dataset.RasterXSize
    red_height = red_dataset.RasterYSize

    # Open the NIR band raster
    nir_dataset = gdal.Open(nir_band_path)
    nir_band = nir_dataset.GetRasterBand(1)
    nir_width = nir_dataset.RasterXSize
    nir_height = nir_dataset.RasterYSize

    # Determine the output dimensions
    output_width = min(red_width, nir_width)
    output_height = min(red_height, nir_height)

    if red_width != nir_width or red_height != nir_height:
        # Calculate the scaling factors in x and y dimensions
        # Resample all input arrays to the same size using the nearest-neighbor method
        if red_width > nir_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            red_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            red_dataset_resampled.SetGeoTransform(red_dataset.GetGeoTransform())
            red_dataset_resampled.SetProjection(red_dataset.GetProjection())
            red_dataset_resampled.GetRasterBand(1).Fill(0)
            red_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(red_dataset, red_dataset_resampled, None, None, resample_method)
            red_band = red_dataset_resampled.GetRasterBand(1)
        else:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            nir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1, gdal.GDT_Float32)
            nir_dataset_resampled.SetGeoTransform(nir_dataset.GetGeoTransform())
            nir_dataset_resampled.SetProjection(nir_dataset.GetProjection())
            nir_dataset_resampled.GetRasterBand(1).Fill(0)
            nir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(nir_dataset, nir_dataset_resampled, None, None, resample_method)
            nir_band = nir_dataset_resampled.GetRasterBand(1)


    # Read the band arrays as NumPy arrays
    red_array = red_band.ReadAsArray().astype(np.float32)
    nir_array = nir_band.ReadAsArray().astype(np.float32)
    # Calculate MSAVI = ((2 * NIR + 1) - sqrt((2 * NIR + 1)^2 - 8 * (NIR - Red))) / 2
    MSAVI = ((2 * nir_array + 1) - np.sqrt((2 * nir_array + 1) ** 2 - 8 * (nir_array - red_array))) / 2

    # Create the output file
    output_path = data_source_path
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_path, output_width, output_height, 1, gdal.GDT_Float32)
    output_dataset.GetRasterBand(1).WriteArray(MSAVI)
    output_dataset.SetGeoTransform(red_dataset.GetGeoTransform())
    output_dataset.SetProjection(red_dataset.GetProjection())
    output_dataset.FlushCache()
    output_dataset = None
    red_dataset = None
    nir_dataset = None

    # Create the output file path
    output_folder = QgsProject.instance().homePath()  # Get the QGIS project's home folder
    output_file_path = os.path.join(output_folder, output_path)

    # Add the output raster to QGIS
    output_layer = QgsRasterLayer(output_file_path, "MSAVI")
    QgsProject.instance().addMapLayer(output_layer)

    print("MSAVI calculated and added to QGIS.")

#MSI = (NIR - Red) / (NIR + Red)
def calculate_msi(data_source_path, red_band_path, nir_band_path):
    # Open the red band raster
    red_dataset = gdal.Open(red_band_path)
    red_band = red_dataset.GetRasterBand(1)
    red_width = red_dataset.RasterXSize
    red_height = red_dataset.RasterYSize

    # Open the NIR band raster
    nir_dataset = gdal.Open(nir_band_path)
    nir_band = nir_dataset.GetRasterBand(1)
    nir_width = nir_dataset.RasterXSize
    nir_height = nir_dataset.RasterYSize

    # Determine the output dimensions
    output_width = min(red_width, nir_width)
    output_height = min(red_height, nir_height)

    if red_width != nir_width or red_height != nir_height:
        # Calculate the scaling factors in x and y dimensions
        # Resample all input arrays to the same size using the nearest-neighbor method
        if red_width > nir_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            red_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            red_dataset_resampled.SetGeoTransform(red_dataset.GetGeoTransform())
            red_dataset_resampled.SetProjection(red_dataset.GetProjection())
            red_dataset_resampled.GetRasterBand(1).Fill(0)
            red_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(red_dataset, red_dataset_resampled, None, None, resample_method)
            red_band = red_dataset_resampled.GetRasterBand(1)
        else:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            nir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1, gdal.GDT_Float32)
            nir_dataset_resampled.SetGeoTransform(nir_dataset.GetGeoTransform())
            nir_dataset_resampled.SetProjection(nir_dataset.GetProjection())
            nir_dataset_resampled.GetRasterBand(1).Fill(0)
            nir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(nir_dataset, nir_dataset_resampled, None, None, resample_method)
            nir_band = nir_dataset_resampled.GetRasterBand(1)


    # Read the band arrays as NumPy arrays
    red_array = red_band.ReadAsArray().astype(np.float32)
    nir_array = nir_band.ReadAsArray().astype(np.float32)
    # MSI = (NIR - Red) / (NIR + Red)
    MSI = (nir_array-red_array)/(nir_array+red_array)

    # Create the output file
    output_path = data_source_path
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_path, output_width, output_height, 1, gdal.GDT_Float32)
    output_dataset.GetRasterBand(1).WriteArray(MSI)
    output_dataset.SetGeoTransform(red_dataset.GetGeoTransform())
    output_dataset.SetProjection(red_dataset.GetProjection())
    output_dataset.FlushCache()
    output_dataset = None
    red_dataset = None
    nir_dataset = None

    # Create the output file path
    output_folder = QgsProject.instance().homePath()  # Get the QGIS project's home folder
    output_file_path = os.path.join(output_folder, output_path)

    # Add the output raster to QGIS
    output_layer = QgsRasterLayer(output_file_path, "MSI")
    QgsProject.instance().addMapLayer(output_layer)

    print("MSI calculated and added to QGIS.")

def calculate_evi(data_source_path, red_band_path, nir_band_path, blue_band_path):
    # Open the blue band raster
    blue_dataset = gdal.Open(blue_band_path)
    blue_band = blue_dataset.GetRasterBand(1)

    # Open the red band raster
    red_dataset = gdal.Open(red_band_path)
    red_band = red_dataset.GetRasterBand(1)

    # Open the NIR band raster
    nir_dataset = gdal.Open(nir_band_path)
    nir_band = nir_dataset.GetRasterBand(1)

    # Get the dimensions of the input bands
    blue_width = blue_dataset.RasterXSize
    blue_height = blue_dataset.RasterYSize
    red_width = red_dataset.RasterXSize
    red_height = red_dataset.RasterYSize
    nir_width = nir_dataset.RasterXSize
    nir_height = nir_dataset.RasterYSize

    # Determine the output dimensions
    output_width = min(blue_width, red_width, nir_width)
    output_height = min(blue_height, red_height, nir_height)
    if red_width != nir_width or red_height != nir_height:
        # Calculate the scaling factors in x and y dimensions
        # Resample all input arrays to the same size using the nearest-neighbor method
        if red_width > nir_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour

            red_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            red_dataset_resampled.SetGeoTransform(red_dataset.GetGeoTransform())
            red_dataset_resampled.SetProjection(red_dataset.GetProjection())
            red_dataset_resampled.GetRasterBand(1).Fill(0)
            red_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(red_dataset, red_dataset_resampled, None, None, resample_method)
            red_band = red_dataset_resampled.GetRasterBand(1)

            blue_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                      gdal.GDT_Float32)
            blue_dataset_resampled.SetGeoTransform(blue_dataset.GetGeoTransform())
            blue_dataset_resampled.SetProjection(blue_dataset.GetProjection())
            blue_dataset_resampled.GetRasterBand(1).Fill(0)
            blue_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(blue_dataset, blue_dataset_resampled, None, None, resample_method)
            blue_band = blue_dataset_resampled.GetRasterBand(1)
        else:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            nir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1, gdal.GDT_Float32)
            nir_dataset_resampled.SetGeoTransform(nir_dataset.GetGeoTransform())
            nir_dataset_resampled.SetProjection(nir_dataset.GetProjection())
            nir_dataset_resampled.GetRasterBand(1).Fill(0)
            nir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(nir_dataset, nir_dataset_resampled, None, None, resample_method)
            nir_band = nir_dataset_resampled.GetRasterBand(1)

    # Read the band arrays as NumPy arrays
    red_array = red_band.ReadAsArray().astype(np.float32)
    blue_array = blue_band.ReadAsArray().astype(np.float32)
    nir_array = nir_band.ReadAsArray().astype(np.float32)

        # Calculate EVI = 2.5 * ((NIR - Red) / (NIR + 6 * Red - 7.5 * Blue + 1))
    EVI = 2.5 * ((nir_array - red_array) / ((nir_array + 6 * red_array - 7.5 * blue_array) + 1))
    # Create the output file
    output_path = data_source_path
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_path, output_width, output_height, 1, gdal.GDT_Float32)
    output_dataset.GetRasterBand(1).WriteArray(EVI)
    output_dataset.SetGeoTransform(red_dataset.GetGeoTransform())
    output_dataset.SetProjection(red_dataset.GetProjection())
    output_dataset.FlushCache()
    output_dataset = None
    red_dataset = None
    nir_dataset = None
    blue_dataset = None

    # Create the output file path
    output_folder = QgsProject.instance().homePath()  # Get the QGIS project's home folder
    output_file_path = os.path.join(output_folder, output_path)

    # Add the output raster to QGIS
    output_layer = QgsRasterLayer(output_file_path, "EVI")
    QgsProject.instance().addMapLayer(output_layer)

    print("EVI calculated and added to QGIS.")

#EVI = 2.5 * ((NIR - Red) / (NIR + 6 * Red - 7.5 * Blue + 1))
def calculate_evi2(data_source_path, red_band_path, nir_band_path, blue_band_path):
    # Open the blue band raster
    blue_dataset = gdal.Open(blue_band_path)
    blue_band = blue_dataset.GetRasterBand(1)
    blue_width = blue_dataset.RasterXSize
    blue_height = blue_dataset.RasterYSize

    # Open the red band raster
    red_dataset = gdal.Open(red_band_path)
    red_band = red_dataset.GetRasterBand(1)
    red_width = red_dataset.RasterXSize
    red_height = red_dataset.RasterYSize

    # Open the NIR band raster
    nir_dataset = gdal.Open(nir_band_path)
    nir_band = nir_dataset.GetRasterBand(1)
    nir_width = nir_dataset.RasterXSize
    nir_height = nir_dataset.RasterYSize

    output_width = min(red_width, blue_width)
    output_height = min(red_height, blue_height)

    print('width red : ',red_width,' | height red : ',red_height)
    print('width blue : ',blue_width,' | height blue : ',blue_height)
    print('width nir : ',nir_width,' | height nir : ',nir_height)

    if red_width != blue_width or red_height != blue_height:
        # Calculate the scaling factors in x and y dimensions
        # Resample all input arrays to the same size using the nearest-neighbor method
        if red_width > blue_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_Bilinear
            red_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            red_dataset_resampled.SetGeoTransform(red_dataset.GetGeoTransform())
            red_dataset_resampled.SetProjection(red_dataset.GetProjection())
            red_dataset_resampled.GetRasterBand(1).Fill(0)
            red_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(red_dataset, red_dataset_resampled, None, None, resample_method)
            red_band = red_dataset_resampled.GetRasterBand(1)
        elif red_width < blue_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            blue_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            blue_dataset_resampled.SetGeoTransform(blue_dataset.GetGeoTransform())
            blue_dataset_resampled.SetProjection(blue_dataset.GetProjection())
            blue_dataset_resampled.GetRasterBand(1).Fill(0)
            blue_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(blue_dataset, blue_dataset_resampled, None, None, resample_method)
            blue_band = blue_dataset_resampled.GetRasterBand(1)

    # Determine the output dimensions
    red_width=red_band.ReadAsArray().astype(np.float32).shape[1]
    red_height=red_band.ReadAsArray().astype(np.float32).shape[0]



    print('width red : ', red_width, ' | height red : ', red_height)
    print('width blue : ', blue_width, ' | height blue : ', blue_height)
    print('width nir : ', nir_width, ' | height nir : ', nir_height)

    # Determine the output dimensions
    output_width = min(red_width, nir_width)
    output_height = min(red_height, nir_height)
    if red_width!=nir_width or red_height!=nir_height:
        # Calculate the scaling factors in x and y dimensions
        # Resample all input arrays to the same size using the nearest-neighbor method
        if red_width > nir_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_Bilinear
            red_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            red_dataset_resampled.SetGeoTransform(red_dataset.GetGeoTransform())
            red_dataset_resampled.SetProjection(red_dataset.GetProjection())
            red_dataset_resampled.GetRasterBand(1).Fill(0)
            red_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(red_dataset, red_dataset_resampled, None, None, resample_method)
            red_band = red_dataset_resampled.GetRasterBand(1)
        elif red_width < nir_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            nir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            nir_dataset_resampled.SetGeoTransform(nir_dataset.GetGeoTransform())
            nir_dataset_resampled.SetProjection(nir_dataset.GetProjection())
            nir_dataset_resampled.GetRasterBand(1).Fill(0)
            nir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(nir_dataset, nir_dataset_resampled, None, None, resample_method)
            nir_band = nir_dataset_resampled.GetRasterBand(1)

    red_width=red_band.ReadAsArray().astype(np.float32).shape[1]
    red_height=red_band.ReadAsArray().astype(np.float32).shape[0]
    blue_width=blue_band.ReadAsArray().astype(np.float32).shape[1]
    blue_height=blue_band.ReadAsArray().astype(np.float32).shape[0]
    nir_width=nir_band.ReadAsArray().astype(np.float32).shape[1]
    nir_height=nir_band.ReadAsArray().astype(np.float32).shape[0]
    print('width red : ', red_width, ' | height red : ', red_height)
    print('width blue : ', blue_width, ' | height blue : ', blue_height)
    print('width nir : ', nir_width, ' | height nir : ', nir_height)

    output_width = min(blue_width, nir_width)
    output_height = min(blue_height, nir_height)
    if blue_width!=nir_width or blue_height!=nir_height:
        # Calculate the scaling factors in x and y dimensions
        # Resample all input arrays to the same size using the nearest-neighbor method
        if blue_width > nir_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_Bilinear
            blue_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            blue_dataset_resampled.SetGeoTransform(blue_dataset.GetGeoTransform())
            blue_dataset_resampled.SetProjection(blue_dataset.GetProjection())
            blue_dataset_resampled.GetRasterBand(1).Fill(0)
            blue_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(blue_dataset, blue_dataset_resampled, None, None, resample_method)
            blue_band = blue_dataset_resampled.GetRasterBand(1)
        elif blue_width < nir_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_Bilinear
            nir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            nir_dataset_resampled.SetGeoTransform(nir_dataset.GetGeoTransform())
            nir_dataset_resampled.SetProjection(nir_dataset.GetProjection())
            nir_dataset_resampled.GetRasterBand(1).Fill(0)
            nir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(nir_dataset, nir_dataset_resampled, None, None, resample_method)
            nir_band = nir_dataset_resampled.GetRasterBand(1)

    red_width=red_band.ReadAsArray().astype(np.float32).shape[1]
    red_height=red_band.ReadAsArray().astype(np.float32).shape[0]
    blue_width=blue_band.ReadAsArray().astype(np.float32).shape[1]
    blue_height=blue_band.ReadAsArray().astype(np.float32).shape[0]
    nir_width=nir_band.ReadAsArray().astype(np.float32).shape[1]
    nir_height=nir_band.ReadAsArray().astype(np.float32).shape[0]
    print('width red : ', red_width, ' | height red : ', red_height)
    print('width blue : ', blue_width, ' | height blue : ', blue_height)
    print('width nir : ', nir_width, ' | height nir : ', nir_height)
    # Read the band arrays as NumPy arrays
    red_array = red_band.ReadAsArray().astype(np.float32)
    nir_array = nir_band.ReadAsArray().astype(np.float32)
    blue_array = blue_band.ReadAsArray().astype(np.float32)
    red_array[np.isinf(red_array)] = 0
    red_array[np.isnan(red_array)] = 0

    nir_array[np.isinf(nir_array)] = 0
    nir_array[np.isnan(nir_array)] = 0

    blue_array[np.isinf(blue_array)] = 0
    blue_array[np.isnan(blue_array)] = 0

    # Calculate EVI = 2.5 * ((NIR - Red) / (NIR + 6 * Red - 7.5 * Blue + 1))
    EVI = 2.5 * ((nir_array - red_array) / (nir_array + 6 * red_array - 7.5 * blue_array + 1))
    # Create the output file
    output_path = data_source_path
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_path, output_width, output_height, 1, gdal.GDT_Float32)
    output_dataset.GetRasterBand(1).WriteArray(EVI)
    output_dataset.SetGeoTransform(red_dataset.GetGeoTransform())
    output_dataset.SetProjection(red_dataset.GetProjection())
    output_dataset.FlushCache()
    output_dataset = None
    red_dataset = None
    nir_dataset = None
    blue_dataset = None

    # Create the output file path
    output_folder = QgsProject.instance().homePath()  # Get the QGIS project's home folder
    output_file_path = os.path.join(output_folder, output_path)

    # Add the output raster to QGIS
    output_layer = QgsRasterLayer(output_file_path, "EVI")
    QgsProject.instance().addMapLayer(output_layer)

    print("EVI calculated and added to QGIS.")


#BSI = (SWIR1 - Red) / (SWIR1 + Red)
def calculate_bsi(data_source_path, red_band_path, swir_band_path):
    # Open the red band raster
    red_dataset = gdal.Open(red_band_path)
    red_band = red_dataset.GetRasterBand(1)
    red_width = red_dataset.RasterXSize
    red_height = red_dataset.RasterYSize

    # Open the NIR band raster
    swir_dataset = gdal.Open(swir_band_path)
    swir_band = swir_dataset.GetRasterBand(1)
    swir_width = swir_dataset.RasterXSize
    swir_height = swir_dataset.RasterYSize

    # Determine the output dimensions
    output_width = min(red_width, swir_width)
    output_height = min(red_height, swir_height)

    if red_width != swir_width or red_height != swir_height:
        # Calculate the scaling factors in x and y dimensions
        # Resample all input arrays to the same size using the nearest-neighbor method
        if red_width > swir_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            red_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            red_dataset_resampled.SetGeoTransform(red_dataset.GetGeoTransform())
            red_dataset_resampled.SetProjection(red_dataset.GetProjection())
            red_dataset_resampled.GetRasterBand(1).Fill(0)
            red_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(red_dataset, red_dataset_resampled, None, None, resample_method)
            red_band = red_dataset_resampled.GetRasterBand(1)
        else:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            swir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1, gdal.GDT_Float32)
            swir_dataset_resampled.SetGeoTransform(swir_dataset.GetGeoTransform())
            swir_dataset_resampled.SetProjection(swir_dataset.GetProjection())
            swir_dataset_resampled.GetRasterBand(1).Fill(0)
            swir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(swir_dataset, swir_dataset_resampled, None, None, resample_method)
            swir_band = swir_dataset_resampled.GetRasterBand(1)


    # Read the band arrays as NumPy arrays
    red_array = red_band.ReadAsArray().astype(np.float32)
    swir_array = swir_band.ReadAsArray().astype(np.float32)
    # Calculate BSI = (SWIR1 - Red) / (SWIR1 + Red)
    BSI = (swir_array - red_array) / (swir_array + red_array)

    # Create the output file
    output_path = data_source_path
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_path, output_width, output_height, 1, gdal.GDT_Float32)
    output_dataset.GetRasterBand(1).WriteArray(BSI)
    output_dataset.SetGeoTransform(red_dataset.GetGeoTransform())
    output_dataset.SetProjection(red_dataset.GetProjection())
    output_dataset.FlushCache()
    output_dataset = None
    red_dataset = None
    swir_dataset = None

    # Create the output file path
    output_folder = QgsProject.instance().homePath()  # Get the QGIS project's home folder
    output_file_path = os.path.join(output_folder, output_path)

    # Add the output raster to QGIS
    output_layer = QgsRasterLayer(output_file_path, "BSI")
    QgsProject.instance().addMapLayer(output_layer)

    print("BSI calculated and added to QGIS.")

#brba = (RED/SWIR)
def calculate_brba(data_source_path, red_band_path, swir_band_path):
    # Open the red band raster
    red_dataset = gdal.Open(red_band_path)
    red_band = red_dataset.GetRasterBand(1)
    red_width = red_dataset.RasterXSize
    red_height = red_dataset.RasterYSize

    # Open the NIR band raster
    swir_dataset = gdal.Open(swir_band_path)
    swir_band = swir_dataset.GetRasterBand(1)
    swir_width = swir_dataset.RasterXSize
    swir_height = swir_dataset.RasterYSize

    # Determine the output dimensions
    output_width = min(red_width, swir_width)
    output_height = min(red_height, swir_height)

    if red_width != swir_width or red_height != swir_height:
        # Calculate the scaling factors in x and y dimensions
        # Resample all input arrays to the same size using the nearest-neighbor method
        if red_width > swir_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            red_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            red_dataset_resampled.SetGeoTransform(red_dataset.GetGeoTransform())
            red_dataset_resampled.SetProjection(red_dataset.GetProjection())
            red_dataset_resampled.GetRasterBand(1).Fill(0)
            red_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(red_dataset, red_dataset_resampled, None, None, resample_method)
            red_band = red_dataset_resampled.GetRasterBand(1)
        else:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            swir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1, gdal.GDT_Float32)
            swir_dataset_resampled.SetGeoTransform(swir_dataset.GetGeoTransform())
            swir_dataset_resampled.SetProjection(swir_dataset.GetProjection())
            swir_dataset_resampled.GetRasterBand(1).Fill(0)
            swir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(swir_dataset, swir_dataset_resampled, None, None, resample_method)
            swir_band = swir_dataset_resampled.GetRasterBand(1)


    # Read the band arrays as NumPy arrays
    red_array = red_band.ReadAsArray().astype(np.float32)
    swir_array = swir_band.ReadAsArray().astype(np.float32)
    # Calculate BSI = (SWIR1 - Red) / (SWIR1 + Red)
    BRBA = (red_array) / (swir_array)

    # Create the output file
    output_path = data_source_path
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_path, output_width, output_height, 1, gdal.GDT_Float32)
    output_dataset.GetRasterBand(1).WriteArray(BRBA)
    output_dataset.SetGeoTransform(red_dataset.GetGeoTransform())
    output_dataset.SetProjection(red_dataset.GetProjection())
    output_dataset.FlushCache()
    output_dataset = None
    red_dataset = None
    swir_dataset = None

    # Create the output file path
    output_folder = QgsProject.instance().homePath()  # Get the QGIS project's home folder
    output_file_path = os.path.join(output_folder, output_path)

    # Add the output raster to QGIS
    output_layer = QgsRasterLayer(output_file_path, "BRBA")
    QgsProject.instance().addMapLayer(output_layer)

    print("BRBA calculated and added to QGIS.")

#BUI = (SWIR - NIR) / (SWIR + NIR) -( NIR - RED)/(NIR+RED)
def calculate_bui(data_source_path, swir_band_path,nir_band_path,red_band_path):

    # Open the NIR band raster
    swir_dataset = gdal.Open(swir_band_path)
    swir_band = swir_dataset.GetRasterBand(1)
    swir_width = swir_dataset.RasterXSize
    swir_height = swir_dataset.RasterYSize

    # Open the NIR band raster
    nir_dataset = gdal.Open(nir_band_path)
    nir_band = nir_dataset.GetRasterBand(1)
    nir_width = nir_dataset.RasterXSize
    nir_height = nir_dataset.RasterYSize

    # Open the red band raster
    red_dataset = gdal.Open(red_band_path)
    red_band = red_dataset.GetRasterBand(1)
    red_width = red_dataset.RasterXSize
    red_height = red_dataset.RasterYSize



    # Determine the output dimensions
    output_width = min(red_width, swir_width,nir_width)
    output_height = min(red_height, swir_height,nir_height)

    if red_width != swir_width or red_height != swir_height:
        # Calculate the scaling factors in x and y dimensions
        # Resample all input arrays to the same size using the nearest-neighbor method
        if red_width > swir_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            red_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            red_dataset_resampled.SetGeoTransform(red_dataset.GetGeoTransform())
            red_dataset_resampled.SetProjection(red_dataset.GetProjection())
            red_dataset_resampled.GetRasterBand(1).Fill(0)
            red_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(red_dataset, red_dataset_resampled, None, None, resample_method)
            red_band = red_dataset_resampled.GetRasterBand(1)

            nir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            nir_dataset_resampled.SetGeoTransform(nir_dataset.GetGeoTransform())
            nir_dataset_resampled.SetProjection(nir_dataset.GetProjection())
            nir_dataset_resampled.GetRasterBand(1).Fill(0)
            nir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(nir_dataset, nir_dataset_resampled, None, None, resample_method)
            nir_band = nir_dataset_resampled.GetRasterBand(1)
        else:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            swir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1, gdal.GDT_Float32)
            swir_dataset_resampled.SetGeoTransform(swir_dataset.GetGeoTransform())
            swir_dataset_resampled.SetProjection(swir_dataset.GetProjection())
            swir_dataset_resampled.GetRasterBand(1).Fill(0)
            swir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(swir_dataset, swir_dataset_resampled, None, None, resample_method)
            swir_band = swir_dataset_resampled.GetRasterBand(1)


    # Read the band arrays as NumPy arrays
    red_array = red_band.ReadAsArray().astype(np.float32)
    nir_array = nir_band.ReadAsArray().astype(np.float32)
    swir_array = swir_band.ReadAsArray().astype(np.float32)
    # Calculate BUI = (SWIR - NIR) / (SWIR + NIR) -( NIR - RED)/(NIR+RED)
    BUI = ((swir_array-nir_array)/(swir_array+nir_array))-((nir_array-red_array)/(nir_array+red_array))

    # Create the output file
    output_path = data_source_path
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_path, output_width, output_height, 1, gdal.GDT_Float32)
    output_dataset.GetRasterBand(1).WriteArray(BUI)
    output_dataset.SetGeoTransform(red_dataset.GetGeoTransform())
    output_dataset.SetProjection(red_dataset.GetProjection())
    output_dataset.FlushCache()
    output_dataset = None
    red_dataset = None
    swir_dataset = None

    # Create the output file path
    output_folder = QgsProject.instance().homePath()  # Get the QGIS project's home folder
    output_file_path = os.path.join(output_folder, output_path)

    # Add the output raster to QGIS
    output_layer = QgsRasterLayer(output_file_path, "BUI")
    QgsProject.instance().addMapLayer(output_layer)
    print("BUI calculated and added to QGIS.")

#CBI = (NIR * (SWIR / Red)) - 1
def calculate_cbi(data_source_path, swir_band_path,nir_band_path,red_band_path):

    # Open the NIR band raster
    swir_dataset = gdal.Open(swir_band_path)
    swir_band = swir_dataset.GetRasterBand(1)
    swir_width = swir_dataset.RasterXSize
    swir_height = swir_dataset.RasterYSize

    # Open the NIR band raster
    nir_dataset = gdal.Open(nir_band_path)
    nir_band = nir_dataset.GetRasterBand(1)
    nir_width = nir_dataset.RasterXSize
    nir_height = nir_dataset.RasterYSize

    # Open the red band raster
    red_dataset = gdal.Open(red_band_path)
    red_band = red_dataset.GetRasterBand(1)
    red_width = red_dataset.RasterXSize
    red_height = red_dataset.RasterYSize



    # Determine the output dimensions
    output_width = min(red_width, swir_width,nir_width)
    output_height = min(red_height, swir_height,nir_height)

    if red_width != swir_width or red_height != swir_height:
        # Calculate the scaling factors in x and y dimensions
        # Resample all input arrays to the same size using the nearest-neighbor method
        if red_width > swir_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            red_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            red_dataset_resampled.SetGeoTransform(red_dataset.GetGeoTransform())
            red_dataset_resampled.SetProjection(red_dataset.GetProjection())
            red_dataset_resampled.GetRasterBand(1).Fill(0)
            red_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(red_dataset, red_dataset_resampled, None, None, resample_method)
            red_band = red_dataset_resampled.GetRasterBand(1)

            nir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            nir_dataset_resampled.SetGeoTransform(nir_dataset.GetGeoTransform())
            nir_dataset_resampled.SetProjection(nir_dataset.GetProjection())
            nir_dataset_resampled.GetRasterBand(1).Fill(0)
            nir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(nir_dataset, nir_dataset_resampled, None, None, resample_method)
            nir_band = nir_dataset_resampled.GetRasterBand(1)
        else:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            swir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1, gdal.GDT_Float32)
            swir_dataset_resampled.SetGeoTransform(swir_dataset.GetGeoTransform())
            swir_dataset_resampled.SetProjection(swir_dataset.GetProjection())
            swir_dataset_resampled.GetRasterBand(1).Fill(0)
            swir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(swir_dataset, swir_dataset_resampled, None, None, resample_method)
            swir_band = swir_dataset_resampled.GetRasterBand(1)


    # Read the band arrays as NumPy arrays
    red_array = red_band.ReadAsArray().astype(np.float32)
    nir_array = nir_band.ReadAsArray().astype(np.float32)
    swir_array = swir_band.ReadAsArray().astype(np.float32)
    # Calculate CBI = (NIR * (SWIR / Red)) - 1
    CBI = (nir_array* (swir_array/ red_array)) -1

    # Create the output file
    output_path = data_source_path
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_path, output_width, output_height, 1, gdal.GDT_Float32)
    output_dataset.GetRasterBand(1).WriteArray(CBI)
    output_dataset.SetGeoTransform(red_dataset.GetGeoTransform())
    output_dataset.SetProjection(red_dataset.GetProjection())
    output_dataset.FlushCache()
    output_dataset = None
    red_dataset = None
    swir_dataset = None

    # Create the output file path
    output_folder = QgsProject.instance().homePath()  # Get the QGIS project's home folder
    output_file_path = os.path.join(output_folder, output_path)

    # Add the output raster to QGIS
    output_layer = QgsRasterLayer(output_file_path, "CBI")
    QgsProject.instance().addMapLayer(output_layer)
    print("CBI calculated and added to QGIS.")
#NBR = (NIR - SWIR) / (NIR + SWIR)
def calculate_nbr(data_source_path, swir_band_path, nir_band_path):
    # swir band 12
    # Open the red band raster
    nir_dataset = gdal.Open(nir_band_path)
    nir_band = nir_dataset.GetRasterBand(1)
    nir_width = nir_dataset.RasterXSize
    nir_height = nir_dataset.RasterYSize

    # Open the NIR band raster
    swir_dataset = gdal.Open(swir_band_path)
    swir_band = swir_dataset.GetRasterBand(1)
    swir_width = swir_dataset.RasterXSize
    swir_height = swir_dataset.RasterYSize

    # Determine the output dimensions
    output_width = min(nir_width, swir_width)
    output_height = min(nir_height, swir_height)

    if nir_width != swir_width or nir_height != swir_height:
        # Calculate the scaling factors in x and y dimensions
        # Resample all input arrays to the same size using the nearest-neighbor method
        if nir_width > swir_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            nir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            nir_dataset_resampled.SetGeoTransform(nir_dataset.GetGeoTransform())
            nir_dataset_resampled.SetProjection(nir_dataset.GetProjection())
            nir_dataset_resampled.GetRasterBand(1).Fill(0)
            nir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(nir_dataset, nir_dataset_resampled, None, None, resample_method)
            nir_band = nir_dataset_resampled.GetRasterBand(1)
        else:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            swir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1, gdal.GDT_Float32)
            swir_dataset_resampled.SetGeoTransform(swir_dataset.GetGeoTransform())
            swir_dataset_resampled.SetProjection(swir_dataset.GetProjection())
            swir_dataset_resampled.GetRasterBand(1).Fill(0)
            swir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(swir_dataset, swir_dataset_resampled, None, None, resample_method)
            swir_band = swir_dataset_resampled.GetRasterBand(1)


    # Read the band arrays as NumPy arrays
    nir_array = nir_band.ReadAsArray().astype(np.float32)
    swir_array = swir_band.ReadAsArray().astype(np.float32)
    # Calculate NBR = (NIR - SWIR) / (NIR + SWIR)
    NBR = (nir_array - swir_array) / (nir_array + swir_array)
    # Create the output file
    output_path = data_source_path
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_path, output_width, output_height, 1, gdal.GDT_Float32)
    output_dataset.GetRasterBand(1).WriteArray(NBR)
    output_dataset.SetGeoTransform(nir_dataset.GetGeoTransform())
    output_dataset.SetProjection(nir_dataset.GetProjection())
    output_dataset.FlushCache()
    output_dataset = None
    swir_dataset = None
    nir_dataset = None

    # Create the output file path
    output_folder = QgsProject.instance().homePath()  # Get the QGIS project's home folder
    output_file_path = os.path.join(output_folder, output_path)

    # Add the output raster to QGIS
    output_layer = QgsRasterLayer(output_file_path, "NBR")
    QgsProject.instance().addMapLayer(output_layer)

    print("NBR calculated and added to QGIS.")

#LSWI = (NIR - SWIR) / (NIR + SWIR)
def calculate_lswi(data_source_path, swir_band_path, nir_band_path):
    #swir band 11
    # Open the red band raster
    nir_dataset = gdal.Open(nir_band_path)
    nir_band = nir_dataset.GetRasterBand(1)
    nir_width = nir_dataset.RasterXSize
    nir_height = nir_dataset.RasterYSize

    # Open the NIR band raster
    swir_dataset = gdal.Open(swir_band_path)
    swir_band = swir_dataset.GetRasterBand(1)
    swir_width = swir_dataset.RasterXSize
    swir_height = swir_dataset.RasterYSize

    # Determine the output dimensions
    output_width = min(nir_width, swir_width)
    output_height = min(nir_height, swir_height)

    if nir_width != swir_width or nir_height != swir_height:
        # Calculate the scaling factors in x and y dimensions
        # Resample all input arrays to the same size using the nearest-neighbor method
        if nir_width > swir_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            nir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            nir_dataset_resampled.SetGeoTransform(nir_dataset.GetGeoTransform())
            nir_dataset_resampled.SetProjection(nir_dataset.GetProjection())
            nir_dataset_resampled.GetRasterBand(1).Fill(0)
            nir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(nir_dataset, nir_dataset_resampled, None, None, resample_method)
            nir_band = nir_dataset_resampled.GetRasterBand(1)
        else:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            swir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1, gdal.GDT_Float32)
            swir_dataset_resampled.SetGeoTransform(swir_dataset.GetGeoTransform())
            swir_dataset_resampled.SetProjection(swir_dataset.GetProjection())
            swir_dataset_resampled.GetRasterBand(1).Fill(0)
            swir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(swir_dataset, swir_dataset_resampled, None, None, resample_method)
            swir_band = swir_dataset_resampled.GetRasterBand(1)


    # Read the band arrays as NumPy arrays
    nir_array = nir_band.ReadAsArray().astype(np.float32)
    swir_array = swir_band.ReadAsArray().astype(np.float32)
    # Calculate LSWI = (NIR - SWIR) / (NIR + SWIR)
    LSWI = (nir_array - swir_array) / (nir_array + swir_array)
    # Create the output file
    output_path = data_source_path
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_path, output_width, output_height, 1, gdal.GDT_Float32)
    output_dataset.GetRasterBand(1).WriteArray(LSWI)
    output_dataset.SetGeoTransform(nir_dataset.GetGeoTransform())
    output_dataset.SetProjection(nir_dataset.GetProjection())
    output_dataset.FlushCache()
    output_dataset = None
    swir_dataset = None
    nir_dataset = None

    # Create the output file path
    output_folder = QgsProject.instance().homePath()  # Get the QGIS project's home folder
    output_file_path = os.path.join(output_folder, output_path)

    # Add the output raster to QGIS
    output_layer = QgsRasterLayer(output_file_path, "LSWI")
    QgsProject.instance().addMapLayer(output_layer)
    print("LSWI calculated and added to QGIS.")




#NDBI = (SWIR - NIR) / (SWIR + NIR)
def calculate_ndbi(data_source_path, swir_band_path, nir_band_path):
    # Open the red band raster
    nir_dataset = gdal.Open(nir_band_path)
    nir_band = nir_dataset.GetRasterBand(1)
    nir_width = nir_dataset.RasterXSize
    nir_height = nir_dataset.RasterYSize

    # Open the NIR band raster
    swir_dataset = gdal.Open(swir_band_path)
    swir_band = swir_dataset.GetRasterBand(1)
    swir_width = swir_dataset.RasterXSize
    swir_height = swir_dataset.RasterYSize

    # Determine the output dimensions
    output_width = min(nir_width, swir_width)
    output_height = min(nir_height, swir_height)

    if nir_width != swir_width or nir_height != swir_height:
        # Calculate the scaling factors in x and y dimensions
        # Resample all input arrays to the same size using the nearest-neighbor method
        if nir_width > swir_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            nir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            nir_dataset_resampled.SetGeoTransform(nir_dataset.GetGeoTransform())
            nir_dataset_resampled.SetProjection(nir_dataset.GetProjection())
            nir_dataset_resampled.GetRasterBand(1).Fill(0)
            nir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(nir_dataset, nir_dataset_resampled, None, None, resample_method)
            nir_band = nir_dataset_resampled.GetRasterBand(1)
        else:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            swir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                        gdal.GDT_Float32)
            swir_dataset_resampled.SetGeoTransform(swir_dataset.GetGeoTransform())
            swir_dataset_resampled.SetProjection(swir_dataset.GetProjection())
            swir_dataset_resampled.GetRasterBand(1).Fill(0)
            swir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(swir_dataset, swir_dataset_resampled, None, None, resample_method)
            swir_band = swir_dataset_resampled.GetRasterBand(1)

    # Read the band arrays as NumPy arrays
    nir_array = nir_band.ReadAsArray().astype(np.float32)
    swir_array = swir_band.ReadAsArray().astype(np.float32)
    # Calculate NDBI = (SWIR - NIR) / (SWIR + NIR)
    NDBI = (swir_array - nir_array) / (swir_array + nir_array)

    # Create the output file
    output_path = data_source_path
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_path, output_width, output_height, 1, gdal.GDT_Float32)
    output_dataset.GetRasterBand(1).WriteArray(NDBI)
    output_dataset.SetGeoTransform(swir_dataset.GetGeoTransform())
    output_dataset.SetProjection(swir_dataset.GetProjection())
    output_dataset.FlushCache()
    output_dataset = None
    swir_dataset = None
    nir_dataset = None

    # Create the output file path
    output_folder = QgsProject.instance().homePath()  # Get the QGIS project's home folder
    output_file_path = os.path.join(output_folder, output_path)

    # Add the output raster to QGIS
    output_layer = QgsRasterLayer(output_file_path, "NDBI")
    QgsProject.instance().addMapLayer(output_layer)

    print("NDBI calculated and added to QGIS.")


#HS = sqrt((RED - GREEN)^2 + (RED - BLUE)(GREEN - BLUE))
def calculate_hs(data_source_path, red_band, green_band, blue_band):
    # Open the red band raster
    red_dataset = gdal.Open(red_band)
    red_band = red_dataset.GetRasterBand(1).ReadAsArray().astype(np.float32)
    red_width = red_dataset.RasterXSize
    red_height = red_dataset.RasterYSize

    # Open the green band raster
    green_dataset = gdal.Open(green_band)
    green_band = green_dataset.GetRasterBand(1).ReadAsArray().astype(np.float32)
    green_width = green_dataset.RasterXSize
    green_height = green_dataset.RasterYSize

    # Open the blue band raster
    blue_dataset = gdal.Open(blue_band)
    blue_band = blue_dataset.GetRasterBand(1).ReadAsArray().astype(np.float32)
    blue_width = blue_dataset.RasterXSize
    blue_height = blue_dataset.RasterYSize
    # Determine the output dimensions
    output_width = min(red_width, green_width, blue_width)
    output_height = min(red_height, green_height, blue_height)

    # Calculate HS = sqrt((RED - GREEN)^2 + (RED - BLUE)(GREEN - BLUE))
    HS = np.sqrt((red_band - green_band) ** 2 + (red_band - blue_band) * (green_band - blue_band))

    # Create the output file
    output_path = data_source_path
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_path, output_width, output_height, 1, gdal.GDT_Float32)
    output_dataset.GetRasterBand(1).WriteArray(HS)
    output_dataset.SetGeoTransform(red_dataset.GetGeoTransform())
    output_dataset.SetProjection(red_dataset.GetProjection())
    output_dataset.FlushCache()
    output_dataset = None
    red_dataset = None
    green_dataset = None
    blue_dataset = None

    # Create the output file path
    output_folder = QgsProject.instance().homePath()  # Get the QGIS project's home folder
    output_file_path = os.path.join(output_folder, output_path)

    # Add the output raster to QGIS
    output_layer = QgsRasterLayer(output_file_path, "HS")
    QgsProject.instance().addMapLayer(output_layer)

    print("HS calculated and added to QGIS.")


#UTFVI = (NIR - R) / (NIR + RED - BLUE)
def calculate_utfvi(data_source_path, red_band, nir_band, blue_band):
    # Open the red band raster
    red_dataset = gdal.Open(red_band)
    red_band = red_dataset.GetRasterBand(1).ReadAsArray().astype(np.float32)
    red_width = red_dataset.RasterXSize
    red_height = red_dataset.RasterYSize

    # Open the NIR band raster
    nir_dataset = gdal.Open(nir_band)
    nir_band = nir_dataset.GetRasterBand(1).ReadAsArray().astype(np.float32)
    nir_width = nir_dataset.RasterXSize
    nir_height = nir_dataset.RasterYSize

    # Open the blue band raster
    blue_dataset = gdal.Open(blue_band)
    blue_band = blue_dataset.GetRasterBand(1).ReadAsArray().astype(np.float32)
    blue_width = blue_dataset.RasterXSize
    blue_height = blue_dataset.RasterYSize
    # Determine the output dimensions
    output_width = min(red_width, nir_width, blue_width)
    output_height = min(red_height, nir_height, blue_height)
    if red_width > nir_width:
        target_width = output_width
        target_height = output_height
        resample_method = gdal.GRA_NearestNeighbour
        red_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                   gdal.GDT_Float32)
        red_dataset_resampled.SetGeoTransform(red_dataset.GetGeoTransform())
        red_dataset_resampled.SetProjection(red_dataset.GetProjection())
        red_dataset_resampled.GetRasterBand(1).Fill(0)
        red_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
        gdal.ReprojectImage(red_dataset, red_dataset_resampled, None, None, resample_method)
        red_band = red_dataset_resampled.GetRasterBand(1)
    else:
        target_width = output_width
        target_height = output_height
        resample_method = gdal.GRA_NearestNeighbour
        nir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1, gdal.GDT_Float32)
        nir_dataset_resampled.SetGeoTransform(nir_dataset.GetGeoTransform())
        nir_dataset_resampled.SetProjection(nir_dataset.GetProjection())
        nir_dataset_resampled.GetRasterBand(1).Fill(0)
        nir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
        gdal.ReprojectImage(nir_dataset, nir_dataset_resampled, None, None, resample_method)
        nir_band = nir_dataset_resampled.GetRasterBand(1)

    # Calculate UTFVI = (NIR - R) / (NIR + RED - BLUE)
    UTFVI = (nir_band - red_band) / (nir_band + red_band - blue_band)

    # Create the output file
    output_path = data_source_path
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_path, output_width, output_height, 1, gdal.GDT_Float32)
    output_dataset.GetRasterBand(1).WriteArray(UTFVI)
    output_dataset.SetGeoTransform(red_dataset.GetGeoTransform())
    output_dataset.SetProjection(red_dataset.GetProjection())
    output_dataset.FlushCache()
    output_dataset = None
    red_dataset = None
    nir_dataset = None
    blue_dataset = None

    # Create the output file path
    output_folder = QgsProject.instance().homePath()  # Get the QGIS project's home folder
    output_file_path = os.path.join(output_folder, output_path)

    # Add the output raster to QGIS
    output_layer = QgsRasterLayer(output_file_path, "UTFVI")
    QgsProject.instance().addMapLayer(output_layer)

    print("UTFVI calculated and added to QGIS.")


#NDWI = (Green - NIR) / (Green + NIR)
def calculate_ndwi(data_source_path, green_band_path, nir_band_path):
    # Open the red band raster
    nir_dataset = gdal.Open(nir_band_path)
    nir_band = nir_dataset.GetRasterBand(1)
    nir_width = nir_dataset.RasterXSize
    nir_height = nir_dataset.RasterYSize

    # Open the NIR band raster
    green_dataset = gdal.Open(green_band_path)
    green_band = green_dataset.GetRasterBand(1)
    green_width = green_dataset.RasterXSize
    green_height = green_dataset.RasterYSize

    # Determine the output dimensions
    output_width = min(nir_width, green_width)
    output_height = min(nir_height, green_height)

    if nir_width != green_width or nir_height != green_height:
        # Calculate the scaling factors in x and y dimensions
        # Resample all input arrays to the same size using the nearest-neighbor method
        if nir_width > green_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            nir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            nir_dataset_resampled.SetGeoTransform(nir_dataset.GetGeoTransform())
            nir_dataset_resampled.SetProjection(nir_dataset.GetProjection())
            nir_dataset_resampled.GetRasterBand(1).Fill(0)
            nir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(nir_dataset, nir_dataset_resampled, None, None, resample_method)
            nir_band = nir_dataset_resampled.GetRasterBand(1)
        else:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            green_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                        gdal.GDT_Float32)
            green_dataset_resampled.SetGeoTransform(green_dataset.GetGeoTransform())
            green_dataset_resampled.SetProjection(green_dataset.GetProjection())
            green_dataset_resampled.GetRasterBand(1).Fill(0)
            green_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(green_dataset, green_dataset_resampled, None, None, resample_method)
            green_band = green_dataset_resampled.GetRasterBand(1)

    # Read the band arrays as NumPy arrays
    nir_array = nir_band.ReadAsArray().astype(np.float32)
    green_array = green_band.ReadAsArray().astype(np.float32)
    # Calculate NDWI = (Green - NIR) / (Green + NIR)
    NDWI = (green_array - nir_array) / (green_array + nir_array)

    # Create the output file
    output_path = data_source_path
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_path, output_width, output_height, 1, gdal.GDT_Float32)
    output_dataset.GetRasterBand(1).WriteArray(NDWI)
    output_dataset.SetGeoTransform(green_dataset.GetGeoTransform())
    output_dataset.SetProjection(green_dataset.GetProjection())
    output_dataset.FlushCache()
    output_dataset = None
    green_dataset = None
    nir_dataset = None

    # Create the output file path
    output_folder = QgsProject.instance().homePath()  # Get the QGIS project's home folder
    output_file_path = os.path.join(output_folder, output_path)

    # Add the output raster to QGIS
    output_layer = QgsRasterLayer(output_file_path, "NDWI")
    QgsProject.instance().addMapLayer(output_layer)

    print("NDWI calculated and added to QGIS.")

#GNDVI=(NIR - GREEN)/( NIR + GREEN )
def calculate_gndvi(data_source_path, green_band_path, nir_band_path):
    # Open the red band raster
    nir_dataset = gdal.Open(nir_band_path)
    nir_band = nir_dataset.GetRasterBand(1)
    nir_width = nir_dataset.RasterXSize
    nir_height = nir_dataset.RasterYSize

    # Open the NIR band raster
    green_dataset = gdal.Open(green_band_path)
    green_band = green_dataset.GetRasterBand(1)
    green_width = green_dataset.RasterXSize
    green_height = green_dataset.RasterYSize

    # Determine the output dimensions
    output_width = min(nir_width, green_width)
    output_height = min(nir_height, green_height)

    if nir_width != green_width or nir_height != green_height:
        # Calculate the scaling factors in x and y dimensions
        # Resample all input arrays to the same size using the nearest-neighbor method
        if nir_width > green_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            nir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            nir_dataset_resampled.SetGeoTransform(nir_dataset.GetGeoTransform())
            nir_dataset_resampled.SetProjection(nir_dataset.GetProjection())
            nir_dataset_resampled.GetRasterBand(1).Fill(0)
            nir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(nir_dataset, nir_dataset_resampled, None, None, resample_method)
            nir_band = nir_dataset_resampled.GetRasterBand(1)
        else:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            green_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                        gdal.GDT_Float32)
            green_dataset_resampled.SetGeoTransform(green_dataset.GetGeoTransform())
            green_dataset_resampled.SetProjection(green_dataset.GetProjection())
            green_dataset_resampled.GetRasterBand(1).Fill(0)
            green_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(green_dataset, green_dataset_resampled, None, None, resample_method)
            green_band = green_dataset_resampled.GetRasterBand(1)

    # Read the band arrays as NumPy arrays
    nir_array = nir_band.ReadAsArray().astype(np.float32)
    green_array = green_band.ReadAsArray().astype(np.float32)
    # Calculate GNDVI= ( NIR - GREEN ) / ( NIR + GREEN )
    GNDVI = (nir_array-green_array) / (nir_array+green_array)

    # Create the output file
    output_path = data_source_path
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_path, output_width, output_height, 1, gdal.GDT_Float32)
    output_dataset.GetRasterBand(1).WriteArray(GNDVI)
    output_dataset.SetGeoTransform(green_dataset.GetGeoTransform())
    output_dataset.SetProjection(green_dataset.GetProjection())
    output_dataset.FlushCache()
    output_dataset = None
    green_dataset = None
    nir_dataset = None

    # Create the output file path
    output_folder = QgsProject.instance().homePath()  # Get the QGIS project's home folder
    output_file_path = os.path.join(output_folder, output_path)

    # Add the output raster to QGIS
    output_layer = QgsRasterLayer(output_file_path, "GNDVI")
    QgsProject.instance().addMapLayer(output_layer)

    print("GNDVI calculated and added to QGIS.")

#MNDWI = (Green - SWIR) / (Green + SWIR)
def calculate_mndwi(data_source_path, green_band_path, swir_band_path):
    # Open the red band raster
    swir_dataset = gdal.Open(swir_band_path)
    swir_band = swir_dataset.GetRasterBand(1)
    swir_width = swir_dataset.RasterXSize
    swir_height = swir_dataset.RasterYSize

    # Open the NIR band raster
    green_dataset = gdal.Open(green_band_path)
    green_band = green_dataset.GetRasterBand(1)
    green_width = green_dataset.RasterXSize
    green_height = green_dataset.RasterYSize

    # Determine the output dimensions
    output_width = min(swir_width, green_width)
    output_height = min(swir_height, green_height)

    if swir_width != green_width or swir_height != green_height:
        # Calculate the scaling factors in x and y dimensions
        # Resample all input arrays to the same size using the nearest-neighbor method
        if output_width > green_width:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            swir_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                       gdal.GDT_Float32)
            swir_dataset_resampled.SetGeoTransform(swir_dataset.GetGeoTransform())
            swir_dataset_resampled.SetProjection(swir_dataset.GetProjection())
            swir_dataset_resampled.GetRasterBand(1).Fill(0)
            swir_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(swir_dataset, swir_dataset_resampled, None, None, resample_method)
            swir_band = swir_dataset_resampled.GetRasterBand(1)
        else:
            target_width = output_width
            target_height = output_height
            resample_method = gdal.GRA_NearestNeighbour
            green_dataset_resampled = gdal.GetDriverByName('MEM').Create('', target_width, target_height, 1,
                                                                         gdal.GDT_Float32)
            green_dataset_resampled.SetGeoTransform(green_dataset.GetGeoTransform())
            green_dataset_resampled.SetProjection(green_dataset.GetProjection())
            green_dataset_resampled.GetRasterBand(1).Fill(0)
            green_dataset_resampled.GetRasterBand(1).SetNoDataValue(0)
            gdal.ReprojectImage(green_dataset, green_dataset_resampled, None, None, resample_method)
            green_band = green_dataset_resampled.GetRasterBand(1)

    # Read the band arrays as NumPy arrays
    swir_array = swir_band.ReadAsArray().astype(np.float32)
    green_array = green_band.ReadAsArray().astype(np.float32)
    # Calculate MNDWI = (Green - SWIR) / (Green + SWIR)
    MNDWI = (green_array - swir_array) / (green_array + swir_array)

    # Create the output file
    output_path = data_source_path
    driver = gdal.GetDriverByName('GTiff')
    output_dataset = driver.Create(output_path, output_width, output_height, 1, gdal.GDT_Float32)
    output_dataset.GetRasterBand(1).WriteArray(MNDWI)
    output_dataset.SetGeoTransform(green_dataset.GetGeoTransform())
    output_dataset.SetProjection(green_dataset.GetProjection())
    output_dataset.FlushCache()
    output_dataset = None
    green_dataset = None
    swir_dataset = None

    # Create the output file path
    output_folder = QgsProject.instance().homePath()  # Get the QGIS project's home folder
    output_file_path = os.path.join(output_folder, output_path)

    # Add the output raster to QGIS
    output_layer = QgsRasterLayer(output_file_path, "MNDWI")
    QgsProject.instance().addMapLayer(output_layer)

    print("MNDWI calculated and added to QGIS.")

