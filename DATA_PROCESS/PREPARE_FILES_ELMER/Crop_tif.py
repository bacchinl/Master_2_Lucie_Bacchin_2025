from osgeo import gdal


years = [2022]
#years=range(1992,2022,1)
for year in years:

    # Charger le fichier TIFF
    #input_tiff = f"../DATA/BACKGROUND_IMAGES/rema_mosaic_100m_v2.0_browse.tif"
    #output_tiff = f"../DATA/BACKGROUND_IMAGES/rema_mosaic_100m_v2.0_browse_crop.tif"

    #input_tiff = f"../DATA/Mask/mask_full_PIG_{year}.tif"
    #output_tiff = f"../DATA/Mask/mask_PIG_{year}_crop.tif"

    input_tiff = f"../DATA/BACKGROUND_IMAGES/background_2010.tif"
    output_tiff = f"../DATA/BACKGROUND_IMAGES/background_2010_crop.tif"



    # Spécifiez les coordonnées de la zone d'intérêt
    xmin, xmax = -1694617, -1550000   # Exemple
    ymin, ymax = -371235, -250000

    #dataset = gdal.Open(input_tiff)
    #geo_transform = dataset.GetGeoTransform()

    #x_min = geo_transform[0]
    #y_max = geo_transform[3]
    #x_max = x_min + dataset.RasterXSize * geo_transform[1]
    #y_min = y_max + dataset.RasterYSize * geo_transform[5]

    #print(f"Limites du raster : x = [{x_min}, {x_max}], y = [{y_min}, {y_max}]")

    # Utiliser GDAL pour lire et extraire la zone
    gdal.Translate(
        output_tiff,
        input_tiff,
        projWin=[xmin, ymax, xmax, ymin]  # Attention à l'ordre
    )
print("Découpage terminé !")

