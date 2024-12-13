""" In order for these scripts to work, the LRS needs to be prepared with the following
    steps applied to both master and overlap versions:
        - Reduce lrs to only routes that will are near the input conflation dataset.
            This will reduce processing time
        - Convert lrs to singlepart geometries
        - Add a field called ID with a unique id for each record
        - Save as a shapefile
"""

import arcpy
import os

arcpy.env.overwriteOutput = True

# Input LRS paths
input_master_lrs = r'C:\Users\daniel.fourquet\Documents\ArcGIS\LRS_23.gdb\SDE_VDOT_RTE_MASTER_LRS_DY'
input_overlap_lrs = r'C:\Users\daniel.fourquet\Documents\ArcGIS\LRS_23.gdb\SDE_VDOT_RTE_OVERLAP_LRS_DY'

output_lrs_path = r'C:\Users\daniel.fourquet\Documents\Tasks\TMC Conflation 2025\Data'
output_master_lrs_filename = 'master_lrs.shp'  # eg master_lrs.shp
output_overlap_lrs_filename = 'overlap_lrs.shp'  # eg overlap_lrs.shp

lrs_versions = [(input_master_lrs, output_master_lrs_filename), (input_overlap_lrs, output_overlap_lrs_filename)]

# Input conflation layer (eg TMCs)
input_conflation_layer = r'C:\Users\daniel.fourquet\Documents\Tasks\TMC Conflation 2025\NPMRDS\data\Virginia_NPMRDS.shp'
lyr_conflation_layer = 'lyr_conflation_layer'
arcpy.MakeFeatureLayer_management(input_conflation_layer, lyr_conflation_layer)

# Output LRS SHP path

for lrs, output_filename in lrs_versions:
    print(lrs)
    # Make lrs layer
    print('    Making Feature Layer')
    arcpy.MakeFeatureLayer_management(lrs, 'lrs')

    # Select nearby routes
    print('    Making Selection')
    arcpy.SelectLayerByLocation_management('lrs', 'WITHIN_A_DISTANCE', lyr_conflation_layer, '200 METERS')

    # Make single-part geometries as SHP
    print('    Single-Part Geometry')
    lrs_singlepart = os.path.join('memory', 'lrs_singlepart')
    arcpy.MultipartToSinglepart_management('lrs', lrs_singlepart)
    
    print('    Save as SHP')
    output_path = os.path.join(output_lrs_path, output_filename)
    arcpy.FeatureClassToFeatureClass_conversion(lrs_singlepart, output_lrs_path, output_filename)

    # Add ID field
    print('    Adding ID Field')
    arcpy.AddField_management(output_path, 'ID', 'LONG')
    with arcpy.da.UpdateCursor(output_path, ('FID', 'ID')) as cur:
        for row in cur:
            row[1] = row[0]
            cur.updateRow(row)