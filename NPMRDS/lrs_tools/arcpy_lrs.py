import arcpy

def get_line_mp(inputPolyline, lrs, rte_nm, check_for_multipart=False):
    """ Locates the begin and end MP values of an input line along the LRS
        ** The spatial reference of the input must match the spatial reference
           of the lrs! **
    Input:
        inputPolyline - an arcpy Polyline object
        lrs - a reference to the lrs layer
        rte_nm - the lrs rte_nm that the polyline will be placed on
    Output:
        (beginMP, endMP)
    """

    try:
        # Get the geometry for the LRS route
        RouteGeom = None
        if lrs.getSelectionSet():
            arcpy.management.SelectLayerByAttribute(lrs, 'CLEAR_SELECTION')
        with arcpy.da.SearchCursor(lrs, "SHAPE@", "RTE_NM = '{}'".format(rte_nm)) as cur:
            for row in cur:
                RouteGeom = row[0]

        if not RouteGeom:
            print(f'Route "{rte_nm}" not found')
            return None, None

        # Check for multipart geometry.  If multipart, find closest part to
        # ensure that the correct MP is returned
        if check_for_multipart and RouteGeom.isMultipart:
            # Get list of parts
            parts = [arcpy.Polyline(RouteGeom[i], has_m=True) for i in range(RouteGeom.partCount)]

            # Get input polyline midpoint
            midPoint = inputPolyline.positionAlongLine(0.5, use_percentage=True)

            # Get distances from inputPolyline's mid-point to each route part
            partDists = {midPoint.distanceTo(part):part for part in parts}

            # Replace RouteGeom with closest polyline part
        #     RouteGeom = partDists[min(partDists)]

        def get_mp_from_point(route, point):
            """ Returns the m-value along the input route geometry
                given an input point geometry """
            rteMeasure = route.measureOnLine(point)
            rtePosition = route.positionAlongLine(rteMeasure)
            mp = rtePosition.firstPoint.M
            return mp

        beginPt = inputPolyline.firstPoint
        beginMP = get_mp_from_point(RouteGeom, beginPt)

        endPt = inputPolyline.lastPoint
        endMP = get_mp_from_point(RouteGeom, endPt)

        return round(beginMP, 3), round(endMP, 3)

    except Exception as e:
        print(e)
        return None, None

        
def update_line_events_known_rte_nm(layer, lrs, rte_nm_field, begin_msr_update, end_msr_update):
    """ Updates the input layer with updated measures based on the input lrs.
        the measures will be updated in the begin_msr_update and end_msr_update
        fields
    
        The input layer and lrs must be in the same projection.
    """

    def get_layer(name):
        prj = arcpy.mp.ArcGISProject("CURRENT")
        map = prj.listMaps()[0]
        layer = map.listLayers(name)[0]
        return layer

    if type(layer) == str:
        layer = get_layer(layer)
    if type(lrs) == str:
        lrs = get_layer(lrs)

    with arcpy.da.UpdateCursor(layer, [rte_nm_field, begin_msr_update, end_msr_update, 'SHAPE@']) as cur:
        for row in cur:
            begin_msr, end_msr = get_line_mp(row[-1], lrs, row[0])
            row[1] = begin_msr
            row[2] = end_msr
            cur.updateRow(row)
