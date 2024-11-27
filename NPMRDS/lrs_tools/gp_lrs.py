import geopandas as gp
import shapefile
from shapely.geometry import Point, LineString
import shapely
import math

class LRS:
    def __init__(self, lrs_shp_path, filter=False, ramps=2, epsg=3968):
        self.geodataframe = self.open_lrs(lrs_shp_path, filter, ramps, epsg)
        self.mDict = self.load_m_values(lrs_shp_path) # List of M-values
        self.spatial_index = self.geodataframe.sindex


    def open_lrs(self, lrs_shp_path, filter, ramps=None, epsg=3968):
        """
        Opens the input lrs shapefile as a GeoPandas GeoDataFrame

        Args:
            lrs_shp_path: Path to the LRS in a shapefile format.
            filter: If true, only include NB, EB, and PA routes and IS routes
                in both directions.
            ramps: 0 = ramps excluded.  1 = only include ramps.  Ignored if no filter.
            epsg: Desired output projection.  Virginia Lambert is default.

        Returns:
            Geopandas LRS object
        """

        lrs = gp.read_file(lrs_shp_path)
        if filter:
            include_interstates = (lrs['RTE_NM'].str.startswith('R-VA') & lrs['RTE_NM'].str.contains('IS'))
            include_ramps = lrs['RTE_NM'].str.contains('RMP')
            exclude_ramps = ~lrs['RTE_NM'].str.contains('RMP')
            include_NB = (lrs['RTE_NM'].str.startswith('R-VA') & lrs['RTE_NM'].str.contains('NB'))
            include_EB = (lrs['RTE_NM'].str.startswith('R-VA') & lrs['RTE_NM'].str.contains('EB'))
            include_PA = (lrs['RTE_NM'].str.startswith('R-VA') & lrs['RTE_NM'].str.contains('PA'))
            include_S_PR = (lrs['RTE_NM'].str.startswith('S-VA') & lrs['RTE_NM'].str.contains('PR'))

            if ramps == 0:
                lrs = lrs.loc[(include_interstates & exclude_ramps) | include_NB | include_EB | include_PA | include_S_PR]
            
            if ramps == 1:
                lrs = lrs.loc[(lrs['RTE_NM'].str.startswith('R-VA') & include_ramps)]

        lrs = lrs.to_crs(epsg=epsg)

        return lrs


    def load_m_values(self, lrsPath):
        """ Reads the LRS as a shapefile and returns an m-value dictionary 
        
        Args:
            lrsPath: Path to the LRS as a shapefile
            
        Returns:
            mValueDict: A python dictionary containing m-values
        """
        mValueDict = {}
        with shapefile.Reader(lrsPath) as shp:
                for row in shp.iterShapeRecords():
                    try:
                        record = row.record
                        shape = row.shape
                        rte_nm = record['RTE_NM']
                        mValues = shape.m
                        mValueDict[rte_nm] = mValues
                    except:
                        # print('Error finding m-value for row')
                        continue

        return mValueDict


def get_nearby_routes(point, distance, lrs_object):
    """ Returns a list of RTE_NMs within the given distance of the
        input point 
        
        Args:
            point: Shapely point object.  Projection must match
                lrs projection, default espg 3968.
            distance: Distance in lrs projection's units
            lrs_object: Instance of LRS
            
        Returns:
            routes: A list of nearby routes by RTE_NM
        """
    
    # Set the lrs and spatial index from the LRS object
    lrs = lrs_object.geodataframe
    lrs_index = lrs_object.spatial_index

    # Locate LRS routes located within the input distance of the input point
    pointBuffer = point.buffer(distance)
    possible_routes_index = list(lrs_index.query(pointBuffer))
    possible_routes = lrs.iloc[possible_routes_index]
    routes = possible_routes[possible_routes.intersects(pointBuffer)]["RTE_NM"].tolist()
    
    return routes


def calculate_azimuth(segment, normalized=False):
    """
    Calculates the azimuth between two points on Earth.

    Args:
        segment: Shapely LineString
        normalized: If true, will return a bearing between 0-180

    Returns:
        The azimuth in degrees (clockwise from North)
    """
    lon1, lat1 = segment.coords[0]
    lon2, lat2 = segment.coords[1]

    dlon = lon2 - lon1
    y = math.sin(math.radians(dlon)) * math.cos(math.radians(lat2))
    x = math.cos(math.radians(lat1)) * math.sin(math.radians(lat2)) - \
        math.sin(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.cos(math.radians(dlon))
    bearing = math.atan2(y, x)
    bearing = math.degrees(bearing)
    bearing = (bearing + 360) % 360

    if normalized:
        return bearing % 180

    return bearing


def get_nearby_route_by_segments(segment, distance, lrs_object, bearing_threshold=10):
    """ This function is similar to get_nearby_routes, but it takes into consideration
        the geometry of the LRS routes.  The input is a segment rather than a point
        and only routes whose nearby segment closely matches the bearing of the input
        segment will be returned.  This takes more time than get_nearby_routes, but
        it will result in more accurate matches.

        Args:
            segment: Shapely LineString object
            distance: Distance in lrs projection's units
            lrs_object: Instance of LRS
            bearing_threshold: The number of degrees that the bearing of the matched
                route should be in relationo to the bearing of the input segment

        Returns:
            routes: The RTE_NM for the closest route with a closely matching bearing
    """
    # Set the lrs and spatial index from the LRS object
    lrs = lrs_object.geodataframe
    lrs_index = lrs_object.spatial_index

    # Get input segment bearing
    segment_bearing = calculate_azimuth(segment, normalized=True)

    # Find nearest LRS routes
    segBuffer = segment.centroid.buffer(distance)
    possible_routes_index = list(lrs_index.query(segBuffer))
    possible_routes = lrs.iloc[possible_routes_index]
    nearby_routes = possible_routes[possible_routes.intersects(segBuffer)]
    print(f'len nearby_routes: {len(nearby_routes)}')
    # Locate nearest segment from each of the located routes


    def get_nearest_route_segment(record):
        print(record)
        route_geom = record.geometry
        ranked_segments = []
        for i in range(len(route_geom.coords) - 1):
            route_segment = LineString([route_geom.coords[i], route_geom.coords[i+1]])
            segment_distance = route_segment.distance(segment.centroid)
            ranked_segments.append((route_segment, segment_distance))
        ranked_segments.sort(key=lambda x: x[1])  # Sort by segment_distance
        route_segment_bearing = calculate_azimuth(ranked_segments[0][0], normalized=True)
        record['closest_segment'] = ranked_segments[0][0]
        record['segment_distance'] = ranked_segments[0][1]
        record['bearing'] = route_segment_bearing
        record['bearing_difference'] = abs(segment_bearing - route_segment_bearing)
        record['suitable_match'] = True if abs(segment_bearing - route_segment_bearing) < bearing_threshold else False
        print(f'get_nearest_route_segment record: {record}')
        return record

    
    nearby_routes = nearby_routes.apply(get_nearest_route_segment, axis=1)
    print(f'nearby_routes: {nearby_routes}')
    suitable_matches = nearby_routes.loc[nearby_routes['suitable_match'] == True]['RTE_NM'].tolist()
    
    return suitable_matches
    # Calculate bearing of nearest segment

    # Compare bearings to input segment



    # Rank segments based on distance and directionality
    # ranked_segments = []
    # for polyline in nearest_polylines.itertuples():
    #     for i in range(len(polyline.geometry.coords) - 1):
    #         segment = LineString([polyline.geometry.coords[i], polyline.geometry.coords[i+1]])
    #         distance = segment.distance(segment.centroid)
    #         # segment_bearing = azimuth(segment.coords[0], segment.coords[-1])
    #         # bearing_diff = abs(tmc_bearing - segment_bearing)
    #         # ranked_segments.append((segment, distance, bearing_diff))
    #         ranked_segments.append((segment, distance))

    # ranked_segments.sort(key=lambda x: x[1] + x[2])  # Sort by distance and bearing difference
    # ranked_segments.sort(key=lambda x: x[1])  # Sort by distance

    # return gp.GeoSeries([ranked_segments[0][0]])
    # return ranked_segments



if __name__ == '__main__':
    lrs_path = r'C:\Users\daniel.fourquet\Documents\Tasks\TMC Conflation 2025\Data\LRS_MASTER_RICHMOND.shp'
    lrs = LRS(lrs_path, filter=False)
    # lrs.geodataframe['geometry'] = shapely.line_merge(lrs.geodataframe.geometry)
    # lrs.geodataframe = lrs.geodataframe.explode()
    # lrs.geodataframe.reset_index(inplace=True)
    # lrs.geodataframe.head(1)[['RTE_NM', 'geometry']].to_clipboard()
    # print(len(lrs.geodataframe))
    print('ugghh', '\n' * 8)
    test_point = Point(185418.89, 74979.19)
    # # test_point = Point(185956.54383137805, 173669.74369517615)
    # # test_point = Point(185956.73, 173669.59)
    results = get_nearby_routes(test_point, 150, lrs)
    
    # test_segment = LineString([[180732.34, 75652.74], [180715.58, 75712.97]])  # 64WB
    # test_segment = LineString([[180416.56, 74055.36], [180402.43, 74070.38]])  # Park Ave
    # test_segment = LineString([[174946.94913422744, 172621.3780545953], [174956.873458169, 172616.04427231924]])  # Cleveland St
    # test_segment = LineString([[185418.89, 74979.19], [185423.62, 74996.99]])  # Hartman St
    # results = get_nearby_route_by_segments(test_segment, 150, lrs)

    print(results)