import geopandas as gp
import pandas as pd
import shapefile
from shapely import intersects
from shapely.geometry import Point, LineString, MultiLineString
import shapely
import math
from pyproj import Transformer
transformer = Transformer.from_crs("EPSG:3968", "EPSG:4326", always_xy=True)

class LRS:
    def __init__(self, lrs_shp_path, filter=False, ramps=2, epsg=3968):
        # Initialize filtered geodataframe
        self.geodataframe = self.open_lrs(lrs_shp_path, filter, ramps, epsg)

        # Initialize unfiltered geodataframe
        if filter:
            self.geodataframe_unfiltered = self.open_lrs(lrs_shp_path)
        else:
            self.geodataframe_unfiltered = self.geodataframe

        # List of M-values
        self.mDict = self.load_m_values(lrs_shp_path)

        # Dictionary of opposite direction route names 
        self.opp_route_dict = self.load_opp_route_dict()

        # Reference to spatial index
        self.spatial_index = self.geodataframe.sindex
        self.spatial_index_unfiltered = self.geodataframe_unfiltered.sindex


    def open_lrs(self, lrs_shp_path, filter=False, ramps=None, epsg=3968):
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
                        mValues = shape.m
                        mValueDict[record.ID] = mValues
                    except Exception as e:
                        print(e)
                        continue

        return mValueDict
    
    def load_opp_route_dict(self):
        opp_route_df = self.geodataframe[['RTE_NM', 'RTE_OPPOSI']].copy()
        opp_route_df.drop_duplicates()
        opp_route_df.set_index('RTE_NM', inplace=True)
        return opp_route_df.to_dict()

    
    def get_opp_route(self, rte_nm):
        return self.opp_route_dict['RTE_OPPOSI'].get(rte_nm)


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

    # Convert to degrees
    begin_point = transformer.transform(segment.coords[0][0], segment.coords[0][1])
    end_point = transformer.transform(segment.coords[1][0], segment.coords[1][1])
    segment = LineString((begin_point, end_point))

    x1, y1 = segment.coords[0]
    x2, y2 = segment.coords[1]

    def get_bearing(lon1, lat1, lon2, lat2):
        dlon = lon2 - lon1
        y = math.sin(math.radians(dlon)) * math.cos(math.radians(lat2))
        x = math.cos(math.radians(lat1)) * math.sin(math.radians(lat2)) - \
            math.sin(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.cos(math.radians(dlon))
        bearing = math.atan2(y, x)
        bearing = math.degrees(bearing)
        bearing = (bearing + 360) % 360

        return bearing

    bearing = get_bearing(x1, y1, x2, y2)
    if normalized:
        bearing_reversed = get_bearing(x2, y2, x1, y1)
        bearing = min(bearing, bearing_reversed)

    return bearing


def get_nearby_route_by_segments(segment, distance, lrs_object, bearing_threshold=15, closest_match=True, expanded_buffer=False, unfiltered=False):
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
                route should be in relation to the bearing of the input segment
            closest_match: If True, returns only the RTE_NM that is the closest
                to the input segment.  If False, returns all potential RTE_NM matches
                within the search distance and bearing threshold.
            unfiltered: If True, then the unfiltered version of the LRS is used to
                match routes by segment.  This is useful on divided highways where
                the segments are a significant distance (eg 50m) from each other
                and they would be missed if non-prime direction is filtered out.

        Returns:
            routes: The RTE_NM for the closest route with a closely matching bearing
    """
    # Set the lrs and spatial index from the LRS object
    lrs = lrs_object.geodataframe
    lrs_index = lrs_object.spatial_index
    
    if unfiltered:
        lrs = lrs_object.geodataframe_unfiltered
        lrs_index = lrs_object.spatial_index_unfiltered


    # Get input segment bearing
    segment_bearing = calculate_azimuth(segment, normalized=True)

    # Find nearest LRS routes
    segBuffer = segment.centroid.buffer(distance)
    possible_routes_index = list(lrs_index.query(segBuffer))
    possible_routes = lrs.iloc[possible_routes_index]
    nearby_routes = possible_routes[possible_routes.intersects(segBuffer)]
    
    
    # Locate nearest segment from each of the located routes
    def get_nearest_route_segment(record):
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
        record['input_segment'] = segment
        record['segment_bearing'] = segment_bearing
        record['bearing_raw'] = calculate_azimuth(ranked_segments[0][0], normalized=False)
        record['bearing'] = route_segment_bearing
        record['bearing_difference'] = abs(segment_bearing - route_segment_bearing)
        record['suitable_match'] = True if abs(segment_bearing - route_segment_bearing) < bearing_threshold else False
        return record

    
    nearby_routes = nearby_routes.apply(get_nearest_route_segment, axis=1)
    if len(nearby_routes) == 0:  # If no nearby routes are found...
        if expanded_buffer == False:  # Try expanding the buffer by 5m
            return get_nearby_route_by_segments(segment, distance+5, lrs_object, expanded_buffer=True)
        if expanded_buffer == True and unfiltered == False:  # If still no match, try using unfiltered LRS
            return get_nearby_route_by_segments(segment, distance+5, lrs_object, expanded_buffer=True, unfiltered=True)
        else:  # If still no match, then return empty result
            return []


    suitable_matches = nearby_routes.loc[nearby_routes['suitable_match'] == True]

    # Return only the closest 
    minimum_distance = suitable_matches['segment_distance'].min()
    suitable_matches_rte_nm = suitable_matches.loc[suitable_matches['segment_distance'] == minimum_distance]['RTE_NM'].tolist()

    return suitable_matches_rte_nm


# def get_nearby_route_and_m_by_segment(segment, distance, lrs_object, bearing_threshold=15, closest_match=True, expanded_buffer=False, unfiltered=False):
#     """ This function is similar to get_nearby_route_by_segments, but it returns
#         a record with RTE_NM and m-values for the nearest segment

#         Args:
#             segment: Shapely LineString object
#             distance: Distance in lrs projection's units
#             lrs_object: Instance of LRS
#             bearing_threshold: The number of degrees that the bearing of the matched
#                 route should be in relation to the bearing of the input segment
#             closest_match: If True, returns only the RTE_NM that is the closest
#                 to the input segment.  If False, returns all potential RTE_NM matches
#                 within the search distance and bearing threshold.
#             unfiltered: If True, then the unfiltered version of the LRS is used to
#                 match routes by segment.  This is useful on divided highways where
#                 the segments are a significant distance (eg 50m) from each other
#                 and they would be missed if non-prime direction is filtered out.

#         Returns:
#             routes: The RTE_NM for the closest route with a closely matching bearing
#     """
#         # Set the lrs and spatial index from the LRS object
#     lrs = lrs_object.geodataframe
#     lrs_index = lrs_object.spatial_index
    
#     if unfiltered:
#         lrs = lrs_object.geodataframe_unfiltered
#         lrs_index = lrs_object.spatial_index_unfiltered



def get_m_value(point, lrs_object, rte_nm):
    """ Returns the m-value of the nearest location along the rte_nm on the input
        lrs to the input point.

        Args:
            point: Shapely Point object
            lrs_object: Instance of LRS
            rte_nm: The RTE_NM in the LRS that will be used to find the m-value

        Returns:
            msr: The m-value along the input RTE_NM closest to the input point
    """

    try:
        lrs = lrs_object.geodataframe_unfiltered
        
        # Get LineString of the closest segment of {rte_nm} to the input point
        lrs_route = lrs.loc[lrs['RTE_NM'] == rte_nm].copy()

        if len(lrs_route) == 0:
            return None

        if len(lrs_route) > 1:
            # Find distance of each segment to the input point
            # Return segment with minimum distance
            lrs_route['dist'] = lrs_route.geometry.distance(point)
            lrs_route = lrs_route.loc[lrs_route['dist'] == lrs_route['dist'].min()]
    
        lrs_route = lrs_route.iloc[0]
        route_geom = lrs_route.geometry  # Geometry of lrs route
        route_id = lrs_route['ID']  # ID used to look up m-values
        mList = lrs_object.mDict[route_id]  # List of m-values
        
        # Get point along route_geom that is closest to the input point
        point_on_geom = route_geom.interpolate(route_geom.project(point))

        # Find route segment that intersects point_on_geom and obtain m-values
        # of begin and end point of segment
        for i in range(len(route_geom.coords)):
            if i+2 > len(route_geom.coords):
                break

            p1 = route_geom.coords[i]
            p2 = route_geom.coords[i+1]
            line = LineString((p1,p2))
            if intersects(line, point_on_geom.buffer(1)):
                m1 = mList[i]
                m2 = mList[i+1]
                p3 = line.project(point_on_geom, normalized=True)
                break

        # Interpolate segment m-values to determine input point measure
        msr = min(m1, m2) + (abs(m2 - m1) * p3)
        
        return round(msr, 3)
    except Exception as e:
        return -1


def get_segment_msr(segment, lrs, rte_nm):
    """ Returns the begin and end measure of the input segment on the lrs.

        Args:
            segment: Shapely LineString object
            lrs_object: Instance of LRS
            rte_nm: The RTE_NM in the LRS that will be used to find the m-values

        Returns:
            begin_msr: The m-value of the first point along the input segment
            end_msr:  The m-value of the last point along the input segment
    """
    if isinstance(segment, MultiLineString):
        begin_point = shapely.get_point(segment.geoms[0], 0)
        end_point = shapely.get_point(segment.geoms[-1], -1)
    else:
        begin_point = shapely.get_point(segment, 0)
        end_point = shapely.get_point(segment, -1)
        
    begin_msr = get_m_value(begin_point, lrs, rte_nm)
    end_msr = get_m_value(end_point, lrs, rte_nm)

    return begin_msr, end_msr


def locate_features_along_routes(geodataframe, record_id, lrs, in_rte_nm, out_rte_nm="RTE_NM", begin_msr="Begin_Msr", end_msr="End_Msr"):
    """ Finds the begin and end m-values for the input geodataframe.  The route name must
        be identified first.

        Args:
            geodataframe: GeoDataFrame with at least a geometry, segment id, and route name
            record_id: The name of the input column containing the ID for the segment,
                eg TMC code
            lrs: Instance of LRS
            in_rte_nm: The name of the input column containing the LRS RTE_NM value that will be
                used to find the m-values
            out_rte_nm: The name of the output column containing the LRS RTE_NM
            begin_msr: The name of the output column containing the begin measure
            end_msr: The name of the output column containing the end measure

        Returns:
            New DataFrame containing only the following columns:
                record_id, out_rte_nm, begin_msr, end_msr
    """
    def get_msrs(record):
        record[begin_msr], record[end_msr] = get_segment_msr(record.geometry, lrs, record[in_rte_nm])
        return record

    geodataframe = geodataframe.apply(get_msrs, axis=1)
    geodataframe.rename(columns={in_rte_nm: out_rte_nm}, inplace=True)
    
    return geodataframe


def flip_routes(geodataframe, rte_nm, begin_msr, end_msr, lrs, ignore_interstates=True):
    """ Locates segments where the begin_msr is greater than the end_msr.
        These are likely non-prime segments that should be flipped to the
        opposite direction of the route.  In these cases, the rte_nm is
        reassigned and the begin and end measures are recalculated.

        Args:
            geodataframe: GeoDataFrame with at least a rte_nm, begin_msr, and end_msr
            rte_nm: The RTE_NM in the LRS that will be used to find the opposite direction
                route
            begin_msr: The name of the column containing the begin measure
            end_msr: The name of the column containing the end measure
            lrs: Instance of LRS
            ignore_interstates: If True, interstate routes will return unflipped

        Returns:
            New DataFrame containing flipped values with new begin and end measures
    """

    def get_opposite_route(record):
        if record[rte_nm] is None:
            return record
            
        if ignore_interstates:
            if record[rte_nm].startswith('R-VA   IS'):
                return record

        if record[begin_msr] < record[end_msr]:
            # Located on correct side of route
            # Return existing record
            return record

        # Locate opposite route from lrs object's dictionary
        opposite_route = lrs.get_opp_route(record[rte_nm])

        # If no opposite route to flip to, return existing record
        if opposite_route is None:
            return record

        record[rte_nm] = opposite_route

        # Recalculate measures
        record[begin_msr], record[end_msr] = get_segment_msr(record.geometry, lrs, record[rte_nm])

        return record
    
    geodataframe = geodataframe.apply(get_opposite_route, axis=1)

    return geodataframe


        





if __name__ == '__main__':
    lrs_path = r'C:\Users\daniel.fourquet\Documents\Tasks\TMC Conflation 2025\Data\LRS_MASTER_RICHMOND_SINGLEPART.shp'
    lrs = LRS(lrs_path, filter=True, ramps=0)

    tmcs_path = r'C:\Users\daniel.fourquet\Documents\Tasks\TMC Conflation 2025\NPMRDS\data\NPMRDS_medium_test.shp'
    tmcs = gp.read_file(tmcs_path)
    tmcs = tmcs[['Tmc', 'RoadNumber', 'RoadName', 'TmcLinear', 'geometry']]
    tmcs = tmcs.to_crs(epsg=3968)

    tmcs = tmcs.loc[tmcs['Tmc'] == '110N04907']
    tmcs['matching_rte_nm'] = 'R-VA   IS00064EB'
    result = locate_features_along_routes(tmcs, 'Tmc', lrs, 'matching_rte_nm')
    result.to_clipboard(index=False)