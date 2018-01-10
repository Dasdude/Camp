from geopy import distance
import numpy as np
def  gps_to_cartesian_dist_angle(lat1_ref,long1_ref,lat2,long2):
	""":returns {distance,angle}"""
	result = {}
	result['distance'] = distance.vincenty((lat1_ref,long1_ref),(lat2,long2)).meters
	# result['angle'] = 0 if result['distance']<.5 else np.sign(lat2-lat1_ref)*np.arccos(np.sign(long2-long1_ref)*distance.vincenty((lat1_ref,
	#                                                                                                                               long1_ref),
	#                                                                                                                   (lat1_ref,
	#                                                                                                                 long2)).meters/result[
	# 	'distance'])
	if result['distance'] < .5:
		result['angle'] = 0
	else:
		dx = distance.vincenty((lat1_ref, long1_ref), (lat1_ref, long2)).meters
		lat_sign = np.sign(lat2 - lat1_ref)
		long_sign= np.sign(long2 - long1_ref)
		if dx>result['distance']:
			result['angle'] = lat_sign * np.arccos(long_sign * 1)
		else:
			result['angle']=lat_sign* np.arccos(long_sign *dx / result['distance'])
	return result


def gps_to_cartesian_xy(lat1_ref, long1_ref, lat2, long2):
	""":returns {x,y} x being east/west, y north/south"""
	result = {}
	result['x'] = np.sign(long2-long1_ref)*distance.vincenty((lat1_ref, long1_ref), (lat1_ref, long2)).meters
	result['y'] = np.sign(lat2 - lat1_ref) * distance.vincenty((lat1_ref, long1_ref), (lat2, long1_ref)).meters
	return result
