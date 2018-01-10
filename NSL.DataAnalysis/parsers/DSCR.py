import numpy as np
import pandas as pd
import re
import os
import pandautil.csv_utils as csv_u
import matplotlib
import gmplot
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from utils.opts_class import Opts
from utils import gps_utils
import geopy
FILENAME_TEMPLATE = '63210002_D051_2017_03_26_09_33_05.csv'
a = 1
html_color_codes = {
	'darkblue'            : '#00008B',
	'darkcyan'            : '#008B8B',
	'darkgoldenrod'       : '#B8860B',
	'darkgray'            : '#A9A9A9',
	'darkgreen'           : '#006400',
	'darkkhaki'           : '#BDB76B',
	'darkmagenta'         : '#8B008B',
	'darkolivegreen'      : '#556B2F',
	'darkorange'          : '#FF8C00',
	'darkorchid'          : '#9932CC',
	'darkred'             : '#8B0000',
	'darksalmon'          : '#E9967A',
	'darkseagreen'        : '#8FBC8F',
	'darkslateblue'       : '#483D8B',
	'darkslategray'       : '#2F4F4F',
	'darkturquoise'       : '#00CED1',

}
SEALBEACH_REF = {'lat_ref': 33.774000,'long_ref': -118.069827}#as longitude increases we move toward east,latitude increase
SEALBEACH_COORD = {'long_min':-118.069363,'long_max':-118.051763}
def get_relative_xy_list(lat_list,long_list):
	# lat is y and long is x
	x_list = []
	y_list = []
	for i,_ in enumerate(lat_list):
		res = gps_utils.gps_to_cartesian_xy(SEALBEACH_REF['lat_ref'], SEALBEACH_REF['long_ref'], lat_list[i], long_list[i])
		x_list+=[res['x']]
		y_list += [res['y']]
	return {'x':x_list,'y':y_list}


def get_relative_xy(lat, long):
	# lat is y and long is x
	res = gps_utils.gps_to_cartesian_xy(SEALBEACH_REF['lat_ref'], SEALBEACH_REF['long_ref'], lat, long)
	return {'x': res['x'], 'y': res['y']}


def get_polar_list_input(lar_lor_lat_lot):
	res = gps_utils.gps_to_cartesian_dist_angle(lar_lor_lat_lot[0],lar_lor_lat_lot[1],lar_lor_lat_lot[2],lar_lor_lat_lot[3])
	return pd.Series(res)
def get_polar(lat_reciever,long_reciever,lat_transmitter,long_transmitter):
	res = gps_utils.gps_to_cartesian_dist_angle(lat_reciever
	                                    , long_reciever, lat_transmitter, long_transmitter)
	return pd.Series({'distance':res['distance'],'angle':res['angle']})
def animate_location(df_list,lat_range,long_range,opts):
	def data_gen(t=0):
		cnt = 0
		update_step = 100
		df_merged = pd.DataFrame(pd.concat(df_list))
		df_ste = df_merged[df_merged.LogRecType=='STE']
		df_ste = df_ste.sort_values('TimeStamp_ms',ascending=True)
		df_ste['UniqueOBE_ID_Alias'] = df_ste['UniqueOBE_ID_Alias'].apply(lambda x: int(int(re.sub("[^0-9]","",x))-51))
		df_ste= df_ste[df_ste.UniqueOBE_ID_Alias<5]
		while cnt < df_merged.last_valid_index():
			lat = df_ste.lat.values[cnt:cnt+update_step]
			long = df_ste.long.values[cnt:cnt+update_step]
			id = df_ste.UniqueOBE_ID_Alias.values[cnt:cnt+update_step]
			time = df_ste.TimeStamp_ms.values[cnt:cnt+update_step]
			cnt += update_step
			xy = get_relative_xy_list(lat, long)
			yield {'lat':xy['y'],'long':xy['x'],'id':id,'time':time}

	def init():
		xy_max = get_relative_xy_list([lat_range['max']],[long_range['max']])
		ax.set_xlim(0, xy_max['x'][0])
		ax.set_ylim(0, xy_max['y'][0])
		for xdata_sample in xdata:
			del xdata_sample
		for ydata_sample in ydata:
			del ydata_sample
		for idx,line_sample in enumerate(line):
			line_sample.set_data(xdata[idx],ydata[idx])
		# line.set_data(xdata, ydata)
		return tuple(line)

	fig, ax = plt.subplots()
	fig.axes[0].set_xlabel('long')
	fig.axes[0].set_ylabel('lat')

	line=[]
	xdata,ydata=[],[]
	for df in df_list:
		line_temp, = ax.plot([],[],lw=2)
		line+=[line_temp]
		xdata+=[[]]
		ydata+=[[]]

	ax.grid()

	def run(data):
		# update the data
		lat_list, long_list,id =data['lat'],data['long'],data['id']
		# print(len(xdata),data['id'])
		for idx,did in enumerate(id):
			ydata[did].append(lat_list[idx])
			xdata[did].append(long_list[idx])
			line[did].set_data(xdata[did],ydata[did])

		return tuple(line)

	# with writer.saving(fig,'./hello.mp4',5):
	ani = animation.FuncAnimation(fig, run, data_gen, blit=True, interval=0, repeat=False, init_func=init)
	plt.show()
def dsrc_crop_range(abs_path):
	path = abs_path

	df_list = []
	max_car = 10
	opts = Opts()
	opts.max_car = max_car
	file_names = [];
	file_names += [each for each in os.listdir(abs_path) if each.endswith('.csv')]
	for file_name in file_names:
		df = pd.read_csv(os.path.join(path, file_name))
		df = df.dropna(axis=1, how='all')
		if 'long' in df.columns:
			df = df[df.long >= SEALBEACH_COORD['long_min']]
			df = df[df.long < SEALBEACH_COORD['long_max']]
			df.to_csv(os.path.join(path, file_name))
		else:
			print('no long attr')
def get_DSRC_dataframes(path_to_experiment,labels=['lat', 'rxlat', 'rxlong', 'LogRecType', 'Range', 'RSS', 'UniqueOBE_ID_Alias', 'long', 'speed',
                                                   'heading', 'TimeStamp_ms', 'MsgSeqNum']):
	df_list = 10*[None]
	file_names = [];
	file_names += [each for each in os.listdir(path_to_experiment) if each.endswith('.csv')]
	for file_name in file_names:
		print(os.path.join(path_to_experiment,file_name))
		df = pd.read_csv(os.path.join(path_to_experiment, file_name))
		df = df[df.LogRecType=='RXE']
		df = df[labels]
		name = re.findall("D0[0-9]*", file_name)
		cid = int(name[0][1:])-51
		df['egoID'] = cid
		df = df.dropna(axis=1, how='all')
		col_remove_list =[]
		for col_str in df.columns:
			if col_str.find('Unnamed') >=0:
				col_remove_list+=[col_str]
		df=df.drop(col_remove_list,axis=1)
		df_list[cid]= df
	return df_list
def animate_gps(abs_path):
	path =abs_path

	df_list = []
	df_lat_range_list = {'max': [], 'min': []}
	df_long_range_list = {'max': [], 'min': []}
	max_car = 5
	opts = Opts()
	opts.max_car = max_car
	file_names = [];file_names += [each for each in os.listdir(abs_path) if each.endswith('.csv')]
	for file_name in file_names:
		df = pd.read_csv(os.path.join(path,file_name))
		df = df.dropna(axis=1, how='all')
		df_list += [df]
		df_lat_range_list['max'] += [df.lat.max()]
		df_lat_range_list['min'] += [df.lat.min()]
		df_long_range_list['max'] += [df.long.max()]
		df_long_range_list['min'] += [df.long.min()]
	df_lat_range = {'min': np.min(df_lat_range_list['min']), 'max': np.max(df_lat_range_list['max'])}
	df_long_range = {'min': np.min(df_long_range_list['min']), 'max': np.max(df_long_range_list['max'])}
	# FFMpegWriter = animation.writers['ffmpeg']
	metadata = dict(title='DSRC sync car path visualization', artist='Ehsan Emad', comment='v1')
	# writer = FFMpegWriter(fps=2, metadata=metadata)
	coordinates_min = str(df_lat_range['min'])+','+str(df_long_range['min'] )
	coordinates_max = str(df_lat_range['max']) + ',' + str(df_long_range['max'])
	print(coordinates_min)
	print(coordinates_max)
	animate_location(df_list, df_lat_range, df_long_range, opts)
def plot_scenario_gps_html(abs_path,experiment_tag,save_path):
	path = abs_path

	df_list = []
	df_lat_range_list = {'max': [], 'min': []}
	df_long_range_list = {'max': [], 'min': []}
	max_car = 10
	opts = Opts()
	opts.max_car = max_car
	file_names = [];
	file_names += [each for each in os.listdir(abs_path) if each.endswith('.csv')]
	for file_name in file_names:
		df = pd.read_csv(os.path.join(path, file_name))
		df = df.dropna(axis=1, how='all')
		# df[['long', 'lat']] = df[['long', 'lat']].apply(lambda x: (x))
		df_list += [df]
		df_lat_range_list['max'] += [df.lat.max()]
		df_lat_range_list['min'] += [df.lat.min()]
		df_long_range_list['max'] += [df.long.max()]
		df_long_range_list['min'] += [df.long.min()]
	df_lat_range = {'min': np.min(df_lat_range_list['min']), 'max': np.max(df_lat_range_list['max'])}
	df_long_range = {'min': np.min(df_long_range_list['min']), 'max': np.max(df_long_range_list['max'])}
	# FFMpegWriter = animation.writers['ffmpeg']
	metadata = dict(title='DSRC sync car path visualization', artist='Ehsan Emad', comment='v1')
	# writer = FFMpegWriter(fps=2, metadata=metadata)
	coordinates_min = str(df_lat_range['min']) + ',' + str(df_long_range['min'])
	coordinates_max = str(df_lat_range['max']) + ',' + str(df_long_range['max'])
	gmap = gmplot.GoogleMapPlotter((df_lat_range['min'] + df_lat_range['max']) / 2, (df_long_range['min'] + df_long_range['max']) / 2, zoom=15)
	# gmap.scatter(df_list[0].lat.values,df_list[0].long.values,color='r',size=100)
	colors = list(html_color_codes)
	for idx,df in enumerate(df_list):

		df_ste = df[df.LogRecType == 'STE']
		gmap.plot(df_ste.lat.values, df_ste.long.values, threshold=1, radius=25,color=colors[idx])
	# gmap.paths(df_list[0].lat.values, df_list[0].long.values)
	gmap.draw(os.path.join(save_path, experiment_tag+'.html'))
def pythogorean_relative_coordinate(lat1_ref,long1_ref,lat2,long2):
	y= np.sin(lat2-lat1_ref)
	x = np.sin(long2-long1_ref)
	return {'x':x,'y':y}
def merge_all_DSRC(experiment_path,max_experiment=-1):
	df_list=[]
	cnt = max_experiment
	for day in os.listdir(experiment_path):
		if not 'March' in day:
			break
		for experiment_num_1 in os.listdir(os.path.join(experiment_path, day)):
			for experiment_num_2 in os.listdir(os.path.join(experiment_path, day, experiment_num_1)):
				if(cnt==0):break
				path = os.path.join(experiment_path, day, experiment_num_1, experiment_num_2)
				file_names = [];
				file_names += [each for each in os.listdir(path) if each.endswith('.csv')]
				for file_name in file_names:
					df = pd.read_csv(os.path.join(path, file_name))
					# df = df.dropna(axis=1, how='all')
					df_list+=[df]
				cnt-=1
				print(path)
	res_df = pd.concat(df_list)
	res_df = res_df.sort_values('TimeStamp_ms', ascending=False)
	return pd.DataFrame(res_df).to_csv(os.path.join(experiment_path,'all.csv'))

if __name__ == '__main__':
	date = 'March 26'
	experiment_path = '../Dataset/DSRC'
	jump=5
	counter = 4
	# merge_all_DSRC(experiment_path,1)
	# animate_gps(os.path.abspath('../Dataset/DSRC/{}/1/2/'.format(date)))
	for day in os.listdir(experiment_path):
		if not 'March' in day:
			break
		for experiment_num_1 in os.listdir(os.path.join(experiment_path,day)):
			for experiment_num_2 in os.listdir(os.path.join(experiment_path, day,experiment_num_1)):
				if jump>0:
					jump-=1
					continue
				if counter==0:break
				# try:
				path = os.path.join(experiment_path, day, experiment_num_1,experiment_num_2)
				print(path)
				animate_gps(path)
				counter-=1
				# dsrc_crop_range(path)
				# plot_scenario_gps_html(path,'{}_{}_{}'.format(day,experiment_num_1,experiment_num_2),save_path=os.path.abspath(
				# os.path.join(experiment_path,'DSRC_GPS_PLOT')))
				# except AttributeError:
				# 	print('f')
